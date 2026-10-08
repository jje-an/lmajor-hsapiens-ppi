from Bio import SeqIO
import matplotlib.pyplot as plt
import torch
import os
import sqlite3
import json
import re
import warnings

"""
Helper functions for ESMFold2 Snakemake workflow
"""


def fasta_to_hash(file_handle, key):
    """
    Add each record in a fasta library into a dictionary with keys as defined
    in config.yaml. Trey, this should be really similar to parts 
    of Bio::Adventure::Structure

    Arguments:
        file_handle - File handle of the opened fasta library.
        key - Type of the key in each key-value pair (ex. transcript id, gene id).
    """

    output_dict = {}

    for record in SeqIO.parse(file_handle, "fasta"):

        # Weird things to accomodate for ensembl library weirdness.
        # Stole this from Trey
        if key:
            seqid = record.description
        else:
            seqid = record.id
        short_id = seqid
        if key:
            match = re.search(rf"{key}\s*[:|=](\S+)", seqid)
        else:
           match = None
        if (match):
            short_id = match.group(1) 
            ## Get rid of terminal ID numbers (assume it doesn't pass 99)
            short_id = re.sub(r"\.\d{1,2}$", "", short_id) 
        output_dict[short_id] = record

    return output_dict


def complex_get_metrics(result, job_id, database_path, output_dir, keep_distogram=False):
    """
    Retrieves metrics from the ESMFold2 result and uploads metrics to a SQLite3
    database. If uploading to database fails, stores metrics in a JSON file instead.

    Arguments:
        result - MolecularComplexResult object from output of fold(). 
        job_id - {hsapiens_transcript_id}_{lmajor_transcript_id} formated string.
        database_path - File path of database to write to.
        keep_distogram - If true, will save distograms to the database. By default,
            distograms will not be saved, and NULL will be uploaded to the column element.
            Warning: distograms are NxNx64 size arrays, and will take up a significant
            amount of space.
    """

    con = sqlite3.connect(database_path)
    cur = con.cursor()

    # round each metric to 2 decimals and convert tensors to lists
    # convert to float64 for clean rounding
    ptm = round(result.ptm, 2)
    iptm = round(result.iptm, 2)
    plddt_mean = torch.round(result.plddt.mean().to(torch.float64), decimals=2).tolist()
    plddt = torch.round(result.plddt.to(torch.float64), decimals=2).tolist()
    pae = torch.round(result.pae.to(torch.float64), decimals=2).tolist()
    pair_chains_iptm = torch.round(result.pair_chains_iptm.to(torch.float64), decimals=2).tolist()
    residue_index = result.residue_index.tolist()
    entity_id = result.entity_id.tolist()
    gpu_type = torch.cuda.get_device_name()

    if keep_distogram:
        distogram = torch.round(result.distogram.to(torch.float64), decimals=2).tolist()
    else:
        distogram = None # will be inserted to the database as NULL object
        
    row_entry = (
        job_id,
        ptm,
        iptm,
        plddt_mean,
        str(plddt),
        str(pae),
        str(distogram) if distogram else None,
        str(pair_chains_iptm),
        str(residue_index),
        str(entity_id),
        gpu_type
    )
    
    cur.execute("BEGIN TRANSACTION;")
    try: 
        cur.execute(f"INSERT INTO esm_summaries VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);",
                    row_entry)
    except Exception as error:
        cur.execute("ROLLBACK;")
        con.close()
        print(error)

        # save to json file instead

        datum = {
            "job_id": job_id,
            "ptm": ptm,
            "iptm": iptm,
            "plddt_mean": plddt_mean,
            "plddt": plddt,
            "pae": pae,
            "distogram": distogram,
            "pair_chains_iptm": pair_chains_iptm,
            "residue_index": residue_index,
            "entity_id": entity_id,
            "gpu_type": gpu_type
        }

        json_str = json.dumps(datum, sort_keys=True, indent=1,separators=(',', ': '))
        json_str = json_str.replace('NaN', 'null')

        # remove newlines from each list entry
        def repl_func(match: re.Match):
            return " ".join(match.group().split())
        json_str = re.sub(r"(?<=\[)[^\[\]]+(?=])", repl_func, json_str)
    
        with open(f"{output_dir}/{job_id}_metrics.json", "w") as f:
            f.write(json_str)
        warnings.simplefilter("always")
        warnings.warn(f"Metrics could not be written to SQL database due to: {error}. \nMetrics have been saved to {output_dir}/{job_id}_metrics.json instead.", RuntimeWarning)
    else:
        cur.execute("COMMIT;")
        con.close()


def plot_pae_matrix(result, job_id, output_dir):
    """
    Retrieves PAE from the ESMFold2 result and saves a plot of the matrix.

    Arguments:
        result - MolecularComplexResult object from output of fold(). 
        job_id - {hsapiens_transcript_id}_{lmajor_transcript_id} formated string.
        database_path - File path of database to write to.
        output_dir - Output directory to save image file to.
    """

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pae = result.pae.numpy()

    plt.imshow(pae, cmap="summer")
    plt.title(f'PAE Matrix for {job_id}')
    plt.xlabel("Scored Residue", labelpad=20)
    plt.ylabel("Aligned Residue", labelpad=20)
    plt.colorbar(label="PAE (Angstroms)")

    #separate graph based on chain
    token_chain_ids = result.entity_id.tolist()
    unique_token_chain_ids = sorted(set(token_chain_ids))
    tokens_per_chain = []
    for i in range(len(unique_token_chain_ids)):
        tokens_per_chain.append(token_chain_ids.count(unique_token_chain_ids[i]))

    for i in range(1, len(unique_token_chain_ids)):
        
        # .5 centers line on boundary of token
        plt.axvline(sum(tokens_per_chain[0:i]) - .5, color='black', linewidth=.8)
        plt.axhline(sum(tokens_per_chain[0:i]) - .5, color='black', linewidth=.8)

    # boundary index where each chain's tokens start/end
    chain_starts = [0]
    for count in tokens_per_chain[:-1]:
        chain_starts.append(chain_starts[-1] + count)
    chain_ends = [start + count for start, count in zip(chain_starts, tokens_per_chain)]

    plt.tick_params(labelbottom=False, labelleft=False) 

    # chain name labels
    for start, end, chain_name in zip(chain_starts, chain_ends, unique_token_chain_ids):
        midpoint = (start + end) / 2 - 0.5
        plt.text(midpoint, len(token_chain_ids) + 2, f"Chain {chain_name}",
                ha='center', va='top', fontsize=10)
        plt.text(-2, midpoint, f"Chain {chain_name}",
                ha='right', va='center', fontsize=10, rotation=90)


    plt.savefig(f"{output_dir}/{job_id}_pae_matrix.png")
    plt.close()




