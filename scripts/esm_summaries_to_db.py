import glob
import json
import sys
import sqlite3
import warnings
import os
from pathlib import Path



def esm_summaries_to_db(input_dir, db_dir):
    """
    Searches for json files in the input directory from esm jobs that
    failed to upload metrics to the database. Adds metrics from each
    json file to database from db_dir. If the metrics were successfully
    added or the job_id already exists in the database, the json file
    is moved to input_dir/../old_json_files.

    Input directory should typically be results/pairs/ when using Jean's 
    snakemake workflow.

    Example usage:
    `python esm_summaries_to_db.py ../scavenger_partition/results/pairs ../lmajor_hsapiens_ppi.db`

    Arguments:
        input_dir - file path of parent directory for all jobs
        db_dir - file path of directory to add to.
    """

    
    # two directories deep to match file structure of snakemake outputs.
    files = sorted(glob.glob(f"{input_dir}/*/*/*.json"))
    print(f"Found {len(files)} .json files")
    if len(files) == 0:
        return

    os.makedirs(f"{input_dir}/../old_json_files", exist_ok=True)
    
    con = sqlite3.connect(db_dir)
    cur = con.cursor()

    cur.execute(
    "CREATE TABLE IF NOT EXISTS esm_summaries" \
        "(job_id TEXT PRIMARY KEY," \
        "ptm REAL," \
        "iptm REAL," \
        "plddt_mean REAL," \
        "plddt TEXT," \
        "pae TEXT," \
        "distogram TEXT," \
        "pair_chain_iptm TEXT," \
        "residue_index TEXT," \
        "entity_id TEXT," \
        "gpu_type TEXT);")
    
    written_count = 0
    for path in files:
        #print(f"reading {path}...") #really loud, uncomment if you want
        with open(path) as f:
            data = json.load(f)
        

        row = (
            data.get("job_id"),
            data.get("ptm"),
            data.get("iptm"),
            data.get("plddt_mean"),
            str(data.get("plddt")) if data.get("plddt") is not None else None,
            str(data.get("pae")) if data.get("pae") is not None else None,
            str(data.get("distogram")) if data.get("distogram") is not None else None,
            str(data.get("pair_chains_iptm")) if data.get("pair_chains_iptm") is not None else None,
            str(data.get("residue_index")) if data.get("residue_index") is not None else None,
            str(data.get("entity_id")) if data.get("entity_id") is not None else None,
            str(data.get("gpu_type")) if data.get("gpu_type") is not None else None
        )
        cur.execute("BEGIN TRANSACTION;")
        try: 
            cur.execute(f"INSERT INTO esm_summaries VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);",
                        row)
            written_count += 1

            # move json file to results/old_json_files
            basename = os.path.basename(path)
            new_path = f"{input_dir}/../old_json_files/{basename}"
            os.replace(path, new_path)
            print(f"{os.path.basename(path)} moved to {Path(new_path)}")

        except Exception as error:
            # return to previous state if error occurred during insertion
            cur.execute("ROLLBACK;")
            warnings.simplefilter("always")

            if str(error) == "UNIQUE constraint failed: esm_summaries.job_id":
                print(f"Metrics from {data.get("job_id")} weren't uploaded because the job_id already exists in the database.")

                # move json file to results/old_json_files
                basename = os.path.basename(path)
                new_path = f"{input_dir}/../old_json_files/{basename}"
                os.replace(path, new_path)
                print(f"{os.path.basename(path)} moved to {Path(new_path)}")

            else:
                warnings.warn(f"Metrics from {data.get("job_id")} could not be written to SQL database due to: {error}.", RuntimeWarning)
        else:
            cur.execute("COMMIT;")

    con.close()
    print(f"Wrote {written_count} rows to {db_dir}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <input_dir> <db_path>")
    esm_summaries_to_db(sys.argv[1], sys.argv[2])
