import sqlite3
import glob
import json
import warnings
import sys
from pathlib import Path

def filter_output_list(beacon_result_dir, scavenger_result_dir, db_dir):
    """
    Takes the result directory paths for both workflows and
    finds job_ids of jobs that completed (created .cif file) but
    did not write to the database. If these job_ids have 
    a json file, metrics will be uploaded to the database. 
    Otherwise, files will be deleted (after user permission) 
    so that the workflow will run those jobs and generate 
    metrics again.


    Example usage:
    'python scripts/filter_output_list.py beacon_partition/results/pairs 
        scavenger_partition/results/pairs lmajor_hsapiens_ppi.db'

    Arguments:
        beacon_result_dir - file path of directory containing results from the
            beacon partition workflow. For Jean's workflow it should be:
            `beacon_partition/results/pairs`

        scavenger_result_dir - file path of directory containing results from the
            scavenger partition workflow. For Jean's workflow it should be:
            `scavenger_partition/results/pairs`

        db_dir - file path of directory to check against. For Jeans workflow
            it is currently `lmajor_hsapiens_ppi.db`
    """

    beacon_files = sorted(glob.glob(f"{beacon_result_dir}/*/*/*/*.cif"))
    scavenger_files = sorted(glob.glob(f"{scavenger_result_dir}/*/*/*/*.cif"))
    print(f"Found {len(beacon_files) + len(scavenger_files)} .cif files")
    if len(beacon_files) + len(scavenger_files) == 0:
        return


    beacon_job_ids = set()
    beacon_lmajor_ids = set()
    scavenger_job_ids = set()
    scavenger_lmajor_ids = set()

    for file in beacon_files:
        beacon_job_ids.add(file.split("/")[-1].split(".cif")[0])
        beacon_lmajor_ids.add(file.split("/")[-1].split(".cif")[0].split("_")[1])
    for file in scavenger_files:
        scavenger_job_ids.add(file.split("/")[-1].split(".cif")[0])
        scavenger_lmajor_ids.add(file.split("/")[-1].split(".cif")[0].split("_")[1])
    con = sqlite3.connect(db_dir)
    cur = con.cursor()

    cur.execute("SELECT job_id FROM esm_summaries")
    db_job_ids = set([item[0] for item in cur.fetchall()])
    unwritten_job_ids = (beacon_job_ids | scavenger_job_ids) - db_job_ids

    print(f"Found {len(unwritten_job_ids)} job_ids from the input file that are missing in the database")

    con.close()

    # separate to find which ids belong to which partition directory
    files_to_delete = []
    for job_id in unwritten_job_ids:
        human_id, lmajor_id = job_id.split("_")

        if lmajor_id in scavenger_lmajor_ids:
            #check if json file exists
            json_file_path = Path(f"{scavenger_result_dir}/{human_id}/{lmajor_id}/{job_id}_metrics.json")

            if json_file_path.exists():
                json_to_db(json_file_path, db_dir)
            else:
                files_to_delete.append(Path(f"{scavenger_result_dir}/{human_id}/{lmajor_id}/{job_id}/{job_id}.cif"))

        elif lmajor_id in beacon_lmajor_ids:
            #check if json file exists
            json_file_path = Path(f"{beacon_result_dir}/{human_id}/{lmajor_id}/{job_id}_metrics.json")

            if json_file_path.exists():
                json_to_db(json_file_path, db_dir)
            else:
                files_to_delete.append(Path(f"{beacon_result_dir}/{human_id}/{lmajor_id}/{job_id}/{job_id}.cif"))

    print(f"{len(files_to_delete)} files were not found in the database and alternate metric files could not be found. ")

    if len(files_to_delete) == 0:
        return

    print(f"Would you like to delete these files? (y/N)")
    choice = input()

    if choice == "y":
        for file in files_to_delete:
            file.unlink(missing_ok=True)
            print(f"Deleted {file}")

        print(f"{len(files_to_delete)} files were deleted.")
    else:
        print("No files were deleted.")

def json_to_db(input_json, db_dir):
    """
    Adds input_json file into the database db_dir. If the metrics 
    were successfully added or the job_id already exists in the database, 
    the json file is moved to input_dir/../old_json_files.


    Arguments:
        input_json - file path of json file containing metrics for one ESMFold job.
        db_dir - file path of directory to add to.
    """

    Path(f"{input_json.parent}/../../../old_json_files").mkdir(exist_ok=True)
    
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
    

    with open(input_json) as f:
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

        # move json file to results/old_json_files
        basename = input_json.name
        new_path = Path(f"{input_json.parent}/../../../old_json_files/{basename}")
        input_json.replace(new_path)
        print(f"{input_json.name} moved to {new_path.parent.resolve()}")

    except Exception as error:
        # return to previous state if error occurred during insertion
        cur.execute("ROLLBACK;")
        con.close()
        warnings.simplefilter("always")

        if str(error) == "UNIQUE constraint failed: esm_summaries.job_id":
            print(f"Metrics from {data.get("job_id")} weren't uploaded because the job_id already exists in the database.")

            # move json file to results/old_json_files
            basename = input_json.name
            new_path = Path(f"{input_json.parent}/../../../old_json_files/{basename}")
            input_json.replace(new_path)
            print(f"{input_json.name} moved to {new_path.parent.resolve()}")

        else:
            warnings.warn(f"Metrics from {data.get("job_id")} could not be written to SQL database due to: {error}.", RuntimeWarning)
    else:
        cur.execute("COMMIT;")
        con.close()
        print(f"Wrote {basename} to {db_dir}")


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} <beacon_result_dir> <scavenger_result_dir> <db_dir>")
    filter_output_list(sys.argv[1], sys.argv[2], sys.argv[3])



