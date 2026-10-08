"""
Rule test code for unit testing of rules generated with Snakemake 9.23.1.
"""

"""
I had to edit this a little to add the script files and input files correctly.
Also, I can't compare ESMFold2 outputs byte by byte because .cif files 
are always slightly different even with the same seed. I have no idea why.
I have config.yaml file set so that it will do a test run with 3 leishmania
and 3 hsapiens transcripts.
Test should pass if:
1. the snakemake command ran successfully and returned an exit code of 0.
2. the resulting database has metrics similar to the expected metrics 
    (not differing by more than 10 for plddt or .1 for the other metrics.)
-Jean
"""

import os
import sys
import shutil
import tempfile
from pathlib import Path
from subprocess import check_output
import subprocess
import glob
sys.path.insert(0, os.path.dirname(__file__))


def test_fold_pairwise(conda_prefix):

    with tempfile.TemporaryDirectory(dir=".tests/unit/tmp") as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        config_path = Path(".tests/unit/fold_pairwise/config")
        data_path = Path(".tests/unit/fold_pairwise/data")
        expected_path = Path(".tests/unit/fold_pairwise/expected")
        script_path = Path(".tests/unit/scripts")
        input_path = Path(".tests/unit/inputs")

        # Copy config to the temporary workdir.
        shutil.copytree(config_path, workdir)

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir, dirs_exist_ok=True)

        # Copy scripts to the temporary workdir.
        shutil.copytree(script_path, workdir / "scripts", dirs_exist_ok=True)

        # Copy inputs to the temporary workdir.
        shutil.copytree(input_path, workdir / "inputs", dirs_exist_ok=True)

        # Run the test job.
        check_output(
            [
                "python",
                "-m",
                "snakemake",
                "--snakefile",
                "snakefile",
                "-f",
                "--notemp",
                "--show-failed-logs",
                "--keep-going",
                "-j1",
                "--target-files-omit-workdir-adjustment",
                "--allowed-rules",
                "fold_pairwise",
                "--configfile",
                "./.tests/unit/fold_pairwise/config/config.yaml",
                "--executor",
                "slurm",
                "--default-resources",
                "slurm_partition=beacon",
                "qos=default",
                "--directory",
                workdir,
            ]
            + conda_prefix
        )

        # delete previous test's database file if exists
        old_db_path = Path('.tests/unit/fold_pairwise/data/lmajor_hsapiens_ppi.db')
        old_db_path.unlink(missing_ok=True)

        #shutil.copy(workdir / "lmajor_hsapiens_ppi.db", data_path)
        #shutil.copytree(workdir / Path(".snakemake/slurm_logs/rule_fold_pairwise"), data_path)

        files = glob.glob('./**/*', recursive=True)
        for file in files:
            print(file)


        # Check the output database file with the expected database file.
        # Checks each metric to make sure it doesn't differ from the expected
        # by more than 10 for plddt and .1 for the other metrics.
        import common
        common.OutputChecker(data_path, expected_path, workdir).check()


