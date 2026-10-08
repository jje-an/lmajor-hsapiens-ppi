# Snakemake workflow for ESMFold2 host-pathogen interaction screening

A Snakemake workflow that uses ESMFold2 to predict the structures of protein pairs at scale, for use on compute clusters with the Slurm workload manager. It was built to screen *Leishmania major* proteins against the human proteome for candidate protein-protein interactions.

## Background

The 100 *L. major* genes used as input are the most upregulated genes during macrophage infection at 4 hours post-infection, taken from [Fernandes et al., mBio 2016](https://journals.asm.org/doi/10.1128/mbio.00027-16). Each one is paired against the Ensembl canonical transcript of human protein-coding genes. Every pair is folded as a complex, and the confidence metrics (ipTM in particular) are used to rank pairs as candidate interactions.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/en/stable/) (version 9.23.1)
- [Snakemake's Slurm executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) (version 2.7.1)
- [ESM](https://github.com/Biohub/esm) with ESMFold2 support, and a `transformers` build that includes `ESMFold2Model`
- The `biohub/ESMFold2-Fast` and `biohub/ESMC-6B` model weights from Hugging Face
- A GPU with at least 48 GB of memory
- [pytest](https://docs.pytest.org/) (only needed to run the unit tests)

## Repository Structure

The workflow is split into two nearly identical directories, one per Slurm partition, so that two instances of Snakemake can run at the same time. Snakemake cannot set a different QoS for each partition in a single run, and it locks its working directory while running, so each partition gets its own copy of the workflow, its own profile, and its own share of the input IDs.

The *L. major* IDs are split 1:3 between the two partitions because the beacon partition runs about 5 jobs at a time and the scavenger partition usually runs about 15.

```
.
├── beacon_partition        # run on the beacon partition (last 25 L. major IDs)
│   ├── config.yaml
│   ├── profile.yaml
│   ├── snakefile
│   ├── inputs
│   ├── scripts
│   └── .tests
├── scavenger_partition     # run on the scavenger partition (top 75 L. major IDs)
│   └── ...
└── scripts                 # helper scripts for metrics and cleanup
```

## Configuration

Configuration options are specified in `config.yaml`, where you can specify the mode, input filepaths, library filepaths, and the path to the metrics database. Snakemake's command line options (executor, partition, QoS, per-job resources, and batching) are specified in `profile.yaml`, so they do not need to be typed out for each run.

In pairwise mode, `idfile1` is the list of human transcript IDs and `idfile2` is the list of *L. major* IDs. Both partitions use `filtered_transcript_ids.txt` for `idfile1`. For `idfile2`, the beacon partition uses `lmajor_last25.txt` and the scavenger partition uses `lmajor_top75.txt`.

### Modes

- Separate - Takes an input of a text file containing Ensembl IDs and submits individual ESMFold2 jobs to fold each individual protein sequence.
- Pairwise - Takes two input text files containing IDs and submits jobs of each pair of proteins from the first and second input file. For example, if file 1 had proteins A, B, and C, and file 2 had proteins D, E, and F, there would be 9 results: pairs AD, AE, AF, BD, BE, BF, and CD, CE, CF.

### Inputs

Input ID lists are located in `inputs/`.

- `lmajor_top100.txt` - The 100 *L. major* gene IDs. `lmajor_top75.txt` and `lmajor_last25.txt` are the same list split between the two partitions.
- `filtered_transcript_ids.txt` - Human canonical transcript IDs (18,157 transcripts). Sequences longer than 750 residues are removed, since they run out of memory on a 48 GB GPU.
- `hsapiens_known_interactors1.txt`, `hsapiens_known_interactors2.txt` - Pairs of human proteins with known interactions, used for testing.
- The smaller files (`lmajor_1.txt`, `lmajor_3.txt`, `lmajor_top5.txt`, `hsapiens_transcriptID_3.txt`, etc.) are subsets for test runs and runtime estimates.

## Usage

It is recommended to run the snakemake command in a terminal multiplexer such as Tmux or screen so that the execution will not be canceled if SSH connection drops.
Run the commands from inside the partition directory you are using. Do a dry run first to check for errors and ensure correct target files.

```bash
cd scavenger_partition
snakemake -n
```

Both workflows are run the same way. `profile.yaml` contains all of the command line options:

```bash
snakemake --profile profile.yaml
```

On the scavenger partition, each job needs about 10 minutes of walltime, 1 CPU, 1 GPU, and 20 GB of RAM.

### Batching

The scavenger partition's run is split into 3 batches with Snakemake's `--batch` option so that the DAG loads faster. After the first batch is done, edit the batch setting in `profile.yaml` (see the comments in that file) and run the second batch, then the third. If the DAG still takes too long to load, try increasing the number of batches.

### Before running both workflows

Start a job of your own on the beacon partition before launching the workflows, so that you have somewhere to test and debug. Otherwise Snakemake will take up all of the available jobs on that partition.

## Logs

Logs will be saved in `.snakemake/log/`. When the Slurm executor plugin is used, full error messages are located in `.snakemake/slurm_logs/`.

## Metrics

Each job writes its confidence metrics to a SQLite database (`lmajor_hsapiens_ppi.db`, table `esm_summaries`) as soon as the fold finishes, so there is no separate step to collect results. If the database write fails, the job raises a warning and saves its metrics to a `.json` file instead.

Do not read from the database while the workflow is running. WAL mode cannot be used on a network file system, so reading locks the database and causes running jobs to fail.

| Column | Description |
| --- | --- |
| `job_id` | ID of the protein pair |
| `ptm` | Predicted TM-score for the full structure, from 0 to 1 |
| `iptm` | Predicted interface TM-score, from 0 to 1. Predicted accuracy of the relative positions of the two chains |
| `plddt_mean` | Mean pLDDT across all residues, from 0 to 100 |
| `plddt` | Per-residue pLDDT |
| `pae` | Predicted aligned error for each pair of residues |
| `pair_chain_iptm` | ipTM restricted to each pair of chains |
| `residue_index` | Index of each residue within its chain |
| `entity_id` | Which chain each residue belongs to |
| `distogram` | Predicted distance distribution for each pair of residues. Off by default, since it makes the database extremely large (about 3 TB for the full run). Set `keep_distogram: true` in `config.yaml` to include it |

Various helper scripts are located in `scripts` for analyzing metrics.
- `filter_output_list.py` - Finds completed runs that are missing from the database, then uploads their `.json` metrics or deletes runs that have no metrics. Run this if the number of `.cif` files does not match the number of entries in the database.
- `esm_summaries_to_db.py` - Adds metrics from `.json` files to the database. It does not check for `.cif` files that are missing database entries.
- `iptm_histogram.py` - Creates a log scale histogram of ipTM values.
- `protein_fold_lib.py` - Helper functions for the Snakefile (located in each partition's `scripts` directory).

## Common Errors

- `srun: error: Unable to confirm allocation for job ...: Socket timed out on send/recv operation` - This is typically an error on the cluster's end and can be ignored.
- `OutOfMemoryError ... CUDA out of memory.` - Long sequences are filtered out of the input, but a few pairs may still run out of GPU memory. These can be ignored, or run separately on a node with more GPU memory (80 GB or 141 GB).

## Tests

Unit tests for the `fold_pairwise` rule are located in `.tests/unit/` in each partition directory. They were generated with Snakemake and edited to work with ESMFold2.

```bash
pytest .tests/unit/
```

ESMFold2's `.cif` files differ slightly between runs, even with the same seed and GPU, so the test does not compare outputs byte by byte. Instead, it compares the metrics written to the database against a reference run. The test fails if `plddt_mean` differs by more than 10, or if any of the 0 to 1 metrics differ by more than 0.1.

## Output Structure

All outputs will be located in `results/`. The output structure is slightly different depending on the mode used.

### Separate

Each sample's results are in their own directory, using the same names as each line in the input file. The directory contains the predicted structure as a `.cif` file.

### Pairwise

Each sample pair's directory is formatted as: `results/pairs/{idfile1_id}/{idfile2_id}/`, with the human transcript ID first and the *L. major* ID second. This directory contains ESMFold2's result directory, `{idfile1_id}_{idfile2_id}/`, which holds the predicted structure as a `.cif` file. If the database write failed, a `{idfile1_id}_{idfile2_id}_metrics.json` file is saved in the pair's directory, one level above the `.cif` file. After `filter_output_list.py` uploads it to the database, the file is moved to `results/old_json_files/`. For example, for an input of two files with proteins A, B, and C, D respectively, the directory would look like:

```
results
├── pairs
│   ├── A
│   │   ├── C
│   │   │   ├── A_C                # ESMFold2's result directory
│   │   │   │   └── A_C.cif
│   │   │   └── A_C_metrics.json   # only if the database write failed
│   │   └── D
│   │       └── ...
│   └── B
│       ├── C
│       │   └── ...
│       └── D
│           └── ...
├── old_json_files                # fallback .json files already uploaded to the database
└── separate
```

For example, a real output path looks like `results/pairs/ENST00000000233/LmjF.08.0670/ENST00000000233_LmjF.08.0670/ENST00000000233_LmjF.08.0670.cif`.