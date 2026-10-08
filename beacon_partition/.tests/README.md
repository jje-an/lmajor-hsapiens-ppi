This contains a unit test that Snakemake generated and I edited to work with these workflows. 

Once in beacon_partition directory, run with:
`pytest .tests/unit`

## Differences from the automatically generated test:

1. Added inputs/scripts to the test's working directory
2. Changed temporary directory to be in this directory instead of in /tmp/
    1. /tmp/ is not shared across nodes, so this was needed to test with Slurm
3. Changed byte by byte output comparison to a comparison of metrics from the database instead.
    1. Test will fail if any of the metrics differ by more than 10 for plddt or.1 for the other metrics.
4. Changed from local execution to Slurm execution on each respective partition.

