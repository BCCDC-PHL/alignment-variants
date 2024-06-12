#!/bin/bash

set -eo pipefail

# Check for a sign that we're in the GitHub Actions environment.
# Prevents these settings from being applied in other environments.
if [ -z "${GITHUB_ACTIONS}" ]; then
    echo "Not running in GitHub Actions environment."
else
    echo "Running in GitHub Actions environment."
    echo "Adjusting nextflow configuration for GitHub Actions environment."
    sed -i 's/cpus = 8/cpus = 4/g' nextflow.config
    sed -i 's/cpus = 12/cpus = 4/g' nextflow.config
    sed -i 's/cpus = 16/cpus = 4/g' nextflow.config
    sed -i 's/cpus = 24/cpus = 4/g' nextflow.config
    sed -i 's/cpus = 24/cpus = 4/g' nextflow.config
    sed -i "s/memory = '36G'/memory = '2G'/g" nextflow.config
fi


nextflow run main.nf \
	 -profile conda \
	 -resume \
	 --cache ${HOME}/.conda/envs \
	 --fastq_input .github/data/fastq \
	 --outdir .github/data/test_output \
	 --min_depth 5 \
	 --ref .github/data/refs/NC_000962.3.fa \
	 --collect_outputs \
	 -with-report .github/data/test_output/nextflow_report.html \
 	 -with-trace .github/data/test_output/nextflow_trace.tsv
