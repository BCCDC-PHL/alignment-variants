#!/bin/bash

set -eo pipefail

sed -i 's/cpus = 8/cpus = 4/g' nextflow.config
sed -i 's/cpus = 12/cpus = 4/g' nextflow.config
sed -i 's/cpus = 16/cpus = 4/g' nextflow.config
sed -i 's/cpus = 24/cpus = 4/g' nextflow.config
sed -i 's/cpus = 24/cpus = 4/g' nextflow.config
sed -i "s/memory = '36G'/memory = '2G'/g" nextflow.config

nextflow run main.nf \
	 -profile conda \
	 --cache ${HOME}/.conda/envs \
	 --fastq_input .github/data/fastq \
	 --outdir .github/data/test_output \
	 --min_depth 5 \
	 --ref .github/data/refs/NC_045512.2.fasta \
	 -with-report .github/data/test_output/nextflow_report.html \
 	 -with-trace .github/data/test_output/nextflow_trace.tsv
