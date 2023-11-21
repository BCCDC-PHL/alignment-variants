#!/bin/bash

set -eo pipefail

export PATH=/opt/miniconda3/bin:$PATH
export PATH=/opt/nextflow/bin:$PATH

# write test log as github Action artifact
echo Nextflow run current PR >> artifacts/test_artifact.log
NXF_VER=20.10.0 nextflow -quiet run ./main.nf \
       -profile conda \
       --cache ~/.conda/envs \
       --fastq_input $PWD/.github/data/fastqs/ \
       --ref $PWD/.github/data/refs/MN908947.3/MN908947.3.fa \
       --outdir $PWD/output \
       -with-trace $PWD/output/trace.tsv

cp .nextflow.log artifacts/pull_request.nextflow.log
cp output/trace.tsv artifacts/pull_request.trace.tsv
cp -r output artifacts/pull_request_results

# run tests against previous previous_release to compare outputs 
git clone https://github.com/dfornika/alignment-variants.git previous_release 
cd previous_release
git checkout 097e79740b58e0b171d8935a60fa13921823541a

# the github runner only has 2 cpus available, so replace for that commit required:
sed -i s'/cpus = 24/cpus = 2/'g nextflow.config
sed -i s'/cpus = 12/cpus = 2/'g nextflow.config
sed -i s'/cpus = 16/cpus = 2/'g nextflow.config
sed -i s"/memory = '36G'/memory = '256M'/"g nextflow.config

echo Nextflow run previous release.. >> ../artifacts/test_artifact.log
NXF_VER=20.10.0 nextflow -quiet run ./main.nf \
       -profile conda \
       --cache ~/.conda/envs \
       --fastq_input $PWD/../.github/data/fastqs/MN908947.3 \
       --ref $PWD/../.github/data/refs/MN908947.3/MN908947.3.fa \
       --outdir $PWD/output \
       -with-trace $PWD/output/trace.tsv

cp .nextflow.log ../artifacts/previous_release.nextflow.log
cp -r output ../artifacts/previous_release_results

cd ..

# exclude files from comparison
# and list differences
echo "Compare ouputs of current PR vs those of previous release.." >> artifacts/test_artifact.log
find results ./previous_release/results \
     -name "*.fq.gz" \
     -o -name "*.bam" \
     -o -name "*.bam.bai" \
     -o -name "*.vcf" \
    | xargs rm -rf
if ! git diff --stat --no-index results ./previous_release/results > diffs.txt ; then
  echo "test failed: differences found between PR and previous release" >> artifacts/test_artifact.log
  echo "see diffs.txt" >> artifacts/test_artifact.log 
  cp diffs.txt artifacts/  
  exit 1
else
  echo "no differences found between PR and previous release" >> artifacts/test_artifact.log
fi

# clean-up for following tests
rm -rf previous_release && rm -rf results && rm -rf work && rm -rf .nextflow*
