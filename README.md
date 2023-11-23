# alignment-variants

[![Tests](https://github.com/dfornika/alignment-variants/actions/workflows/pull_request.yml/badge.svg)](https://github.com/dfornika/alignment-variants/actions/workflows/pull_request.yml)

A pipeline to align (map) reads against a reference genome, and call variants based on the alignment.

## Usage

```
nextflow run dfornika/alignment-variants \
  --ref /path/to/ref.fa \
  --fastq_input /path/to/fastqs \
  --outdir /path/to/outputs
```

## Pipeline

```mermaid
flowchart TD
  ref[ref.fa]
  fastq_short[fastq_short]
  fastq_long[fastq_long]
  fastq_short --> fastp(fastp)
  fastq_long --> nanoq_pre_filter("nanoq (pre-filter)")
  fastq_long --> filtlong(filtlong)
  filtlong --> nanoq_post_filter("nanoq (post-filter)")
  ref --> bwa(bwa)
  fastp --> bwa
  ref --> minimap2(minimap2)
  filtlong --> minimap2
  bwa --> qualimap(qualimap_bamqc)
  minimap2 --> qualimap
```

## Outputs

TBD.

