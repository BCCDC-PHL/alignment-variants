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
  ref --> bwa(bwa_mem)
  fastp --> bwa
  bwa --> alignments[alignments]
  ref --> minimap2(minimap2)
  filtlong --> minimap2
  minimap2 --> alignments
  alignments --> qualimap(qualimap_bamqc)
  alignments --> freebayes(freebayes)
  alignments --> samtools_mpileup(samtools_mpileup)
  samtools_mpileup --> generate_low_coverage_bed(generate_low_coverage_bed)
  samtools_mpileup --> percent_coverage_by_depth(percent_coverage_by_depth)
```

## Outputs

TBD.

