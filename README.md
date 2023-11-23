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

If long reads are available, they can be included with the `--fastq_input_long` flag:

```
nextflow run dfornika/alignment-variants \
  --ref /path/to/ref.fa \
  --fastq_input /path/to/fastqs \
  --fastq_input_long /path/to/long_fastqs \
  --outdir /path/to/outputs
```

Alternatively, a `samplesheet.csv` file can be provided:

```
nextflow run dfornika/alignment-variants \
  --ref /path/to/ref.fa \
  --samplesheet_input /path/to/samplesheet.csv \
  --outdir /path/to/outputs
```

The fields should include:

```
ID
R1
R2
```

...and if long reads are included, use the field:

```
LONG
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

```
output
├── SAMPLE-1
│   ├── SAMPLE-1_20231123140142_provenance.yml
│   ├── SAMPLE-1_fastp.csv
│   ├── SAMPLE-1_fastp.json
│   ├── SAMPLE-1_long.bam
│   ├── SAMPLE-1_long.bam.bai
│   ├── SAMPLE-1_short.bam
│   ├── SAMPLE-1_short.bam.bai
│   ├── SAMPLE-1_short_depths.tsv
│   ├── SAMPLE-1_short_freebayes.vcf
│   ├── SAMPLE-1_short_low_coverage_regions.bed
│   └── SAMPLE-1_short_percent_coverage_by_depth.csv
├── SAMPLE-2
│   ├── SAMPLE-2_20231123140112_provenance.yml
│   ├── SAMPLE-2_fastp.csv
│   ├── SAMPLE-2_fastp.json
│   ├── SAMPLE-2_long.bam
│   ├── SAMPLE-2_long.bam.bai
│   ├── SAMPLE-2_short.bam
│   ├── SAMPLE-2_short.bam.bai
│   ├── SAMPLE-2_short_depths.tsv
│   ├── SAMPLE-2_short_freebayes.vcf
│   ├── SAMPLE-2_short_low_coverage_regions.bed
│   └── SAMPLE-2_short_percent_coverage_by_depth.csv
├── collected_fastp.csv
└── collected_qualimap_bamqc.csv
```

