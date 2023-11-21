# alignment-variants

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
  fastq_long --> nanoq_pre_filter(nanoq)
  fastq_long --> filtlong(filtlong)
  filtlong --> nanoq_post_filter(nanoq)
  ref --> bwa(bwa)
  fastp --> bwa
  filtlong --> minimap2(minimap2)
```

## Outputs

TBD.

