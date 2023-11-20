#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { fastp }                     from './modules/align_variants.nf'
include { index_ref }                 from './modules/align_variants.nf'
include { bwa_mem }                   from './modules/align_variants.nf'
include { minimap2 }                  from './modules/align_variants.nf'
include { qualimap_bamqc }            from './modules/align_variants.nf'
include { mpileup }                   from './modules/align_variants.nf'
include { generate_low_coverage_bed } from './modules/align_variants.nf'
include { calculate_gene_coverage }   from './modules/align_variants.nf'
include { pipeline_provenance }       from './modules/provenance.nf'
include { collect_provenance }        from './modules/provenance.nf'

workflow {

    ch_start_time = Channel.of(LocalDateTime.now())
    ch_pipeline_name = Channel.of(workflow.manifest.name)
    ch_pipeline_version = Channel.of(workflow.manifest.version)

    ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))

    if (params.samplesheet_input != 'NO_FILE') {
	ch_illumina_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
	ch_nanopore_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['LONG']] }
    } else {
	ch_illumina_fastq = Channel.fromFilePairs( params.fastq_illumina_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
	ch_nanopore_fastq = Channel.fromPath( params.fastq_nanopore_search_path ).map{ it -> [it.split('_')[0], it] }.unique{ it -> it[0] }
    }

    ch_ref = Channel.fromPath(params.ref)

    main:
    fastp(ch_illumina_fastq)

    ch_indexed_ref = index_ref(ch_ref)
    
    bwa_mem(fastp.out.reads.combine(ch_indexed_ref))

    ch_bwa_alignment = bwa_mem.out.alignment

    qualimap_bamqc(ch_bwa_alignment)

    ch_depths = mpileup(ch_bwa_alignment.combine(ch_ref))

    generate_low_coverage_bed(ch_depths)

    if (params.collect_outputs) {
	fastp.out.fastp_csv.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_fastp.csv", storeDir: "${params.outdir}")
	qualimap_bamqc.out.genome_results.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_qualimap_bamqc.csv", storeDir: "${params.outdir}")
    }
    
    ch_provenance = fastp.out.provenance
    ch_provenance = ch_provenance.join(bwa_mem.out.provenance).map{ it -> [it[0], [it[1], it[2]]] }

    ch_provenance = ch_provenance.join(ch_illumina_fastq.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] ] }

    collect_provenance(ch_provenance)
  
}
