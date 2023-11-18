#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { fastp }                     from './modules/align_variants.nf'
include { index_ref }                 from './modules/align_variants.nf'
include { bwa_mem }                   from './modules/align_variants.nf'
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
	ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
    } else {
	ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    }

    ch_ref = Channel.fromPath(params.ref)

    main:
    fastp(ch_fastq)

    ch_indexed_ref = index_ref(ch_ref)
    
    bwa_mem(fastp.out.reads.combine(ch_indexed_ref))

    ch_alignment = bwa_mem.out.alignment

    qualimap_bamqc(ch_alignment)

    ch_depths = mpileup(ch_alignment.combine(ch_ref))

    generate_low_coverage_bed(ch_depths)

    ch_provenance = fastp.out.provenance
    ch_provenance = ch_provenance.join(bwa_mem.out.provenance).map{ it -> [it[0], [it[1], it[2]]] }

    ch_provenance = ch_provenance.join(ch_fastq.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] ] }

    collect_provenance(ch_provenance)
  
}
