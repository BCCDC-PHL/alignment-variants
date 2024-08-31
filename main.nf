#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { hash_files as hash_ref }         from './modules/hash_files.nf'
include { hash_files as hash_fastq_short } from './modules/hash_files.nf'
include { hash_files as hash_fastq_long }  from './modules/hash_files.nf'
include { fastp }                          from './modules/short_read_qc.nf'
include { filtlong }                       from './modules/long_read_qc.nf'
include { nanoq as nanoq_pre_filter }      from './modules/long_read_qc.nf'
include { nanoq as nanoq_post_filter }     from './modules/long_read_qc.nf'
include { merge_nanoq_reports }            from './modules/long_read_qc.nf'
include { index_ref }                      from './modules/alignment_variants.nf'
include { bwa_mem }                        from './modules/alignment_variants.nf'
include { minimap2 }                       from './modules/alignment_variants.nf'
include { freebayes }                      from './modules/alignment_variants.nf'
include { qualimap_bamqc }                 from './modules/alignment_variants.nf'
include { samtools_mpileup }               from './modules/alignment_variants.nf'
include { samtools_stats }                 from './modules/alignment_variants.nf'
include { combine_alignment_qc }           from './modules/alignment_variants.nf'
include { generate_low_coverage_bed }      from './modules/alignment_variants.nf'
include { percent_coverage_by_depth }      from './modules/alignment_variants.nf'
include { plot_coverage }                  from './modules/alignment_variants.nf'
include { pipeline_provenance }            from './modules/provenance.nf'
include { collect_provenance }             from './modules/provenance.nf'

workflow {
    ch_workflow_metadata = Channel.value([
	workflow.sessionId,
	workflow.runName,
	workflow.manifest.name,
	workflow.manifest.version,
	workflow.start,
    ])

    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)

    if (params.samplesheet_input != 'NO_FILE') {
	ch_illumina_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }.filter{ it -> it[1] != null || it[2] != null }
	ch_nanopore_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['LONG']] }.filter{ it -> it[1] != null }
	ch_ref = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['REF']] }
    } else {
	ch_illumina_fastq = Channel.fromFilePairs( params.fastq_illumina_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
	ch_nanopore_fastq = Channel.fromPath( params.fastq_nanopore_search_path ).map{ it -> [it.getName().split('_')[0], it] }.unique{ it -> it[0] }
    }


    main:
    ch_illumina_sample_ids = ch_illumina_fastq.map{ it -> it[0] }
    ch_nanopore_sample_ids = ch_nanopore_fastq.map{ it -> it[0] }
    ch_sample_ids = ch_illumina_sample_ids.concat(ch_nanopore_sample_ids).unique()

    ch_provenance = ch_sample_ids

    if (params.ref != 'NO_FILE') {
	ch_ref = ch_sample_ids.combine(Channel.fromPath(params.ref))
    }

    hash_ref(ch_ref.combine(Channel.of("ref-fasta")))
    hash_fastq_short(ch_illumina_fastq.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq-input-short")))
    hash_fastq_long(ch_nanopore_fastq.combine(Channel.of("fastq-input-long")))
    
    ch_indexed_ref = index_ref(ch_ref)

    fastp(ch_illumina_fastq)

    if (! params.align_untrimmed_reads) {	
	ch_illumina_reads_to_align = fastp.out.trimmed_reads
    } else {
	ch_illumina_reads_to_align = ch_illumina_fastq
    }

    nanoq_pre_filter(ch_nanopore_fastq.combine(Channel.of("pre_filter")))

    filtlong(ch_nanopore_fastq)

    nanoq_post_filter(filtlong.out.filtered_reads.combine(Channel.of("post_filter")))

    merge_nanoq_reports(nanoq_pre_filter.out.report.join(nanoq_post_filter.out.report))

    if (! params.align_unfiltered_reads) {
	ch_nanopore_reads_to_align = filtlong.out.filtered_reads
    } else {
	ch_nanopore_reads_to_align = ch_nanopore_fastq
    }    

    bwa_mem(ch_illumina_reads_to_align.join(ch_indexed_ref))

    ch_bwa_alignment = bwa_mem.out.alignment

    minimap2(ch_nanopore_reads_to_align.join(ch_indexed_ref))

    ch_minimap2_alignment = minimap2.out.alignment

    ch_alignments = ch_bwa_alignment.concat(ch_minimap2_alignment)

    qualimap_bamqc(ch_alignments)

    freebayes(ch_alignments.join(ch_ref))

    samtools_mpileup(ch_alignments.join(ch_ref))

    samtools_stats(ch_alignments)

    combine_alignment_qc(qualimap_bamqc.out.alignment_qc.join(samtools_stats.out.stats_summary_csv, by: [0, 1]))

    ch_depths = samtools_mpileup.out.depths

    generate_low_coverage_bed(ch_depths)

    percent_coverage_by_depth(ch_depths)

    plot_coverage(ch_depths.join(ch_ref))

    // Collect multi-sample outputs
    if (params.collect_outputs) {
	fastp.out.fastp_csv.map{ it -> it[1] }.collectFile(
	    keepHeader: true,
	    sort: { it.text },
	    name: "${params.collected_outputs_prefix}_fastp.csv",
	    storeDir: "${params.outdir}"
	)
	qualimap_bamqc.out.alignment_qc.map{ it -> it[2] }.collectFile(
	    keepHeader: true,
	    sort: { it.text },
	    name: "${params.collected_outputs_prefix}_qualimap_alignment_qc.csv",
	    storeDir: "${params.outdir}"
	)
	samtools_stats.out.stats_summary_csv.map{ it -> it[2] }.collectFile(
	    keepHeader: true,
	    sort: { it.text },
	    name: "${params.collected_outputs_prefix}_samtools_stats_summary.csv",
	    storeDir: "${params.outdir}"
	)
	combine_alignment_qc.out.map{ it -> it[2] }.collectFile(
	    keepHeader: true,
	    sort: { it.text },
	    name: "${params.collected_outputs_prefix}_combined_alignment_qc.csv",
	    storeDir: "${params.outdir}"
	)
    }

    // Collect Provenance
    // The basic idea is to build up a channel with the following structure:
    // [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
    // At each step, we add another provenance file to the list using the << operator...
    // ...and then concatenate them all together in the 'collect_provenance' process.
    ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it ->            [it[0], [it[1]]] }
    ch_provenance = ch_provenance.join(hash_ref.out.provenance).map{ it ->              [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(hash_fastq_short.out.provenance).map{ it ->      [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it ->                 [it[0], it[1] << it[2]] }
    if (params.align_long_reads) {
	ch_provenance = ch_provenance.join(hash_fastq_long.out.provenance).map{ it ->   [it[0], it[1] << it[2]] }
    }
    ch_provenance = ch_provenance.join(bwa_mem.out.provenance).map{ it ->               [it[0], it[1] << it[2]] }
    if (params.align_long_reads) {
	ch_provenance = ch_provenance.join(nanoq_pre_filter.out.provenance).map{ it ->  [it[0], it[1] << it[2]] }
	ch_provenance = ch_provenance.join(filtlong.out.provenance).map{ it ->          [it[0], it[1] << it[2]] }
	ch_provenance = ch_provenance.join(nanoq_post_filter.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
	ch_provenance = ch_provenance.join(minimap2.out.provenance).map{ it ->          [it[0], it[1] << it[2]] }
    }
    ch_provenance = ch_provenance.join(qualimap_bamqc.out.provenance).map{ it ->        [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(samtools_mpileup.out.provenance).map{ it ->      [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(samtools_stats.out.provenance).map{ it ->        [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(freebayes.out.provenance).map{ it ->             [it[0], it[1] << it[2]] }

    collect_provenance(ch_provenance)
  
}
