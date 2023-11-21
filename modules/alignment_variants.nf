process index_ref {

    tag { ref_filename }

    input:
    path(ref)

    output:
    tuple path('ref.fa'), path('ref.fa.*')

    script:
    ref_filename = ref.getName()
    """
    ln -s ${ref} ref.fa
    bwa index ref.fa
    """
}


process bwa_mem {

    tag { sample_id }
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}.{bam,bam.bai}"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(ref), path(ref_index)

    output:
    tuple val(sample_id), path("${sample_id}.{bam,bam.bai}"), emit: alignment
    tuple val(sample_id), path("${sample_id}_bwa_mem_provenance.yml"), emit: provenance
    
    script:
    bwa_threads = task.cpus - 8
    """
    printf -- "- process_name: bwa_mem\\n" >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "  tools:\\n"                    >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "    - tool_name: bwa\\n"        >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      tool_version: \$(bwa 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n"      >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "    - tool_name: samtools\\n"   >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      tool_version: \$(samtools 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n" >> ${sample_id}_bwa_mem_provenance.yml

    bwa mem \
	-t ${bwa_threads} \
	-R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
	${ref} \
	${reads_1} \
	${reads_2} \
	| samtools view -@ 2 -h -F 1540 \
	| samtools sort -@ 2 -l 0 -m 1000M -n \
        | samtools fixmate -mr - - \
	| samtools sort -@ 2 -l 0 -m 1000M \
	| samtools markdup -r - - \
	> ${sample_id}.bam

    samtools index ${sample_id}.bam    
    """
}


process minimap2 {

    tag { sample_id }
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}.{bam,bam.bai}"

    input:
    tuple val(sample_id), path(reads), path(ref), path(ref_index)

    output:
    tuple val(sample_id), path("${sample_id}.{bam,bam.bai}"), emit: alignment
    tuple val(sample_id), path("${sample_id}_minimap2_provenance.yml"), emit: provenance
    
    script:
    """
    printf -- "- process_name: \"minimap2\"\\n" >> ${sample_id}_minimap2_provenance.yml
    printf -- "  tools:\\n"                     >> ${sample_id}_minimap2_provenance.yml
    printf -- "    - tool_name: minimap2\\n"    >> ${sample_id}_minimap2_provenance.yml
    printf -- "      tool_version: \$(minimap2 --version)\\n"  >> ${sample_id}_minimap2_provenance.yml
    
    minimap2 \
	-t ${task.cpus} \
	-ax map-ont \
	${ref} \
	${reads} \
	> ${sample_id}.bam
    """
}


process qualimap_bamqc {

    tag { sample_id }

    publishDir  "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_qualimap_alignment_qc.csv"
    publishDir  "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_qualimap_report.pdf"

    input:
    tuple val(sample_id), file(alignment)

    output:
    tuple val(sample_id), path("${sample_id}_qualimap_alignment_qc.csv"), emit: genome_results
    tuple val(sample_id), path("${sample_id}_qualimap_report.pdf"), emit: report
    
    script:
    """
    qualimap bamqc \
	--paint-chromosome-limits \
	--collect-overlap-pairs \
	--cov-hist-lim ${params.qualimap_coverage_histogram_limit} \
	--output-genome-coverage ${sample_id}_genome_coverage.txt \
	-nt ${task.cpus} \
	-bam ${alignment[0]} \
	-outformat PDF \
	--outdir ${sample_id}_bamqc

    qualimap_bamqc_genome_results_to_csv.py \
	-s ${sample_id} \
	${sample_id}_bamqc/genome_results.txt \
	> ${sample_id}_qualimap_alignment_qc.csv

    cp ${sample_id}_bamqc/report.pdf ${sample_id}_qualimap_report.pdf
    """
}


process mpileup {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_depths.tsv"

    input:
    tuple val(sample_id), path(alignment), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}_depths.tsv")

    script:
    """
    samtools faidx ${ref}

    printf "chrom\tpos\tref\tdepth\n" > ${sample_id}_depths.tsv

    samtools mpileup -a \
      --fasta-ref ${ref} \
      --min-BQ 0 \
      --count-orphans \
      ${alignment[0]} | cut -f 1-4 >> ${sample_id}_depths.tsv
    """
}


process generate_low_coverage_bed {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_low_coverage_regions.bed"

    input:
    tuple val(sample_id), path(depths)

    output:
    tuple val(sample_id), path("${sample_id}_low_coverage_regions.bed")

    script:
    """
    generate_low_coverage_bed.py \
      --input ${depths} \
      --threshold ${params.min_depth} \
      > ${sample_id}_low_coverage_regions.bed
    """
}


process freebayes {
    
    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_freebayes.vcf"

    input:
    tuple val(sample_id), path(alignment), path(ref), path(ref_index)

    output:
    tuple val(sample_id), path("${sample_id}_freebayes.vcf")

    script:
    """
    freebayes \
	--fasta-reference ${ref} \
	--bam ${alignment[0]} \
	--ploidy 1 \
	--min-base-quality ${params.min_base_qual_for_variant_calling} \
	--min-mapping-quality ${params.min_mapping_qual_for_variant_calling} \
	--min-alternate-count 2 \
	--min-alternate-fraction ${params.min_af_for_variant_calling} \
	--report-genotype-likelihood-max \
	> ${sample_id}_freebayes.vcf
    """
}
