process fastp {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp.{json,csv}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: fastp_json
    tuple val(sample_id), path("${sample_id}_fastp.csv"), emit: fastp_csv
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}_fastp_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: fastp\\n" > ${sample_id}_fastp_provenance.yml
    printf -- "  tool_name: fastp\\n  tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastp_provenance.yml

    fastp \
	--cut_tail \
	-i ${reads_1} \
	-I ${reads_2} \
	-o ${sample_id}_trimmed_R1.fastq.gz \
	-O ${sample_id}_trimmed_R2.fastq.gz

    mv fastp.json ${sample_id}_fastp.json
    fastp_json_to_csv.py -s ${sample_id} ${sample_id}_fastp.json > ${sample_id}_fastp.csv
    """
}

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
    """

    bwa mem \
	-t ${task.cpus} \
	-R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
	${ref} \
	${reads_1} \
	${reads_2} \
	| samtools sort -n \
        | samtools fixmate -mr - - \
	| samtools sort - \
	| samtools markdup -r - - \
	> ${sample_id}.bam

    samtools index ${sample_id}.bam

    printf -- "- process_name: \"bwa mem\"\\n" >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "  tools:\\n"                    >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "    - tool_name: bwa\\n"        >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      tool_version: \$(bwa 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n"      >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "    - tool_name: samtools\\n"   >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      tool_version: \$(samtools 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n" >> ${sample_id}_bwa_mem_provenance.yml
    """
}


process qualimap_bamqc {

    tag { sample_id }

    publishDir  "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_qualimap_alignment_qc.csv"

    input:
    tuple val(sample_id), file(alignment)

    output:
    tuple val(sample_id), path("${sample_id}_qualimap_alignment_qc.csv"), emit: genome_results
    
    script:
    """
    qualimap bamqc -bam ${alignment[0]} --outdir ${sample_id}_bamqc
    qualimap_bamqc_genome_results_to_csv.py -s ${sample_id} ${sample_id}_bamqc/genome_results.txt > ${sample_id}_qualimap_alignment_qc.csv
    """
}


process mpileup {

    tag { sample_id }

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


process calculate_gene_coverage {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_resistance_gene_coverage.csv"

    input:
    tuple val(sample_id), path(depths), path(resistance_genes_bed)

    output:
    tuple val(sample_id), path("${sample_id}_resistance_gene_coverage.csv")

    script:
    """
    calculate_res_gene_depth_qc.py \
      --bed ${resistance_genes_bed} \
      --depth ${depths} \
      --threshold ${params.min_depth} \
      --output ${sample_id}_resistance_gene_coverage.csv
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
	--min-mapping-quality 30 \
	--min-base-quality 20 \
	--min-alternate-count 2 \
	--min-alternate-fraction 0.1 \
	--report-genotype-likelihood-max \
	> ${sample_id}_freebayes.vcf
    """
}
