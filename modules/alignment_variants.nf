process index_ref {

    tag { sample_id + ' / ' + ref_filename }

    input:
    tuple val(sample_id), path(ref)

    output:
    tuple val(sample_id), path('ref.fa'), path('ref.fa.*')

    script:
    ref_filename = ref.getName()
    """
    ln -s ${ref} ref.fa
    bwa index ref.fa
    """
}


process bwa_mem {

    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${short_long}.{bam,bam.bai}"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "ref.fa"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(ref), path(ref_index)

    output:
    tuple val(sample_id), val(short_long), path("${sample_id}_${short_long}.{bam,bam.bai}"), emit: alignment
    tuple val(sample_id), path(ref), emit: ref
    tuple val(sample_id), path("${sample_id}_bwa_mem_provenance.yml"), emit: provenance
    
    script:
    bwa_threads = task.cpus - 8
    short_long = "short"
    samtools_view_filter_flags = params.skip_alignment_cleaning ? "0" : "1548"
    samtools_fixmate_remove_secondary_and_unmapped = params.skip_alignment_cleaning ? "" : "-r"
    samtools_markdup_remove_duplicates = params.skip_alignment_cleaning ? "" : "-r"
    """
    printf -- "- process_name: bwa_mem\\n"     >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "  tools:\\n"                    >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "    - tool_name: bwa\\n"        >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      tool_version: \$(bwa 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n"      >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "    - tool_name: samtools\\n"   >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      tool_version: \$(samtools 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n" >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      subcommand: view\\n"      >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      parameters:\\n"           >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "        - parameter: -F\\n"     >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "          value: ${samtools_view_filter_flags}\\n" >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "    - tool_name: samtools\\n"   >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      tool_version: \$(samtools 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n" >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      subcommand: fixmate\\n"   >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      parameters:\\n"           >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "        - parameter: -m\\n"     >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "          value: null\\n"       >> ${sample_id}_bwa_mem_provenance.yml
    if [ ! "${params.skip_alignment_cleaning}" == "true" ]; then
    printf -- "        - parameter: -r\\n"     >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "          value: null\\n"       >> ${sample_id}_bwa_mem_provenance.yml
    fi
    printf -- "    - tool_name: samtools\\n"   >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      tool_version: \$(samtools 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n" >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      subcommand: markdup\\n"      >> ${sample_id}_bwa_mem_provenance.yml
    if [ ! "${params.skip_alignment_cleaning}" == "true" ]; then
    printf -- "      parameters:\\n"           >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "        - parameter: -r\\n"     >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "          value: null\\n"       >> ${sample_id}_bwa_mem_provenance.yml
    fi

    bwa mem \
	-t ${bwa_threads} \
	-R "@RG\\tID:${sample_id}-ILLUMINA\\tSM:${sample_id}\\tPL:ILLUMINA" \
	${ref} \
	${reads_1} \
	${reads_2} \
	| samtools view -@ 2 -h -F ${samtools_view_filter_flags} \
        | samtools sort -@ 2 -l 0 -m 1000M -n \
        | samtools fixmate -m ${samtools_fixmate_remove_secondary_and_unmapped} - - \
	| samtools sort -@ 2 -l 0 -m 1000M \
	| samtools markdup ${samtools_markdup_remove_duplicates} - - \
	> ${sample_id}_${short_long}.bam

    samtools index ${sample_id}_${short_long}.bam
    """
}


process minimap2 {

    tag { sample_id }
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${short_long}.{bam,bam.bai}"

    input:
    tuple val(sample_id), path(reads), path(ref), path(ref_index)

    output:
    tuple val(sample_id), val(short_long), path("${sample_id}_${short_long}.{bam,bam.bai}"), emit: alignment
    tuple val(sample_id), path("${sample_id}_minimap2_provenance.yml"), emit: provenance
    
    script:
    short_long = "long"
    minimap2_threads = task.cpus - 4
    samtools_view_filter_flags = params.skip_alignment_cleaning ? "0" : "1540"
    """
    printf -- "- process_name: \"minimap2\"\\n" >> ${sample_id}_minimap2_provenance.yml
    printf -- "  tools:\\n"                     >> ${sample_id}_minimap2_provenance.yml
    printf -- "    - tool_name: minimap2\\n"    >> ${sample_id}_minimap2_provenance.yml
    printf -- "      tool_version: \$(minimap2 --version)\\n"  >> ${sample_id}_minimap2_provenance.yml
    printf -- "      parameters:\\n"            >> ${sample_id}_minimap2_provenance.yml
    printf -- "        - parameter: -a\\n"      >> ${sample_id}_minimap2_provenance.yml
    printf -- "          value: null\\n"        >> ${sample_id}_minimap2_provenance.yml
    printf -- "        - parameter: -x\\n"      >> ${sample_id}_minimap2_provenance.yml
    printf -- "          value: map-ont\\n"     >> ${sample_id}_minimap2_provenance.yml
    printf -- "        - parameter: -MD\\n"     >> ${sample_id}_minimap2_provenance.yml
    printf -- "          value: null\\n"        >> ${sample_id}_minimap2_provenance.yml
    printf -- "    - tool_name: samtools\\n"    >> ${sample_id}_minimap2_provenance.yml
    printf -- "      tool_version: \$(samtools 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n"  >> ${sample_id}_minimap2_provenance.yml
    printf -- "      subcommand: view\\n"       >> ${sample_id}_minimap2_provenance.yml
    printf -- "      parameters:\\n"            >> ${sample_id}_minimap2_provenance.yml
    printf -- "        - parameter: -F\\n"      >> ${sample_id}_minimap2_provenance.yml
    printf -- "          value: ${samtools_view_filter_flags}\\n"        >> ${sample_id}_minimap2_provenance.yml

    minimap2 \
	-t ${minimap2_threads} \
	-ax map-ont \
	-R "@RG\\tID:${sample_id}-ONT\\tSM:${sample_id}\\tPL:ONT" \
	-MD \
	${ref} \
	${reads} \
	| samtools view -@ 2 -h -F ${samtools_view_filter_flags} \
	| samtools sort -@ 2 -l 0 -m 1000M -O bam \
	> ${sample_id}_${short_long}.bam

    samtools index ${sample_id}_${short_long}.bam
    """
}


process qualimap_bamqc {

    tag { sample_id + ' / ' + short_long }

    publishDir  "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${short_long}_qualimap_alignment_qc.csv"
    publishDir  "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${short_long}_qualimap_genome_results.txt"
    publishDir  "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${short_long}_qualimap_report.pdf"

    input:
    tuple val(sample_id), val(short_long), file(alignment)

    output:
    tuple val(sample_id), val(short_long), path("${sample_id}_${short_long}_qualimap_alignment_qc.csv"), emit: alignment_qc
    tuple val(sample_id), val(short_long), path("${sample_id}_${short_long}_qualimap_report.pdf"), emit: report, optional: true
    tuple val(sample_id), val(short_long), path("${sample_id}_${short_long}_qualimap_genome_results.txt"), emit: genome_results, optional: true
    tuple val(sample_id), path("${sample_id}_${short_long}_qualimap_bamqc_provenance.yml"), emit: provenance
    
    script:
    """
    printf -- "- process_name: \"qualimap_bamqc\"\\n"  >> ${sample_id}_${short_long}_qualimap_bamqc_provenance.yml
    printf -- "  tools:\\n"                            >> ${sample_id}_${short_long}_qualimap_bamqc_provenance.yml
    printf -- "    - tool_name: qualimap\\n"           >> ${sample_id}_${short_long}_qualimap_bamqc_provenance.yml
    printf -- "      tool_version: \$(qualimap bamqc | head | grep QualiMap | cut -d ' ' -f 2)\\n" >> ${sample_id}_${short_long}_qualimap_bamqc_provenance.yml
    printf -- "      parameters:\\n"                   >> ${sample_id}_${short_long}_qualimap_bamqc_provenance.yml
    printf -- "        - parameter: --collect-overlap-pairs\\n" >> ${sample_id}_${short_long}_qualimap_bamqc_provenance.yml
    printf -- "          value: null\\n"               >> ${sample_id}_${short_long}_qualimap_bamqc_provenance.yml
    printf -- "        - parameter: --cov-hist-lim\\n" >> ${sample_id}_${short_long}_qualimap_bamqc_provenance.yml
    printf -- "          value: ${params.qualimap_coverage_histogram_limit}\\n" >> ${sample_id}_${short_long}_qualimap_bamqc_provenance.yml

    # Assume qualimap exits successfully
    # If it fails we will re-assign the exit code
    # and generate an empty qualimap alignment qc
    qualimap_exit_code=0

    qualimap \
	--java-mem-size=${params.qualimap_memory} \
	bamqc \
	--paint-chromosome-limits \
	--collect-overlap-pairs \
	--cov-hist-lim ${params.qualimap_coverage_histogram_limit} \
	--output-genome-coverage ${sample_id}_${short_long}_genome_coverage.txt \
	-nt ${task.cpus} \
	-bam ${alignment[0]} \
	-outformat PDF \
	--outdir ${sample_id}_${short_long}_bamqc \
	|| qualimap_exit_code=\$?

    if [ \${qualimap_exit_code} -ne 0 ]; then
    echo "Qualimap failed with exit code \${qualimap_exit_code}"
        echo "Generating empty qualimap alignment qc"
        qualimap_bamqc_genome_results_to_csv.py \
	-s ${sample_id} \
	--read-type ${short_long} \
	--failed \
	> ${sample_id}_${short_long}_qualimap_alignment_qc.csv
    else
	qualimap_bamqc_genome_results_to_csv.py \
	-s ${sample_id} \
	--read-type ${short_long} \
	--qualimap-bamqc-genome-results ${sample_id}_${short_long}_bamqc/genome_results.txt \
	> ${sample_id}_${short_long}_qualimap_alignment_qc.csv
        cp ${sample_id}_${short_long}_bamqc/report.pdf ${sample_id}_${short_long}_qualimap_report.pdf
        cp ${sample_id}_${short_long}_bamqc/genome_results.txt ${sample_id}_${short_long}_qualimap_genome_results.txt
    fi
    """
}


process samtools_stats {

    tag { sample_id + ' / ' + short_long }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${short_long}_samtools_stats*.{txt,tsv,csv}"

    input:
    tuple val(sample_id), val(short_long), path(alignment)

    output:
    tuple val(sample_id), val(short_long), path("${sample_id}_${short_long}_samtools_stats.txt"), emit: stats
    tuple val(sample_id), val(short_long), path("${sample_id}_${short_long}_samtools_stats_summary.txt"), emit: stats_summary
    tuple val(sample_id), val(short_long), path("${sample_id}_${short_long}_samtools_stats_summary.csv"), emit: stats_summary_csv
    tuple val(sample_id), val(short_long), path("${sample_id}_${short_long}_samtools_stats_insert_sizes.tsv"), emit: insert_sizes
    tuple val(sample_id), val(short_long), path("${sample_id}_${short_long}_samtools_stats_coverage_distribution.tsv"), emit: coverage_distribution
    tuple val(sample_id),  path("${sample_id}_${short_long}_samtools_stats_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: samtools_stats\\n" >> ${sample_id}_${short_long}_samtools_stats_provenance.yml
    printf -- "  tools:\\n"                       >> ${sample_id}_${short_long}_samtools_stats_provenance.yml
    printf -- "    - tool_name: samtools\\n"      >> ${sample_id}_${short_long}_samtools_stats_provenance.yml
    printf -- "      tool_version: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_${short_long}_samtools_stats_provenance.yml
    printf -- "      subcommand: stats\\n"        >> ${sample_id}_${short_long}_samtools_stats_provenance.yml

    samtools stats \
	--threads ${task.cpus} \
	${alignment[0]} > ${sample_id}_${short_long}_samtools_stats.txt

    grep '^SN' ${sample_id}_${short_long}_samtools_stats.txt | cut -f 2-  > ${sample_id}_${short_long}_samtools_stats_summary.txt

    parse_samtools_stats_summary.py -i ${sample_id}_${short_long}_samtools_stats_summary.txt -s ${sample_id} > ${sample_id}_${short_long}_samtools_stats_summary.csv

    echo "insert_size,pairs_total,inward_oriented_pairs,outward_oriented_pairs,other_pairs" | tr ',' '\t' > ${sample_id}_${short_long}_samtools_stats_insert_sizes.tsv
    grep '^IS' ${sample_id}_${short_long}_samtools_stats.txt | cut -f 2-  >> ${sample_id}_${short_long}_samtools_stats_insert_sizes.tsv

    echo "coverage,depth" | tr ',' '\t' > ${sample_id}_${short_long}_samtools_stats_coverage_distribution.tsv
    grep '^COV' ${sample_id}_${short_long}_samtools_stats.txt | cut -f 2- >> ${sample_id}_${short_long}_samtools_stats_coverage_distribution.tsv	
    """
}


process combine_alignment_qc {

    tag { sample_id + ' / ' + short_long }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${short_long}_combined_alignment_qc.csv"

    input:
    tuple val(sample_id), val(short_long), path(qualimap_genome_results_csv), path(samtools_stats_summary_csv)

    output:
    tuple val(sample_id), val(short_long), path("${sample_id}_${short_long}_combined_alignment_qc.csv")

    script:
    """
    combine_alignment_qc.py \
	--sample-id ${sample_id} \
	--read-type ${short_long} \
	--qualimap-bamqc-genome-results ${qualimap_genome_results_csv} \
	--samtools-stats-summary ${samtools_stats_summary_csv} \
	> ${sample_id}_${short_long}_combined_alignment_qc.csv
    """
}


process samtools_mpileup {

    tag { sample_id + ' / ' + short_long }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${short_long}_depths.tsv"

    input:
    tuple val(sample_id), val(short_long), path(alignment), path(ref)

    output:
    tuple val(sample_id), val(short_long), path("${sample_id}_${short_long}_depths.tsv"), emit: depths
    tuple val(sample_id), path("${sample_id}_${short_long}_samtools_mpileup_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: samtools_mpileup\\n" >> ${sample_id}_${short_long}_samtools_mpileup_provenance.yml
    printf -- "  tools:\\n"                         >> ${sample_id}_${short_long}_samtools_mpileup_provenance.yml
    printf -- "    - tool_name: samtools\\n"        >> ${sample_id}_${short_long}_samtools_mpileup_provenance.yml
    printf -- "      tool_version: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_${short_long}_samtools_mpileup_provenance.yml
    printf -- "      subcommand: mpileup\\n"        >> ${sample_id}_${short_long}_samtools_mpileup_provenance.yml
    printf -- "      parameters:\\n"                >> ${sample_id}_${short_long}_samtools_mpileup_provenance.yml
    printf -- "        - parameter: -a\\n"          >> ${sample_id}_${short_long}_samtools_mpileup_provenance.yml
    printf -- "          value: null\\n"            >> ${sample_id}_${short_long}_samtools_mpileup_provenance.yml
    printf -- "        - parameter: --min-BQ\\n"    >> ${sample_id}_${short_long}_samtools_mpileup_provenance.yml
    printf -- "          value: 0\\n"               >> ${sample_id}_${short_long}_samtools_mpileup_provenance.yml
    printf -- "        - parameter: --count-orphans\\n" >> ${sample_id}_${short_long}_samtools_mpileup_provenance.yml
    printf -- "          value: null\\n"                >> ${sample_id}_${short_long}_samtools_mpileup_provenance.yml

    samtools faidx ${ref}

    printf "chrom\tpos\tref\tdepth\n" > ${sample_id}_${short_long}_depths.tsv

    samtools mpileup -a \
	--fasta-ref ${ref} \
	--min-BQ 0 \
	--count-orphans \
	${alignment[0]} | cut -f 1-4 >> ${sample_id}_${short_long}_depths.tsv
    """
}


process generate_low_coverage_bed {

    tag { sample_id + ' / ' + short_long }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${short_long}_low_coverage_regions.bed"

    input:
    tuple val(sample_id), val(short_long), path(depths)

    output:
    tuple val(sample_id), path("${sample_id}_${short_long}_low_coverage_regions.bed")

    script:
    """
    generate_low_coverage_bed.py \
	--input ${depths} \
	--threshold ${params.min_depth} \
	> ${sample_id}_${short_long}_low_coverage_regions.bed
    """
}

process percent_coverage_by_depth {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${short_long}_percent_coverage_by_depth.csv"

    input:
    tuple val(sample_id), val(short_long), path(depths)

    output:
    tuple val(sample_id), path("${sample_id}_${short_long}_percent_coverage_by_depth.csv")

    script:
    """
    percent_coverage_by_depth.py \
	--sample-id ${sample_id} \
	--input ${depths} \
	--max-depth ${params.coverage_by_depth_limit} \
	> ${sample_id}_${short_long}_percent_coverage_by_depth.csv
    """
}


process freebayes {
    
    tag { sample_id + ' / ' + short_long }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${short_long}_freebayes.vcf"

    input:
    tuple val(sample_id), val(short_long), path(alignment), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}_${short_long}_freebayes.vcf"), emit: variants
    tuple val(sample_id), path("${sample_id}_${short_long}_freebayes_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: freebayes\\n" >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "  tools:\\n"                         >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "    - tool_name: freebayes\\n"       >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "      tool_version: \$(freebayes --version | head -n 1 | cut -d ' ' -f 3)\\n" >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "      parameters:\\n"                >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "        - parameter: --ploidy\\n"    >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "          value: 1\\n"               >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "        - parameter: --min-base-quality\\n" >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "          value: ${params.min_base_qual_for_variant_calling}\\n" >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "        - parameter: --min-mapping-quality\\n" >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "          value: ${params.min_mapping_qual_for_variant_calling}\\n" >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "        - parameter: --min-alternate-count\\n" >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "          value: 2\\n"               >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "        - parameter: --min-alternate-fraction\\n" >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "          value: ${params.min_alternate_fraction_for_variant_calling}\\n" >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "        - parameter: --report-genotype-likelihood\\n" >> ${sample_id}_${short_long}_freebayes_provenance.yml
    printf -- "          value: null\\n"            >> ${sample_id}_${short_long}_freebayes_provenance.yml

    samtools faidx ${ref}

    freebayes \
	--fasta-reference ${ref} \
	--bam ${alignment[0]} \
	--ploidy 1 \
	--min-base-quality ${params.min_base_qual_for_variant_calling} \
	--min-mapping-quality ${params.min_mapping_qual_for_variant_calling} \
	--min-alternate-count 2 \
	--min-alternate-fraction ${params.min_alternate_fraction_for_variant_calling} \
	--report-genotype-likelihood-max \
	> ${sample_id}_${short_long}_freebayes.vcf
    """
}


process plot_coverage {

    tag { sample_id + ' / ' + short_long }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${short_long}_coverage.png"

    input:
    tuple val(sample_id), val(short_long), path(depths), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}_${short_long}_coverage.png"), optional: true

    script:
    """
    plot-coverage.py \
	--sample-id ${sample_id} \
	--ref ${ref} \
	--depths ${depths} \
	--threshold ${params.min_depth} \
	--output ${sample_id}_${short_long}_coverage.png
	
    """
}
