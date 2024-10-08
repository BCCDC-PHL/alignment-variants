#!/usr/bin/env python

import argparse
import csv
import json
import sys


def parse_qualimap_bamqc_genome_results(qualimap_bamqc_genome_results):
    """
    Parse the Qualimap BAMQC genome results file.

    :param qualimap_bamqc_genome_results: Path to the Qualimap BAMQC genome results file
    :type qualimap_bamqc_genome_results: str
    :return: The parsed Qualimap BAMQC genome results
    :rtype: dict
    """
    qualimap_bamqc_genome_results_data = {}
    with open(qualimap_bamqc_genome_results, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('number of reads'):
                number_of_reads = line.split('=')[1].strip().replace(',', '')
                qualimap_bamqc_genome_results_data['num_reads'] = int(number_of_reads)
            if line.startswith('number of mapped reads'):
                num_mapped_reads = line.split('=')[1].strip().split(' ')[0].replace(',', '')
                qualimap_bamqc_genome_results_data['num_mapped_reads'] = int(num_mapped_reads)
                percent_mapped_reads = line.split('=')[1].strip().split(' ')[1].strip().replace('(', '').replace(')', '').replace('%', '')
                qualimap_bamqc_genome_results_data['percent_mapped_reads'] = round(float(percent_mapped_reads), 6)
            if line.startswith('number of supplementary alignments'):
                num_supplementary_alignments = int(line.split('=')[1].strip().split()[0].replace(',', ''))
                qualimap_bamqc_genome_results_data['num_supplementary_alignments'] = num_supplementary_alignments
                percent_supplementary_alignments = line.split('=')[1].strip().split(' ')[1].strip().replace('(', '').replace(')', '').replace('%', '')
                qualimap_bamqc_genome_results_data['percent_supplementary_alignments'] = round(float(percent_supplementary_alignments), 6)
            if line.startswith('number of secondary alignments'):
                num_secondary_alignments = int(line.split('=')[1].strip().replace(',', ''))
                qualimap_bamqc_genome_results_data['num_secondary_alignments'] = num_secondary_alignments
            if line.startswith('number of mapped bases'):
                num_mapped_bases = line.split('=')[1].strip().split()[0].replace(',', '')
                qualimap_bamqc_genome_results_data['num_mapped_bases'] = int(num_mapped_bases)
            if line.startswith('number of sequenced bases'):
                num_sequenced_bases = line.split('=')[1].strip().split()[0].replace(',', '')
                qualimap_bamqc_genome_results_data['num_sequenced_bases'] = int(num_sequenced_bases)
            if line.startswith('number of duplicated reads'):
                num_duplicated_reads = line.split('=')[1].strip().replace(',', '')
                qualimap_bamqc_genome_results_data['num_duplicated_reads'] = int(num_duplicated_reads)
                num_mapped_reads = qualimap_bamqc_genome_results_data['num_mapped_reads']
                duplication_rate_percent = (int(num_duplicated_reads) / int(num_mapped_reads)) * 100
                qualimap_bamqc_genome_results_data['duplication_rate_percent'] = round(duplication_rate_percent, 6)
            if line.startswith('mean coverageData'):
                mean_coverage = line.split('=')[1].strip().strip('X').replace(',', '')
                qualimap_bamqc_genome_results_data['mean_depth_coverage'] = round(float(mean_coverage), 6)
            if line.startswith('std coverageData'):
                stdev_coverage = line.split('=')[1].strip().strip('X').replace(',', '')
                qualimap_bamqc_genome_results_data['stdev_depth_coverage'] = round(float(stdev_coverage), 6)
            if line.startswith('mean mapping quality'):
                mean_mapping_quality = line.split('=')[1].strip()
                qualimap_bamqc_genome_results_data['mean_mapping_quality'] = round(float(mean_mapping_quality), 6)
            if line.startswith('general error rate'):
                general_error_rate = line.split('=')[1].strip()
                qualimap_bamqc_genome_results_data['error_rate'] = round(float(general_error_rate), 6)
            if line.startswith('number of mismatches'):
                number_of_mismatches = line.split('=')[1].strip().replace(',', '')
                qualimap_bamqc_genome_results_data['num_mismatches'] = int(number_of_mismatches)
            if line.startswith('number of insertions'):
                number_of_insertions = line.split('=')[1].strip().replace(',', '')
                qualimap_bamqc_genome_results_data['num_insertions'] = int(number_of_insertions)
            if line.startswith('mapped reads with insertion percentage'):
                mapped_reads_with_insertion_percentage = line.split('=')[1].strip().replace('%', '')
                qualimap_bamqc_genome_results_data['mapped_reads_with_insertion_percent'] = round(float(mapped_reads_with_insertion_percentage), 6)
            if line.startswith('number of deletions'):
                number_of_deletions = line.split('=')[1].strip().replace(',', '')
                qualimap_bamqc_genome_results_data['num_deletions'] = int(number_of_deletions)
            if line.startswith('mapped reads with deletion percentage'):
                mapped_reads_with_deletion_percentage = line.split('=')[1].strip().replace('%', '')
                qualimap_bamqc_genome_results_data['mapped_reads_with_deletion_percent'] = round(float(mapped_reads_with_deletion_percentage), 6)
            if 'reference with a coverageData >= 5X' in line:
                proportion_genome_covered_over_5x = float(line.split(' ')[3].strip('%')) / 100
                qualimap_bamqc_genome_results_data['proportion_genome_covered_over_5x'] = round(proportion_genome_covered_over_5x, 6)
            if 'reference with a coverageData >= 10X' in line:
                proportion_genome_covered_over_10x = float(line.split(' ')[3].strip('%')) / 100
                qualimap_bamqc_genome_results_data['proportion_genome_covered_over_10x'] = round(proportion_genome_covered_over_10x, 6)
            if 'reference with a coverageData >= 20X' in line:
                proportion_genome_covered_over_20x = float(line.split(' ')[3].strip('%')) / 100
                qualimap_bamqc_genome_results_data['proportion_genome_covered_over_20x'] = round(proportion_genome_covered_over_20x, 6)
            if 'reference with a coverageData >= 30X' in line:
                proportion_genome_covered_over_30x = float(line.split(' ')[3].strip('%')) / 100
                qualimap_bamqc_genome_results_data['proportion_genome_covered_over_30x'] = round(proportion_genome_covered_over_30x, 6)
            if 'reference with a coverageData >= 40X' in line:
                proportion_genome_covered_over_40x = float(line.split(' ')[3].strip('%')) / 100
                qualimap_bamqc_genome_results_data['proportion_genome_covered_over_40x'] = round(proportion_genome_covered_over_40x, 6)
            if 'reference with a coverageData >= 50X' in line:
                proportion_genome_covered_over_50x = float(line.split(' ')[3].strip('%')) / 100
                qualimap_bamqc_genome_results_data['proportion_genome_covered_over_50x'] = round(proportion_genome_covered_over_50x, 6)

    return qualimap_bamqc_genome_results_data


def main(args):

    parsed_qualimap_bamqc_genome_results_fieldnames = [
        'mean_depth_coverage',
        'stdev_depth_coverage',
        'num_reads',
        'num_mapped_reads',
        'percent_mapped_reads',
        'mean_mapping_quality',
        'num_sequenced_bases',
        'num_mapped_bases',
        'num_mismatches',
        'num_insertions',
        'num_deletions',
        'error_rate',
        'mapped_reads_with_insertion_percent',
        'mapped_reads_with_deletion_percent',
        'num_secondary_alignments',
        'num_supplementary_alignments',
        'percent_supplementary_alignments',
        'num_duplicated_reads',
        'duplication_rate_percent',
        'proportion_genome_covered_over_5x',
        'proportion_genome_covered_over_10x',
        'proportion_genome_covered_over_20x',
        'proportion_genome_covered_over_30x',
        'proportion_genome_covered_over_40x',
        'proportion_genome_covered_over_50x',
    ]

    output_data = {
        'sample_id': args.sample_id,
        'read_type': args.read_type,
    }
    
    if args.failed:
        for field in parsed_qualimap_bamqc_genome_results_fieldnames:
            output_data[field] = None
    else:
        qualimap_bamqc_genome_results_data = parse_qualimap_bamqc_genome_results(args.qualimap_bamqc_genome_results)
        output_data = {
            'sample_id': args.sample_id,
            'read_type': args.read_type,
        }
        for field in parsed_qualimap_bamqc_genome_results_fieldnames:
            output_data[field] = qualimap_bamqc_genome_results_data.get(field, None)

    output_fieldnames = ['sample_id', 'read_type'] + parsed_qualimap_bamqc_genome_results_fieldnames

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='unix', extrasaction='ignore', quoting=csv.QUOTE_MINIMAL)
    writer.writeheader()
    writer.writerow(output_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--qualimap-bamqc-genome-results', type=str, help='Path to the Qualimap BAMQC genome results file. May be omitted when using the --failed flag to generate a row with all stats set to None (blank)')
    parser.add_argument('-s', '--sample-id', type=str, help='Sample ID')
    parser.add_argument('-t', '--read-type', choices=['short', 'long'], default='short', help='Read type')
    parser.add_argument('-f', '--failed', action='store_true', help='Flag to indicate that the sample failed QC and no input data was provided')
    args = parser.parse_args()
    main(args)
