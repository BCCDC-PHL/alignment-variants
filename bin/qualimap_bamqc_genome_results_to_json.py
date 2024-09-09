#!/usr/bin/env python

import argparse
import json
import os
import sys


def get_section_lines(report_path, section_name):
    """
    """
    section_lines = []

    with open(report_path, 'r') as f:
        # skip lines until we reach the input section
        # indicated by the line starting with '>>>>>>> Input'
        for line in f:
            if line.startswith(f'>>>>>>> {section_name}'):
                break
        for line in f:
            if line.startswith('>>>>>>>'):
                break
            if line.strip() != '':
                section_lines.append(line.strip())
        
    return section_lines


def parse_input_section(report_path):
    """
    """
    section_lines = get_section_lines(report_path, 'Input')
    input_section = {}
    for line in section_lines:
        key, value = line.split('=')
        key = key.strip().replace(' ', '_').lower()
        value = value.strip()
        input_section[key] = value

    return input_section


def parse_reference_section(report_path):
    """
    """
    section_lines = get_section_lines(report_path, 'Reference')
    reference_section = {}
    keys_with_units = [
        'number_of_bases',
    ]
    int_keys = [
        'number_of_bases',
        'number_of_contigs',
    ]
    for line in section_lines:
        key, value = line.split('=')
        key = key.strip().replace(' ', '_').lower()
        value = value.strip()
        if key in keys_with_units:
            value, unit = value.split()
        if key in int_keys:
            value = value.replace(',', '')
            value = int(value)
        reference_section[key] = value

    return reference_section


def parse_int_with_percentage(value):
    """
    Convert a string like "20,096,173 (99.89%)" to a tuple (20096173, 99.89)
    """
    value, percentage = value.split()
    value = value.replace(',', '')
    value = int(value)
    percentage = percentage.strip('()').replace('%', '')
    percentage = float(percentage)

    return value, percentage

def parse_globals_section(report_path):
    """
    """
    section_lines = get_section_lines(report_path, 'Globals')
    globals_section = {}
    int_keys = [
        'number_of_windows',
        'number_of_reads',
        'number_of_secondary_alignments',
        'number_of_mapped_paired_reads_first_in_pair',
        'number_of_mapped_paired_reads_second_in_pair',
        'number_of_mapped_paired_reads_both_in_pair',
        'number_of_mapped_paired_reads_singletons',
        'number_of_overlapping_read_pairs',
        'number_of_duplicated_reads_flagged',
    ]
    int_with_unit_keys = [
        'number_of_mapped_bases',
        'number_of_sequenced_bases',
        'number_of_aligned_bases',
    ]
    int_with_percentage_keys = [
        'number_of_mapped_reads',
        'number_of_supplementary_alignments',
    ]
    for line in section_lines:
        key, value = line.split('=')
        key = key.strip().replace(' ', '_').replace('(', '').replace(')', '').lower()
        value = value.strip()
        if key in int_keys:
            value = value.replace(',', '')
            value = int(value)
        elif key in int_with_percentage_keys:
            value, percentage = parse_int_with_percentage(value)
            percentage_key = key.replace('number_of_', 'percent_')
            globals_section[percentage_key] = percentage
        elif key in int_with_unit_keys:
            value, unit = value.split()
            value = value.replace(',', '')
            value = int(value)
        globals_section[key] = value

    return globals_section


def parse_insert_size_section(report_path):
    """
    """
    insert_size_section = {}
    section_lines = get_section_lines(report_path, 'Insert size')
    float_keys = [
        'mean_insert_size',
        'std_insert_size',   
    ]
    int_keys = [
        'median_insert_size',
    ]
    for line in section_lines:
        key, value = line.split('=')
        key = key.strip().replace(' ', '_').lower()
        value = value.strip()
        if key in float_keys:
            value = value.replace(',', '')
            value = float(value)
        elif key in int_keys:
            value = value.replace(',', '')
            value = int(value)
        insert_size_section[key] = value

    return insert_size_section


def parse_actg_content_section(report_path):
    """
    """
    actg_content_section = {}
    section_lines = get_section_lines(report_path, 'ACTG content')
    int_with_units_and_percentage_keys = [
        'number_of_a_bases',
        'number_of_c_bases',
        'number_of_g_bases',
        'number_of_t_bases',
        'number_of_n_bases',
    ]
    percentage_keys = [
        'gc_percentage',
    ]
    for line in section_lines:
        key, value = line.split('=')
        key = key.strip().replace(' ', '_').replace("'", "_base").lower()
        value = value.strip()
        if key in int_with_units_and_percentage_keys:
            value_split = value.split()
            value_without_units = ' '.join([value_split[0], value_split[2]])
            value, percentage = parse_int_with_percentage(value_without_units)
            percentage_key = key.replace('number_of_', 'percent_')
            actg_content_section[percentage_key] = percentage
        elif key in percentage_keys:
            value = value.replace('%', '')
            value = float(value)
        actg_content_section[key] = value

    return actg_content_section


def parse_mismatches_and_indels_section(report_path):
    """
    """
    mismatches_and_indels_section = {}
    section_lines = get_section_lines(report_path, 'Mismatches and indels')
    int_keys = [
        'number_of_mismatches',
        'number_of_insertions',
        'number_of_deletions',
    ]
    float_keys = [
        'general_error_rate',
    ]
    percentage_keys = [
        'mapped_reads_with_insertion_percentage',
        'mapped_reads_with_deletion_percentage',
        'homopolymer_indels',
    ]
    for line in section_lines:
        key, value = line.split('=')
        key = key.strip().replace(' ', '_').lower()
        value = value.strip()
        if key in int_keys:
            value = value.replace(',', '')
            value = int(value)
        elif key in float_keys:
            value = float(value)
        elif key in percentage_keys:
            value = value.replace('%', '')
            value = float(value)
            if key == 'homopolymer_indels':
                key = 'homopolymer_indels_percentage'
        mismatches_and_indels_section[key] = value

    return mismatches_and_indels_section

def parse_percent_coverage_line(line):
    """
    Convert lines that look like:

    There is a 100% of reference with a coverageData >= 1X

    ...into a tuple (1, 100)
    """
    line_split = line.split()
    coverage = line_split[-1]
    percentage = line_split[3]
    coverage = int(coverage.strip('X'))
    percentage = float(percentage.strip('%'))

    return coverage, percentage

def parse_coverage_section(report_path):
    """
    """
    coverage_section = {
        'coverage_profile': [],
    }
    section_lines = get_section_lines(report_path, 'Coverage')
    float_x_keys = [
        'mean_coverage',
        'std_coverage',
        'paired-end_adapted_mean_coverage',
    ]
    for line in section_lines:
        if line.startswith('There is a'):
            coverage, value = parse_percent_coverage_line(line)
            coverage_profile_record = {
                'depth_threshold': coverage,
                'percent_reference_covered': value,
            }
            coverage_section['coverage_profile'].append(coverage_profile_record)
        else:
            key, value = line.split('=')
            key = key.strip().replace(' ', '_').replace('coverageData', 'coverage').lower()
            value = value.strip()
            if key in float_x_keys:
                value = value.replace('X', '')
                value = float(value)
        coverage_section[key] = value

    return coverage_section


def parse_coverage_per_contig_section(report_path):
    """
    """
    coverage_per_contig_section = []
    section_lines = get_section_lines(report_path, 'Coverage per contig')
    section_headers = [
        'name',
        'length',
        'mapped_bases',
        'mean_coverage',
        'standard_deviation',
    ]
    int_keys = [
        'length',
        'mapped_bases',
    ]
    float_keys = [
        'mean_coverage',
        'standard_deviation',
    ]
    for line in section_lines:
        line_split = line.split('\t')
        contig_data = {}
        for i, header in enumerate(section_headers):
            value = line_split[i]
            if header in int_keys:
                value = int(value)
            elif header in float_keys:
                value = float(value)
            contig_data[header] = value
        coverage_per_contig_section.append(contig_data)

    return coverage_per_contig_section


def main(args):

    parsed_genome_results = {}
    input_section = parse_input_section(args.input)
    parsed_genome_results['input'] = input_section

    reference_section = parse_reference_section(args.input)
    parsed_genome_results['reference'] = reference_section

    globals_section = parse_globals_section(args.input)
    parsed_genome_results['globals'] = globals_section

    insert_size_section = parse_insert_size_section(args.input)
    parsed_genome_results['insert_size'] = insert_size_section

    actg_content_section = parse_actg_content_section(args.input)
    parsed_genome_results['actg_content'] = actg_content_section

    mismatches_and_indels_section = parse_mismatches_and_indels_section(args.input)
    parsed_genome_results['mismatches_and_indels'] = mismatches_and_indels_section

    coverage_section = parse_coverage_section(args.input)
    parsed_genome_results['coverage'] = coverage_section

    coverage_per_contig_section = parse_coverage_per_contig_section(args.input)
    parsed_genome_results['coverage_per_contig'] = coverage_per_contig_section

    print(json.dumps(parsed_genome_results, indent=4))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Description of your script')
    parser.add_argument('--input', type=str, required=True, help='Path to the input file')
    args = parser.parse_args()
    main(args)
    
