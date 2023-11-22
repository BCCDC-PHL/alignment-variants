#!/usr/bin/env python3

import argparse
import csv
import json
import os
import sys


def parse_depths(depths_file_path):
    depths = {}
    with open(depths_file_path, 'r') as depths_file:
        reader = csv.DictReader(depths_file, delimiter='\t')
        for row in reader:
            chrom = row['chrom']
            if chrom not in depths:
                depths[chrom] = {}
                depths[chrom]['depths'] = []
                depths[chrom]['num_positions'] = 0
            depths[chrom]['depths'].append(int(row['depth']))
            depths[chrom]['num_positions'] += 1

    return depths

def main(args):
    depths_by_contig_id = parse_depths(args.input)

    output_rows = []
    for contig_id, contig_depths in depths_by_contig_id.items():
        for depth_threshold in range(1, args.max_depth + 1):
            output = {
                'contig_id': contig_id,
                'depth_threshold': depth_threshold,
                'num_positions_above_threshold': 0,
                'num_positions': contig_depths['num_positions'],
            }
            for depth in contig_depths['depths']:
                if depth >= depth_threshold:
                    output['num_positions_above_threshold'] += 1
            output['percent_positions_above_threshold'] = round(output['num_positions_above_threshold'] / output['num_positions'] * 100, 2)
            output_rows.append(output)

    output_fieldnames = [
    ]
    if args.sample_id is not None:
        output_fieldnames.append('sample_id')
        
    output_fieldnames += [
        'contig_id',
        'depth_threshold',
        'num_positions_above_threshold',
        'num_positions',
        'percent_positions_above_threshold'
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, delimiter=',', extrasaction='ignore')
    writer.writeheader()
    for output_row in output_rows:
        output_row['sample_id'] = args.sample_id
        writer.writerow(output_row)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate genome fraction coverage')
    parser.add_argument('-s', '--sample-id', help='Sample ID')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    max_depth = parser.add_argument('-d', '--max-depth', type=int, default=100, help='Maximum depth', required=False)
    args = parser.parse_args()
    main(args)
