#!/usr/bin/env python3

import argparse
import csv
import json
import os
import sys

def parse_samplesheet(samplesheet_path):
    samples = []
    with open(samplesheet_path, 'r') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            samples.append(row)

    return samples
        

def main(args):
    base_dir_abspath = os.path.abspath(args.base_dir)
    samples = parse_samplesheet(args.input)

    for sample in samples:
        for field in ['R1', 'R2', 'LONG', 'REF']:
            if field in sample:
                sample[field] = os.path.join(base_dir_abspath, sample[field])

    output_fields = [
        'ID',
        'R1',
        'R2',
        'LONG',
        'REF',
    ]

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fields, delimiter=',', extrasaction='ignore')
    writer.writeheader()
    for sample in samples:
        writer.writerow(sample)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input file')
    parser.add_argument('-b', '--base-dir', default='.', help='')
    args = parser.parse_args()
    main(args)
