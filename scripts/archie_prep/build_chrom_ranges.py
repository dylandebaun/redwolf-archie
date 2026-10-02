#!/usr/bin/env python3
import argparse
import glob
import os
import re


def natural_key(value: str):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r'(\d+)', value)]


def first_last_pos(snp_path: str):
    first = None
    last = None
    with open(snp_path, 'r') as fin:
        for line in fin:
            pos = line.rstrip('\n').split('\t')[3]
            if first is None:
                first = pos
            last = pos
    if first is None or last is None:
        raise SystemExit(f'No SNP rows found in {snp_path}')
    return first, last


def main() -> None:
    parser = argparse.ArgumentParser(description='Build chromosome begin/end range TSV from per-chromosome .snp files.')
    parser.add_argument('--chr-dir', required=True, help='Directory containing <chrom>.snp files.')
    parser.add_argument('--out-tsv', required=True, help='Output TSV path with columns: chrom begin end.')
    args = parser.parse_args()

    snp_files = sorted(glob.glob(os.path.join(args.chr_dir, '*.snp')), key=natural_key)
    if not snp_files:
        raise SystemExit(f'No .snp files found in {args.chr_dir}')

    with open(args.out_tsv, 'w') as fout:
        for snp_file in snp_files:
            chrom = os.path.basename(snp_file)[:-4]
            begin, end = first_last_pos(snp_file)
            fout.write(f'{chrom}\t{begin}\t{end}\n')

    print(f'Wrote chromosome ranges: {args.out_tsv}')


if __name__ == '__main__':
    main()
