#!/usr/bin/env python3
import argparse
import os


def main() -> None:
    parser = argparse.ArgumentParser(description='Split shared.snp and panel .geno files by chromosome.')
    parser.add_argument('--split-prefix', required=True, help='Prefix with .shared.snp/.admixed.geno/.reference.geno.')
    parser.add_argument('--chr-dir', required=True, help='Output directory for per-chromosome files.')
    args = parser.parse_args()

    os.makedirs(args.chr_dir, exist_ok=True)

    snp_in = args.split_prefix + '.shared.snp'
    adm_in = args.split_prefix + '.admixed.geno'
    ref_in = args.split_prefix + '.reference.geno'

    handles = {}
    chroms = set()

    def get_handles(chrom):
        if chrom not in handles:
            snp_h = open(f'{args.chr_dir}/{chrom}.snp', 'w')
            adm_h = open(f'{args.chr_dir}/{chrom}.admixed.geno', 'w')
            ref_h = open(f'{args.chr_dir}/{chrom}.reference.geno', 'w')
            handles[chrom] = (snp_h, adm_h, ref_h)
        return handles[chrom]

    with open(snp_in, 'r') as fsnp, open(adm_in, 'r') as fadm, open(ref_in, 'r') as fref:
        for snp_line, adm_line, ref_line in zip(fsnp, fadm, fref):
            chrom = snp_line.split('\t')[1]
            snp_h, adm_h, ref_h = get_handles(chrom)
            snp_h.write(snp_line)
            adm_h.write(adm_line)
            ref_h.write(ref_line)
            chroms.add(chrom)

    for snp_h, adm_h, ref_h in handles.values():
        snp_h.close()
        adm_h.close()
        ref_h.close()

    print(f'Wrote per-chromosome panel files for {len(chroms)} chromosomes in {args.chr_dir}')


if __name__ == '__main__':
    main()
