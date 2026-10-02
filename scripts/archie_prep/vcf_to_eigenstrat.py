#!/usr/bin/env python3
import argparse
import gzip


def open_maybe_gz(path: str):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def gt_to_haps(gt_field: str) -> str:
    gt = gt_field.split(':')[0].replace('/', '|')
    if gt in {'.', '.|.'}:
        return '99'
    if '|' not in gt:
        if gt == '0':
            return '0'
        if gt == '1':
            return '1'
        return '9'
    a, b = gt.split('|')
    if a == '.' or b == '.':
        return '99'
    if a in {'0', '1'} and b in {'0', '1'}:
        return a + b
    return '99'


def main() -> None:
    parser = argparse.ArgumentParser(description='Convert phased VCF to phased Eigenstrat files.')
    parser.add_argument('--vcf', required=True, help='Input VCF path (.vcf or .vcf.gz).')
    parser.add_argument('--out-prefix', required=True, help='Output prefix for .snp/.ind/.geno files.')
    args = parser.parse_args()

    snp_path = args.out_prefix + '.snp'
    ind_path = args.out_prefix + '.ind'
    geno_path = args.out_prefix + '.geno'

    with open_maybe_gz(args.vcf) as fin, open(snp_path, 'w') as fsnp, open(ind_path, 'w') as find, open(geno_path, 'w') as fgeno:
        written = 0
        removed = 0
        for line in fin:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                parts = line.rstrip('\n').split('\t')
                for sample in parts[9:]:
                    find.write(f'{sample}\tU\tPOP\n')
                continue

            parts = line.rstrip('\n').split('\t')
            chrom, pos, vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            if ',' in alt or len(ref) != 1 or len(alt) != 1:
                removed += 1
                continue
            if vid == '.':
                vid = f'{chrom}:{pos}'
            fsnp.write(f'{vid}\t{chrom}\t0.0\t{pos}\t{ref}\t{alt}\n')
            fgeno.write(''.join(gt_to_haps(gt) for gt in parts[9:]) + '\n')
            written += 1

    print(f'Wrote {written} SNPs to {args.out_prefix}.[snp/ind/geno], removed {removed} filtered variants.')


if __name__ == '__main__':
    main()
