#!/usr/bin/env python3
import argparse


def load_samples(ind_path: str):
    samples = []
    with open(ind_path, 'r') as fin:
        for line in fin:
            samples.append(line.rstrip('\n').split('\t')[0])
    return samples


def main() -> None:
    parser = argparse.ArgumentParser(description='Split Eigenstrat .geno into admixed and reference panels.')
    parser.add_argument('--eigen-prefix', required=True, help='Input Eigenstrat prefix (.snp/.ind/.geno).')
    parser.add_argument('--out-prefix', required=True, help='Output prefix for split panel files.')
    parser.add_argument('--admixed-csv', required=True, help='Comma-separated admixed sample IDs.')
    parser.add_argument('--reference-csv', required=True, help='Comma-separated reference sample IDs.')
    args = parser.parse_args()

    admixed_samples = [s for s in args.admixed_csv.split(',') if s]
    reference_samples = [s for s in args.reference_csv.split(',') if s]

    ind_path = args.eigen_prefix + '.ind'
    snp_path = args.eigen_prefix + '.snp'
    geno_path = args.eigen_prefix + '.geno'

    adm_ind_path = args.out_prefix + '.admixed.ind'
    adm_geno_path = args.out_prefix + '.admixed.geno'
    ref_geno_path = args.out_prefix + '.reference.geno'
    shared_snp_path = args.out_prefix + '.shared.snp'

    samples = load_samples(ind_path)
    sample_to_idx = {s: i for i, s in enumerate(samples)}

    for sample in admixed_samples + reference_samples:
        if sample not in sample_to_idx:
            raise SystemExit(f'Sample not found in .ind: {sample}')

    adm_ix = [sample_to_idx[s] for s in admixed_samples]
    ref_ix = [sample_to_idx[s] for s in reference_samples]

    adm_haps = []
    for i in adm_ix:
        adm_haps.extend([2 * i, 2 * i + 1])
    ref_haps = []
    for i in ref_ix:
        ref_haps.extend([2 * i, 2 * i + 1])

    with open(adm_ind_path, 'w') as fout:
        for sample in admixed_samples:
            fout.write(f'{sample}\tU\tADMIXED\n')

    with open(snp_path, 'r') as fin, open(shared_snp_path, 'w') as fout:
        for line in fin:
            fout.write(line)

    with open(geno_path, 'r') as fgeno, open(adm_geno_path, 'w') as fadm, open(ref_geno_path, 'w') as fref:
        for line_no, line in enumerate(fgeno, start=1):
            g = line.rstrip('\n')
            if line_no == 1:
                expected = 2 * len(samples)
                if len(g) != expected:
                    raise SystemExit(
                        f'Unexpected .geno line length {len(g)} at line 1; expected {expected}. '
                        'This script assumes diploid phased calls for all samples.'
                    )
            fadm.write(''.join(g[i] for i in adm_haps) + '\n')
            fref.write(''.join(g[i] for i in ref_haps) + '\n')

    print('Wrote:')
    print(f'  {shared_snp_path}')
    print(f'  {adm_ind_path}')
    print(f'  {adm_geno_path}')
    print(f'  {ref_geno_path}')


if __name__ == '__main__':
    main()
