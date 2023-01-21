#!/usr/bin/env python3
import argparse
import gzip
from collections import Counter

from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(
        description="Output a TSV file of the long-read trasncript expresion per cell")
    parser.add_argument("-r",
                        "--ref",
                        type=str,
                        required=True,
                        help="Reference cDNA FASTA file to get the transcript IDs")
    parser.add_argument("-m",
                        "--lr_tsv",
                        type=str,
                        required=True,
                        help="TSV output of the long-reads barcode matching from scTager (e.g. S1.lr_matches.tsv.gz)")
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="PAF file of the mapping of the long-reads to the reference transcriptome.")
    parser.add_argument("-o",
                        "--prefix",
                        type=str,
                        required=True,
                        help="Prefix path for the output files: <prefix>.npy.gz, <prefix>.bc.gz, and <prefix>.tids.gz")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    rid_to_bc = dict()
    print('Reading LR barcode matches...')
    for l in tqdm(gzip.open(args.lr_tsv,'rt')):
        l = l.rstrip('\n').split('\t')
        rid,d,c,s,bc = l
        if c != '1':
            continue
        rid_to_bc[rid] = bc
    
    bc_to_idx = dict()
    bcs = list()
    print('Collecting unique barcodes...')
    for bc in tqdm(rid_to_bc.values(), total=len(rid_to_bc)):
        if not bc in bc_to_idx:
            bcs.append(bc)
            bc_to_idx[bc] = len(bc_to_idx)

    tid_to_idx = dict()
    tids = list()
    print('Reading transcript IDs from cDNA ref...')
    for l in tqdm(open(args.ref)):
        if not l.startswith('>'):
            continue
        tid = l[1:].split()[0].split('.')[0]
        if not tid in tid_to_idx:
            tids.append(tid)
            tid_to_idx[tid] = len(tid_to_idx)

    A = Counter()
    print('Counting expression from PAF file...')
    for l in tqdm(open(args.paf)):
        l = l.rstrip('\n').split('\t')
        rid = l[0]
        if not rid in rid_to_bc or not 'tp:A:P' in l[12:]:
            continue
        bc = rid_to_bc[rid]
        tid = l[5].split('.')[0]
        A[(bc_to_idx[bc],tid_to_idx[tid])] += 1

    print("Writing to NPY matrix...")
    outfile = gzip.open(f'{args.prefix}.npy.gz','wt+')
    for (bc_idx,tid_idx),v in tqdm(A.items()):
        outfile.write(f'{bc_idx}\t{tid_idx}\t{v}\n')
    outfile.close()

    print("Writing transcript IDs...")
    outfile = gzip.open(f'{args.prefix}.tids.gz','wt+')
    for tid in tids:
        outfile.write(f'{tid}\n')
    outfile.close()

    print("Writing barcodes...")
    outfile = gzip.open(f'{args.prefix}.bc.gz','wt+')
    for bc in bcs:
        outfile.write(f'{bc}\n')
    outfile.close()


if __name__ == "__main__":
    main()


