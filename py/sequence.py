#!/usr/bin/env python3
import argparse
import functools
from multiprocessing import Pool
import gzip

from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(
        description="Sequencer module of RNAInFuser. Generates FASTQ file from MDF (molecule description file) and reference FASTA files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i",
                        "--input",
                        type=str,
                        required=True,
                        help="MDF file.")
    parser.add_argument("-r",
                        "--reference",
                        type=str,
                        nargs="+",
                        required=True,
                        help="Reference FASTA files, space separated.")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        help="FASTQ file.")
    parser.add_argument('-t',
                        '--threads',
                        type=int,
                        default=1,
                        help="Number of threads.",
                        )
    args = parser.parse_args()
    return args

def generate_fasta(fasta):
    name = ""
    seq = list()
    if fasta.endswith('.gz'):
        f = gzip.open(fasta, 'rt')
    else:
        f = open(fasta, 'r')
    for idx, l in enumerate(f):
        l = l.rstrip('\n')
        if l[0]=='>':
            if len(seq) == 0:
                name = l[1:].split(' ')[0]
                continue
            yield (name, ''.join(seq))
            seq = list()
            name = l[1:].split(' ')[0]
        else:
            seq.append(l)
    yield (name, ''.join(seq))

def get_reference_seqs(reference):
    reference_seqs = dict()
    for ref in reference:
        print(f"Loading reference {ref}...")
        reference_seqs.update({
            name: seq for name, seq in generate_fasta(ref)
        })
    return reference_seqs

def mdf_generator(f):
    # MDF format:
    # Each entry starts with a header: +<mol_id>\t<depth>\t<comments>; 
    #   where <mol_id> is the molecule ID, <depth> is the depth of the molecule, and each <comments> entry is key=value pair terminated by a semicolon.
    # Each line after the header, and until the next header or end of file, is an interval description: <ref_id>\t<start:int>\t<end:int>\t<strand>\t<modifications>
    read_id = None
    intervals = list()

    for line in f:
        line = line.strip('\n').split('\t')
        if line[0][0] == '+':
            if read_id is not None:
                for _ in range(depth):
                    yield read_id, intervals
            read_id = line[0][1:]
            depth = int(line[1])
            intervals = list()
        else:
            chrom, start, end, strand, modifications = line
            start, end = int(start), int(end)
            intervals.append(
                (chrom, start, end, strand, modifications)
            )
    if read_id is not None:
        for _ in range(depth):
            yield read_id, intervals

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def apply_modifications(seq, modifications):
    if modifications == '':
        return seq
    else:
        modifications = modifications.split(',')
        seq = list(seq)
        for mod in modifications:
            char = mod[-1]
            pos = int(mod[:-1])
            seq[pos] = char
        return ''.join(seq)

def apply_badread(seq):
    return seq

def mdf_to_seq(mdf):
    read_id, intervals = mdf
    seq = list()
    for chrom, start, end, strand, modifications in intervals:
        segment = reference_seqs[chrom][start:end].upper()
        segment = apply_modifications(segment, modifications)
        if strand == '+':
            seq.append(segment)
        else:
            seq.append(reverse_complement(segment))
    seq = ''.join(seq)
    seq = apply_badread(seq)

    result = list()
    result.append(f"@{read_id} length={len(seq)}")
    result.append(seq)
    result.append("+")
    result.append("K" * len(seq))
    return '\n'.join(result)

if __name__ == "__main__":
    args = parse_args()
    global reference_seqs
    reference_seqs = get_reference_seqs(args.reference)
    with open(args.input, 'r') as f, open(args.output, 'w+') as o:
        mdg = mdf_generator(f)
        if args.threads > 1:
            mapper = functools.partial(Pool(args.threads).imap_unordered, chunksize=1000)
        else:
            mapper = map
        
        for read in tqdm(mapper(mdf_to_seq, mdg), desc="Sequencing"):
            print(read, file=o)
        if args.threads > 1:
            Pool(args.threads).close()