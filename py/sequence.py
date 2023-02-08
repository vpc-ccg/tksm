#!/usr/bin/env python3

import uuid
import random
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
                        help="FASTQ file.")
    parser.add_argument("--perfect-reads",
                        type=str,
                        help="Generate perfect reads.")
    parser.add_argument('-t',
                        '--threads',
                        type=int,
                        default=1,
                        help="Number of threads.",
                        )
    args = parser.parse_args()
    if not args.output and not args.perfect_reads:
        parser.error("Must specify either --output or --perfect-reads.")
    return args

def generate_fasta(fasta):
    name = ""
    seq = list()
    if fasta.endswith('.gz'):
        f = gzip.open(fasta, 'rt')
    else:
        f = open(fasta, 'r')
    for _, l in enumerate(f):
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


def badread_fastq(molecule_id, seq):
    read_id = uuid.UUID(int=random.getrandbits(128))
    result = list()
    info = f"length={len(seq)} molecule_id={molecule_id}"
    result.append(f"@{read_id} {info}")
    result.append(seq)
    result.append("+")
    result.append("K" * len(seq))
    return "\n".join(result)

def perfect_fastq(molecule_id, seq):
    read_id = uuid.UUID(int=random.getrandbits(128))
    result = list()
    info = f"length={len(seq)} molecule_id={molecule_id}"
    result.append(f"@{read_id} {info}")
    result.append(seq)
    result.append("+")
    result.append("K" * len(seq))
    return "\n".join(result)

def mdf_to_seq(mdf, targets=dict()):
    molecule_id, intervals = mdf
    seq = list()
    for chrom, start, end, strand, modifications in intervals:
        segment = reference_seqs[chrom][start:end].upper()
        segment = apply_modifications(segment, modifications)
        if strand == '+':
            seq.append(segment)
        else:
            seq.append(reverse_complement(segment))
    seq = ''.join(seq)

    results = dict()
    for k,v in targets.items():
        results[k] = v(molecule_id, seq)
    return results

if __name__ == "__main__":
    args = parse_args()
    reference_seqs = get_reference_seqs(args.reference)

    targets = dict()
    target_outfiles = dict()
    if args.output:
        targets['badread_fastq']=badread_fastq
        target_outfiles['badread_fastq'] = open(args.output, 'w+')
    if args.perfect_reads:
        targets['perfect_fastq']=perfect_fastq
        target_outfiles['perfect_fastq'] = open(args.perfect_reads, 'w+')
    mdf_to_seq_targeted = functools.partial(mdf_to_seq, targets=targets)
    with open(args.input, 'r') as f:
        mdg = mdf_generator(f)
        if args.threads > 1:
            mapper = functools.partial(Pool(args.threads).imap_unordered, chunksize=1000)
        else:
            mapper = map
        for read_dict in tqdm(mapper(mdf_to_seq_targeted, mdg), desc="Sequencing"):
            for k,v in read_dict.items():
                print(v, file=target_outfiles[k])

        if args.threads > 1:
            Pool(args.threads).close()

    for v in target_outfiles.values():
        v.close()
