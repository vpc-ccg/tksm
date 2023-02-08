#!/usr/bin/env python3

import sys
import uuid
import random
import argparse
import functools
from multiprocessing import Pool
import gzip

from badread_scripts.error_model import ErrorModel
from badread_scripts.qscore_model import QScoreModel
from badread_scripts.tail_noise_model import KDE_noise_generator
from badread_scripts.identities import Identities
from badread_scripts import settings
from badread_scripts.simulate import sequence_fragment

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
                        "--badread-fastq",
                        type=str,
                        help="FASTQ file.")
    parser.add_argument("--perfect-fastq",
                        type=str,
                        help="Generate perfect reads.")
    parser.add_argument('-t',
                        '--threads',
                        type=int,
                        default=1,
                        help="Number of threads.",
                        )

    parser.add_argument('--badread-identity', type=str, default='87.5,97.5,5',
                          help='Sequencing identity distribution (mean, max and stdev, '
                               'default: DEFAULT)')
    parser.add_argument('--badread-error_model', type=str, default='nanopore2020',
                          help='Can be "nanopore2018", "nanopore2020", "pacbio2016", "random" or '
                               'a model filename')
    parser.add_argument('--badread-qscore_model', type=str, default='nanopore2020',
                          help='Can be "nanopore2018", "nanopore2020", "pacbio2016", "random", '
                               '"ideal" or a model filename')
    parser.add_argument('--badread-tail_noise_model', type=str, default='none',
                          help='Can be "nanopore", "pacbio", "none" or a model filename')

    args = parser.parse_args()

    # Process arguments and check for errors
    try:
        identity_parameters = [float(x) for x in args.badread_identity.split(',')]
        assert len(identity_parameters) == 3, "Must specify 3 values for --badread-identity"
        args.badread_mean_identity = identity_parameters[0]
        args.badread_max_identity = identity_parameters[1]
        args.badread_identity_stdev = identity_parameters[2]
    except (ValueError, IndexError):
        sys.exit('Error: could not parse --identity values')
    if args.badread_mean_identity > 100.0:
        sys.exit('Error: mean read identity cannot be more than 100')
    if args.badread_max_identity > 100.0:
        sys.exit('Error: max read identity cannot be more than 100')
    if args.badread_mean_identity <= settings.MIN_MEAN_READ_IDENTITY:
        sys.exit(f'Error: mean read identity must be at least {settings.MIN_MEAN_READ_IDENTITY}')
    if args.badread_max_identity <= settings.MIN_MEAN_READ_IDENTITY:
        sys.exit(f'Error: max read identity must be at least {settings.MIN_MEAN_READ_IDENTITY}')
    if args.badread_mean_identity > args.badread_max_identity:
        sys.exit(f'Error: mean identity ({args.badread_mean_identity}) cannot be larger than max '
                 f'identity ({args.badread_max_identity})')
    if args.badread_identity_stdev < 0.0:
        sys.exit('Error: read identity stdev cannot be negative')

    if not args.badread_fastq and not args.perfect_fastq:
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


def badread_fastq(molecule_id, raw_seq):
    target_identity = identities.get_identity()
    seq, quals, actual_identity, identity_by_qscores = sequence_fragment(raw_seq, target_identity, error_model, qscore_model, tail_model)
    read_id = uuid.UUID(int=random.getrandbits(128))
    result = list()
    info = f"length={len(seq)} error_free_length={len(raw_seq)} read_identity={actual_identity * 100.0:.2f}% molecule_id={molecule_id}"
    result.append(f"@{read_id} {info}")
    result.append(seq)
    result.append("+")
    result.append(quals)
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
    if args.badread_fastq:
        identities = Identities(args.badread_mean_identity, args.badread_identity_stdev, args.badread_max_identity, sys.stderr)
        error_model = ErrorModel(args.badread_error_model, sys.stderr)
        qscore_model = QScoreModel(args.badread_qscore_model, sys.stderr)
        tail_model = KDE_noise_generator.load(args.badread_tail_noise_model)

        targets['badread_fastq']=badread_fastq
        target_outfiles['badread_fastq'] = open(args.badread_fastq, 'w+')
    if args.perfect_fastq:
        targets['perfect_fastq']=perfect_fastq
        target_outfiles['perfect_fastq'] = open(args.perfect_fastq, 'w+')

    mdf_to_seq_targeted = functools.partial(mdf_to_seq, targets=targets)
    with open(args.input, 'r') as f:
        mdg = mdf_generator(f)
        if args.threads > 1:
            mapper = functools.partial(Pool(args.threads).imap_unordered, chunksize=10)
        else:
            mapper = map
        for read_dict in tqdm(mapper(mdf_to_seq_targeted, mdg), desc="Sequencing"):
            for k,v in read_dict.items():
                print(v, file=target_outfiles[k])

        if args.threads > 1:
            Pool(args.threads).close()

    for v in target_outfiles.values():
        v.close()
