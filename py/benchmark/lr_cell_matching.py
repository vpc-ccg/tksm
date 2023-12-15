#!/usr/bin/env python3
from collections import Counter
import sys
import argparse
from multiprocessing import Pool
import gzip

import edlib
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Extract long-read barcodes using short-read adapter alignment",
    )
    parser.add_argument(
        "-i",
        "--reads",
        type=str,
        nargs="+",
        required=True,
        help="FASTQ/A files of long-read (can be gzipped)",
    )
    parser.add_argument(
        "-b",
        "--short-read-barcodes",
        type=str,
        required=True,
        help="Short-read barcode list TSV file (first column is barcode)",
    )
    parser.add_argument(
        "-t", "--threads", default=1, type=int, help="Number of threads."
    )
    parser.add_argument(
        "-s",
        "--start",
        default=25,
        type=int,
        help="Start of segment. Will do seq[s:e] and seq[-e:-s].",
    )
    parser.add_argument(
        "-e",
        "--end",
        default=70,
        type=int,
        help="End of segment. Will do seq[s:e] and seq[-e:-s].",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        default=None,
        help="Path to output file. Output file is gzipped. STDOUT is in normal text.",
    )
    args = parser.parse_args()
    assert args.threads > 0
    return args


sr_barcodes = set()

rev_compl_l = [chr(i) for i in range(128)]
rev_compl_l[ord("A")] = "T"
rev_compl_l[ord("C")] = "G"
rev_compl_l[ord("G")] = "C"
rev_compl_l[ord("T")] = "A"


def rev_compl(s):
    return "".join(rev_compl_l[ord(c)] for c in reversed(s))


def get_matches(lr_segment):
    matches = list()
    d = float("inf")

    for idx, sb in enumerate(sr_barcodes):
        for lb in lr_segment:
            if len(lb) == 0:
                continue
            aln = edlib.align(sb, lb, "HW", "distance")
            if aln["editDistance"] < d:
                matches = list()
                d = aln["editDistance"]
            if aln["editDistance"] == d:
                matches.append(idx)
    return len(matches), d


def get_lr_segments(reads, start, end, sample_size=50_000):
    lr_segments = list()
    seg_len = end - start
    for f in reads:
        if f.endswith("gz"):
            f = gzip.open(f, "rt")
        else:
            f = open(f, "r")
        mod = -1
        for idx, l in tqdm(enumerate(f), desc=f"Reading {f}"):
            if mod == -1:
                if l.startswith("@"):
                    mod = 2
                elif l.startswith(">"):
                    mod = 4
                else:
                    raise ValueError(f"Unknown read format: {l}")
            if idx % mod == 1:
                seq = l.rstrip()
                if len(seq) < 2 * seg_len:
                    lr_segments.append((seq, ""))
                else:
                    lr_segments.append((seq[start:end], seq[-end:-start]))
    sample_idxs = np.random.choice(
        np.arange(len(lr_segments)), min(sample_size, len(lr_segments)), replace=False
    )
    lr_segments = [lr_segments[i] for i in sample_idxs]
    return lr_segments


def run_get_matches(lr_segments, threads):
    lr_counter = Counter()
    with Pool(threads) as p:
        for match_count, d in tqdm(
            p.imap(get_matches, lr_segments, chunksize=10), total=len(lr_segments)
        ):
            lr_counter[(match_count, d)] += 1
    return lr_counter


def get_sr_barcodes(barcodes_tsv):
    gz = barcodes_tsv.endswith("gz")
    global sr_barcodes
    sr_barcodes = set()
    if gz:
        f = gzip.open(barcodes_tsv)
    else:
        f = open(barcodes_tsv)
    for l in tqdm(f, desc=f"Reading {barcodes_tsv}"):
        if gz:
            l = l.decode()
        l = l.rstrip().split("\t")
        b = l[0]
        sr_barcodes.add(b)
        sr_barcodes.add(rev_compl(b))
    sr_barcodes = list(sr_barcodes)


def main():
    args = parse_args()

    get_sr_barcodes(args.short_read_barcodes)
    lr_segments = get_lr_segments(args.reads, args.start, args.end)
    lr_counter = run_get_matches(lr_segments, args.threads)
    new_lr_counter = Counter()
    for (m, d), v in lr_counter.items():
        if d > 2:
            d = ">2"
        else:
            d = str(d)
        if m > 1:
            m = "ambig"
        else:
            m = "uniq"
        new_lr_counter[(m, d)] += v
    if args.outfile:
        args.outfile = open(args.outfile, "wt")
    else:
        args.outfile = sys.stdout
    args.outfile.write(f"Match\tDist\tCount\t%\n")
    T = sum(new_lr_counter.values())
    for k, v in sorted(new_lr_counter.items()):
        args.outfile.write(f"{k[0]}\t{k[1]}\t{v}\t{v/T:.1%}\n")
    args.outfile.close()


if __name__ == "__main__":
    main()
