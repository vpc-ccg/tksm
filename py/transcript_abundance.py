#!/usr/bin/env python3
import argparse
import gzip
from collections import defaultdict
from typing import List, Tuple
import numpy as np

from tqdm import tqdm

# Includes EM code by Jared Simpson from
# https://github.com/jts/nanopore-rna-analysis/blob/master/nanopore_transcript_abundance.py

IUPAC_nts = {
    "A": np.array(["A"], dtype=str),
    "C": np.array(["C"], dtype=str),
    "G": np.array(["G"], dtype=str),
    "T": np.array(["T"], dtype=str),
    "R": np.array(["A", "G"], dtype=str),
    "Y": np.array(["C", "T"], dtype=str),
    "K": np.array(["G", "T"], dtype=str),
    "M": np.array(["A", "C"], dtype=str),
    "S": np.array(["C", "G"], dtype=str),
    "W": np.array(["A", "T"], dtype=str),
    "B": np.array(["C", "G", "T"], dtype=str),
    "D": np.array(["A", "G", "T"], dtype=str),
    "H": np.array(["A", "C", "T"], dtype=str),
    "V": np.array(["A", "C", "G"], dtype=str),
    "N": np.array(["A", "C", "G", "T"], dtype=str),
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Output a TSV file of the long-read trasncript expresion.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-p",
        "--paf",
        type=str,
        required=True,
        help="PAF file of the mapping of the long-reads to the reference transcriptome.",
    )
    parser.add_argument(
        "-m",
        "--lr-br",
        type=str,
        default="",
        help="TSV output of the long-reads barcode matching from scTager (e.g. S1.lr_matches.tsv.gz)",
    )
    parser.add_argument(
        "--cb-count",
        type=int,
        default=0,
        help="Number of cell barcodes to simulate."
        + " If >0, then the abundance will be split by cell barcodes according to a normal distribution."
        + " Cannot be used with --lr-br parameter.",
    )
    parser.add_argument(
        "--cb-lognorm-params",
        type=str,
        default="10,1",
        help="Parameters for the lognormal distribution of cell barcodes abundance. "
        + " Takes two comma-separated values: mean and standard deviation.",
    )
    parser.add_argument(
        "--cb-pattern",
        type=str,
        default="NNNNNNNNNNNN",
        help="FASTA pattern for the cell barcode. If used, --cb-count parameter must also be used.",
    )
    parser.add_argument(
        "--cb-dropout",
        type=float,
        default=0.2,
        help="Fraction of reads without cell barcodes. If used, --cb-count parameter must also be used.",
    )
    parser.add_argument(
        "--cb-txt",
        type=str,
        default="",
        help="Path to a text file of cell barcodes whitelist, one per line."
        + "If used, --cb-count parameter must also be used."
        + "Do not use with --cb-pattern parameter.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path for the output file (e.g cDNA.abundance.tsv.tsv.gz or cDNA.abundance.tsv)."
        + "Will be gzipped if it ends with .gz extenstion.",
    )
    parser.add_argument(
        "-em",
        "--em-iterations",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=42,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        type=int,
        default=0,
    )

    class ListPrinter(argparse.Action):
        def __call__(self, parser, namespace, values, option_string):
            txt = "\n".join([getattr(k, "dest") for k in parser._actions])
            print(txt)
            parser.exit()

    parser.add_argument("--list", nargs=0, action=ListPrinter)

    args = parser.parse_args()
    if args.cb_count > 0:
        if args.lr_br != "":
            parser.error("--lr-br must not be set with --cb-count")
        if args.cb_pattern == "":
            parser.error("--cb-pattern must be set with --cb-count")
        for c in args.cb_pattern:
            if c not in IUPAC_nts.keys():
                parser.error(
                    "--cb-pattern must contain only valid IUPAC nucleotide letters: "
                    + f"<{c}> not in {','.join(IUPAC_nts.keys())}"
                )
        if args.cb_dropout < 0 or args.cb_dropout > 1:
            parser.error("--cb-dropout must be between 0 and 1")
        args.cb_lognorm_params = tuple(
            float(x) for x in args.cb_lognorm_params.split(",")
        )
        assert len(args.cb_lognorm_params) == 2
        assert args.cb_lognorm_params[1] > 0
    return args


def generate_barcodes_from_pattern(pattern: str, count: int):
    barcodes: List[str] = list()
    for _ in range(count):
        barcode: List[str] = list()
        for p in pattern:
            c: str = np.random.choice(IUPAC_nts[p])
            barcode.append(c)
        barcodes.append("".join(barcode))
    return barcodes


def parse_barcodes_txt(barcode_txt: str):
    barcodes: List[str] = list()
    if barcode_txt.endswith(".gz"):
        infile = gzip.open(barcode_txt, "rt")
    else:
        infile = open(barcode_txt, "r")
    for l in infile:
        l: str = l
        barcode = l.rstrip("\n")
        barcodes.append(barcode)
    return barcodes


def parse_lr_bc_matches(lr_br_tsv):
    rid_to_bc: defaultdict[str, str] = defaultdict(lambda: ".")
    if lr_br_tsv.endswith(".gz"):
        infile = gzip.open(lr_br_tsv, "rt")
    else:
        infile = open(lr_br_tsv, "r")
    print("Parsing LR barcode matches TSV...")
    for l in tqdm(infile):  # type: ignore
        l: str = l
        rid, _, c, _, bc = l.rstrip("\n").split("\t")
        if c != "1":
            continue
        rid_to_bc[rid] = bc
    return rid_to_bc


def parse_paf(paf):
    alignments = defaultdict(list)
    tname_to_tid = dict()
    tid_to_tname = dict()
    print("Parsing PAF file...")
    for line in tqdm(open(paf)):
        alignment = dict()
        line = line.rstrip("\n").split("\t")
        rid = line[0]
        alignment["query_length"] = int(line[1])
        tname = line[5]
        if not tname in tname_to_tid:
            tid = len(tname_to_tid)
            tname_to_tid[tname] = tid
            tid_to_tname[tid] = tname
        alignment["tid"] = tname_to_tid[tname]
        alignment["target_start"] = int(line[7])
        alignment["num_matches"] = int(line[9])
        alignment["alignment_block_length"] = int(line[10])

        alignments[rid].append(alignment)
    return tid_to_tname, alignments


# Process the alignment records for a single read and fractionally assign
# that read to the set of transcripts the read is "compatible" with.
# Compatibility is defined by the length of an alignment relative
# to the longest alignment for the read.
def get_compatibility(alignments):
    transcript_compatibility = defaultdict(list)

    full_length_min_distance = 20
    min_read_length = 0
    # All records within threshold of the best score are considered to be compatible with the read
    threshold = 0.95

    def is_full_length(p):
        return p < full_length_min_distance

    print("Computing compatibility of different alignments for each read...")
    for rid, records in tqdm(alignments.items(), total=len(alignments)):
        # Determine best match
        read_length = records[0]["query_length"]
        best_match_align_len = 0
        best_num_matches = 0
        best_is_full_length = False

        for r in records:
            fl = is_full_length(r["target_start"])
            if r["num_matches"] > best_num_matches or (
                r["num_matches"] == best_num_matches and fl
            ):
                best_match_align_len = r["alignment_block_length"]
                best_num_matches = r["num_matches"]
                best_is_full_length = fl

        fraction_aligned = best_match_align_len / float(read_length)
        if fraction_aligned < 0.5 or read_length < min_read_length:
            continue

        def is_equivalent_hit(x):
            f = float(x["num_matches"]) / best_num_matches
            l = is_full_length(x["target_start"])
            return f > threshold and l == best_is_full_length

        # Count equivalent hits
        num_hits = 0
        for r in records:
            if is_equivalent_hit(r):
                num_hits += 1

        for r in records:
            if is_equivalent_hit(r):
                transcript_compatibility[rid].append((r["tid"], 1.0 / num_hits))
    return transcript_compatibility


# Calculate the abundance of the transcript set based on read-transcript compatibilities
def calculate_abundance(compatibility, verbose=0):
    abundance = defaultdict(float)
    total = 0
    for read in compatibility:
        if verbose > 1:
            print(f"[compatibility] {read}: {compatibility[read]}")

        for t in compatibility[read]:
            abundance[t[0]] += t[1]
            total += t[1]

    for transcript in abundance:
        abundance[transcript] = abundance[transcript] / total
        if verbose > 0:
            print(f"[abundance] {transcript}: {abundance[transcript]}")
    return abundance


# Update read-transcript compatibility based on transcript abundances
def update_compatibility(compatibility, abundance):
    for read in compatibility:
        ids = list()
        total = 0
        for t in compatibility[read]:
            total += abundance[t[0]]
            ids.append(t[0])

        compatibility[read] = list()
        for i in ids:
            compatibility[read].append((i, abundance[i] / total))


def calculate_split_abundance(compatibility, rid_to_bc):
    abundance = defaultdict(float)
    total = 0
    for read in compatibility:
        for tid, count in compatibility[read]:
            abundance[(tid, rid_to_bc[read])] += count
            total += count

    for transcript in abundance:
        abundance[transcript] = abundance[transcript] / total
    return abundance


def generate_rid_to_bc(
    barcodes: List[str], dropout: float, lognorm_params: Tuple[float, float]
):
    weights = np.random.lognormal(*lognorm_params, size=len(barcodes))
    total_with_dropout: float = weights.sum() / (1 - dropout)
    dropout_weight: float = total_with_dropout * dropout
    barcodes.append(".")
    weights = np.append(weights, np.array([dropout_weight]))
    weights = weights / weights.sum()

    rid_to_bc: defaultdict[str, str] = defaultdict(
        lambda: np.random.choice(
            a=np.array(barcodes, dtype=str),
            size=1,
            p=weights,
            replace=True,
        )[0]
    )
    return rid_to_bc


def main():
    args = parse_args()
    np.random.seed(args.random_seed)
    if args.lr_br == "":
        if args.cb_count <= 0:
            rid_to_bc: defaultdict[str, str] = defaultdict(lambda: ".")
        else:
            if args.cb_txt != "":
                barcodes = parse_barcodes_txt(args.cb_txt)
                assert len(barcodes) >= args.cb_count
                barcodes = np.random.choice(
                    barcodes,
                    size=args.cb_count,
                    replace=True,
                )
            else:
                barcodes = generate_barcodes_from_pattern(
                    args.cb_pattern,
                    args.cb_count,
                )
            rid_to_bc = generate_rid_to_bc(
                barcodes, args.cb_dropout, args.cb_lognorm_params
            )
    else:
        rid_to_bc = parse_lr_bc_matches(args.lr_br)
    tid_to_tname, alignments = parse_paf(args.paf)
    transcript_compatibility = get_compatibility(alignments)
    del alignments

    print("Running EM...")
    for _ in tqdm(range(args.em_iterations)):
        # Calculate abundance from compatibility assignments
        abundance = calculate_abundance(transcript_compatibility)
        # Update compatibility assignments
        update_compatibility(transcript_compatibility, abundance)

    # Split transcript abundances by cellular barcode
    abundance = calculate_split_abundance(transcript_compatibility, rid_to_bc)
    # Write results as a TSV file
    if args.output.endswith(".gz"):
        outfile = gzip.open(args.output, "wt+")
    else:
        outfile = open(args.output, "w+")
    total_reads = len(transcript_compatibility)
    print(f"Parsed alignments for {total_reads} reads")
    outfile.write("target_id\ttpm\tcell\n")  # type: ignore
    for (tid, cell), a in abundance.items():
        tpm = a * 1_000_000
        # if you need >100M reads to see this transcript, then skip
        if tpm < 0.001 or f"{tpm:.3f}" == "0.000":
            continue
        outfile.write(
            "\t".join(
                [
                    f"{tid_to_tname[tid]}",
                    f"{tpm:.3f}",
                    f"{cell}",
                ]
            )  # type: ignore
        )
        outfile.write("\n")  # type: ignore
    outfile.close()


if __name__ == "__main__":
    main()
