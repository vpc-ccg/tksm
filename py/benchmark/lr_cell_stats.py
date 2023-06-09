import argparse
from collections import Counter
import gzip
import pickle
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get LR stats",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-b",
        "--barcodes",
        type=str,
        required=True,
        help="scTagger detected barcodes TSV file.",
    )
    parser.add_argument(
        "-r",
        "--reads",
        type=str,
        nargs="+",
        required=True,
        help="FASTQ/FASTA file(s).",
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
        "--lr-matches",
        type=str,
        required=True,
        help="scTagger detected LR matches TSV file.",
    )
    parser.add_argument(
        "-l",
        "--lr-br",
        type=str,
        required=True,
        help="scTagger detected LR barcodes TSV file.",
    )
    parser.add_argument(
        "-o",
        "--outpath",
        type=str,
        required=True,
        help="Output path to pickle files which will be: <outpath>.<object>.pickle with <object> being barcodes, reads, seqs.",
    )
    args = parser.parse_args()
    return args


def get_lr_states(
    barcodes_path,
    reads_paths,
    paf_path,
    lr_matches_path,
    lr_br_path,
    pickle_path,
):
    # Get barcodes
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    barcode_to_bid = dict()
    barcodes = list()
    for bid, l in tqdm(
        enumerate(gzip.open(barcodes_path, "rt")), desc=f"Reading {barcodes_path}"
    ):
        barcode = dict(
            bid=bid,
            seq=l.strip().split("\t")[0],
            rids=list(),
        )
        rc = "".join(complement.get(base, base) for base in reversed(barcode["seq"]))

        assert barcode["seq"] not in barcode_to_bid and rc not in barcode_to_bid
        barcode_to_bid[barcode["seq"]] = bid
        barcode_to_bid[rc] = bid
        barcodes.append(barcode)

    # Get reads
    rname_to_rid = dict()
    reads = list()
    seqs = list()
    for fastq in reads_paths:
        if fastq.endswith(".gz"):
            opener = gzip.open(fastq, "rt")
        else:
            opener = open(fastq)
        for idx, line in enumerate(tqdm(opener, desc=f"Processing {fastq}")):
            if idx == 0:
                if line[0] == "@":  # fastq
                    mod = 4
                elif line[0] == ">":  # fasta
                    mod = 2
                else:
                    raise ValueError("Unknown fastq/a format")
            if idx % mod == 0:
                read = dict(
                    rid=len(rname_to_rid),
                    name=line[1:].split(" ")[0],
                    dist=-1,
                    bids=list(),
                    br_seg=(
                        -1,  # dist
                        0,  # loc
                        "",  # seg
                    ),
                    mappings=Counter(),
                    length=-1,
                )
                assert not read["name"] in rname_to_rid
                rname_to_rid[read["name"]] = read["rid"]
                reads.append(read)
            elif idx % mod == 1:
                read["length"] = len(line.strip())
                seqs.append(line.strip())
    # Add mappings
    for line in tqdm(open(paf_path), desc=f"Processing {paf_path}"):
        line = line.rstrip("\n").split("\t")
        if not "tp:A:P" in line:
            continue
        strand = line[4]
        tid = line[5]
        rid = rname_to_rid[line[0]]
        qlen = int(line[1])
        qstart = int(line[2])
        qend = int(line[3])
        tlen = int(line[6])
        tstart = int(line[7])
        tend = int(line[8])
        reads[rid]["mappings"][
            (tid, strand, qstart, qend, qlen, tstart, tend, tlen)
        ] += 1

    # Add barcode segment
    for line in tqdm(gzip.open(lr_br_path, "rt"), desc=f"Processing {lr_br_path}"):
        line = line.rstrip("\n").split("\t")
        rid = rname_to_rid[line[0]]
        if line[2] == "NA":
            continue
        reads[rid]["br_seg"] = (
            int(line[1]),
            int(line[2]),
            line[3],
        )

    # Add matches
    if lr_matches_path.endswith(".gz"):
        opener = gzip.open(lr_matches_path, "rt")
    else:
        opener = open(lr_matches_path)
    for line in tqdm(opener, desc=f"Processing {lr_matches_path}"):
        line = line.rstrip("\n").split("\t")
        rid = rname_to_rid[line[0]]
        reads[rid]["dist"] = int(line[1])
        reads[rid]["bids"] = [barcode_to_bid[b] for b in line[4].split(",")]
        for bid in reads[rid]["bids"]:
            barcodes[bid]["rids"].append(
                (
                    rid,
                    reads[rid]["dist"],
                    len(reads[rid]["bids"]),
                )
            )
    pickle.dump(
        barcodes,
        open(f"{pickle_path}.barcodes.pickle", "wb+"),
    )
    pickle.dump(
        reads,
        open(f"{pickle_path}.reads.pickle", "wb+"),
    )
    pickle.dump(
        seqs,
        open(f"{pickle_path}.seqs.pickle", "wb+"),
    )


def main():
    args = parse_args()
    get_lr_states(
        barcodes_path=args.barcodes,
        reads_paths=args.reads,
        paf_path=args.paf,
        lr_matches_path=args.lr_matches,
        lr_br_path=args.lr_br,
        pickle_path=args.outpath,
    )


if __name__ == "__main__":
    main()
