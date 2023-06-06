import pickle
from collections import defaultdict, Counter
import argparse

from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Compute gene fusion intersections",
    )
    parser.add_argument(
        "--tid_to_gid",
        type=str,
        required=True,
        help="Path to pickle file of tid_to_gid",
    )
    parser.add_argument(
        "--gname_to_gid",
        type=str,
        required=True,
        help="Path to pickle file of gname_to_gid",
    )
    parser.add_argument("--tsb", type=str, required=True, help="Path to Tsb MDF file")
    parser.add_argument("--mrg", type=str, required=True, help="Path to Mrg MDF file")
    parser.add_argument("--uns", type=str, required=True, help="Path to Uns MDF file")
    parser.add_argument(
        "--genion_pass",
        type=str,
        required=True,
        help="Path to Genion pass TSV file",
    )
    parser.add_argument(
        "--genion_fail",
        type=str,
        required=True,
        help="Path to Genion fail TSV file",
    )
    parser.add_argument(
        "--longGF",
        type=str,
        required=True,
        help="Path to LongGF TSV file",
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Path pickle output file"
    )
    args = parser.parse_args()
    return args


tid_to_gid = None
gname_to_gid = None


def get_mid_to_tid_from_mrg_mdf(mrg_mdf):
    mid_to_tid = dict()
    for l in tqdm(open(mrg_mdf), desc=f"Reading {mrg_mdf}"):
        if l[0] != "+":
            continue
        l = l[1:].rstrip("\n").split("\t")
        mid = l[0]
        tid = "NA"
        for t in l[2].split(";"):
            if t.startswith("tid="):
                tid = t[4:]
                break
        mid_to_tid[mid] = tid
    return mid_to_tid


def add_longGF_tsv(longGF_tsv, gene_fusions):
    for line in open(longGF_tsv, "r"):
        if "SumGF" in line:
            gnames = line.split()[1].split(":")
            if len(gnames) != 2:
                continue
            gids = [gname_to_gid[g] for g in gnames]
            gids = sorted(gids)
            gene_fusions[(gids[0], gids[1])].add("LongGF")


def add_genion_tsv(genion_tsv, gene_fusions):
    for l in tqdm(open(genion_tsv), desc=f"Reading {genion_tsv}"):
        l = l.rstrip("\n").split("\t")
        gids = [g for g in l[0].split("::")]
        if len(gids) != 2:
            continue
        if int(l[4]) < 3:
            continue
        for g1, g2 in zip(gids[:-1], gids[1:]):
            gene_fusions[(g1, g2)].add(l[6])


def add_tsb_mdf(tsb_mdf, gene_fusions):
    gf_counter = Counter()
    for l in tqdm(open(tsb_mdf), desc=f"Reading {tsb_mdf}"):
        if l[0] != "+":
            continue
        l = l[1:].rstrip("\n").split("\t")
        if "::" not in l[2]:
            continue
        mid = l[0]
        tid = "NA"
        for t in l[2].split(";"):
            if t.startswith("tid="):
                tids = t[4:].split(":")
                tids = [t for t in tids if t != ""]
                break
        for tid_1, tid_2 in zip(tids[:-1], tids[1:]):
            gid_1 = tid_to_gid[tid_1]
            gid_2 = tid_to_gid[tid_2]
            gid_1, gid_2 = sorted((gid_1, gid_2))
            gf_counter[(gid_1, gid_2)] += int(l[1])
    for k, v in gf_counter.items():
        if v < 3:
            continue
        gene_fusions[k].add("Truth")


def add_uns_mdf(uns_mdf, mid_to_tid, gene_fusions):
    gf_counter = Counter()
    for l in tqdm(open(uns_mdf), desc=f"Reading {uns_mdf}"):
        if l[0] != "+" or "Cat=" not in l:
            continue
        l = l[1:].rstrip("\n").split("\t")
        mid_1 = l[0]
        for t in l[2].split(";"):
            if t.startswith("Cat="):
                mids = [mid_1] + t[4:].split(",")
                for mid_1, mid_2 in zip(mids[:1], mids[1:]):
                    tid_1 = mid_to_tid[mid_1].split(":")[-1]
                    gid_1 = tid_to_gid[tid_1]
                    tid_2 = mid_to_tid[mid_2].split(":")[0]
                    gid_2 = tid_to_gid[tid_2]
                    if gid_1 == gid_2:
                        continue
                    gid_1, gid_2 = sorted((gid_1, gid_2))
                    gf_counter[(gid_1, gid_2)] += int(l[1])
    for k, v in gf_counter.items():
        if v < 3:
            continue
        gene_fusions[k].add("Glue")


def main():
    args = parse_args()
    global tid_to_gid, gname_to_gid
    tid_to_gid = pickle.load(open(args.tid_to_gid, "rb"))
    gname_to_gid = pickle.load(open(args.gname_to_gid, "rb"))
    gene_fusions = defaultdict(set)

    mid_to_tid = get_mid_to_tid_from_mrg_mdf(args.mrg)

    add_tsb_mdf(args.tsb, gene_fusions)
    add_uns_mdf(args.uns, mid_to_tid, gene_fusions)
    add_genion_tsv(args.genion_pass, gene_fusions)
    add_genion_tsv(args.genion_fail, gene_fusions)
    pickle.dump(gene_fusions, open(args.output, "wb"))
    return 0


if __name__ == "__main__":
    exit(main())
