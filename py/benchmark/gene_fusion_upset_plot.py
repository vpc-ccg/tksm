import pickle
from collections import defaultdict, Counter
import argparse

from tqdm import tqdm
import upsetplot
from matplotlib import pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Upset plot for gene fusion",
    )
    parser.add_argument(
        "--pickle",
        type=str,
        help="Pickle file",
        required=True,
    )
    parser.add_argument(
        "--sample",
        type=str,
        help="Sample name",
        required=True,
    )
    parser.add_argument(
        "--outpath",
        type=str,
        help="Output path",
        required=True,
    )
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    S = args.sample
    try:
        if S.endswith("_trc"):
            trc = True
            S = S.rstrip("_trc")
        else:
            trc = False
        S = S.lstrip("TKSM_gene_fusion")
        rate = float(S) / 100
        rate = f"{rate:.0%}"
    except Exception as e:
        print(e)
        rate = "NA"
        trc = False

    gene_fusions = pickle.load(open(args.pickle, "rb"))
    gf_counter = Counter()
    for k, v in gene_fusions.items():
        KEY = set()
        for x in v:
            if x.startswith("PASS"):
                KEY.add("Genion pass")
            elif x.startswith("FAIL"):
                KEY.add("Genion fail")
            else:
                KEY.add(x)
        KEY = tuple(sorted(KEY))
        gf_counter[KEY] += 1
    UPSET_DATA = gf_counter.items()
    UPSET_DATA = upsetplot.from_memberships(
        list(zip(*gf_counter.items()))[0],
        list(zip(*gf_counter.items()))[1],
    )
    fig = plt.figure(figsize=(10, 10))
    upsetplot.plot(
        UPSET_DATA,
        fig=fig,
        show_counts=True,
    )
    fig.suptitle(f"Rate = {rate}; Trc = {trc}")
    fig.savefig(f"{args.outpath}.pdf", dpi=500)
    fig.savefig(f"{args.outpath}.png", dpi=1000)


if __name__ == "__main__":
    main()
