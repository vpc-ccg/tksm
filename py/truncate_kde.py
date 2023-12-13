import argparse
import functools
from multiprocessing import Pool

from tqdm import tqdm
import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import sys


kde_vals = None


def parse_args():
    parser = argparse.ArgumentParser(
        description="KDE computation module of truncation of transcriptomic long-reads using their cDNA mapping.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="PAF file of the mapping of the long-reads to the reference transcriptome.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path prefix for the output files: <>.X_idxs.npy, <>.Y_idxs.npy, <>.grid.npy, <>.sider.tsv",
    )
    parser.add_argument(
        "-b",
        "--bandwidth",
        type=float,
        default=100.0,
        help="Bandwidth value for the KDE. If set to -1, script will use cross validation to compute it.",
    )
    parser.add_argument(
        "--grid-start",
        type=int,
        default=0,
        help="Read/transcript length start of the KDE grid.",
    )
    parser.add_argument(
        "--grid-end",
        type=int,
        default=10000,
        help="Read/transcript length end of the KDE grid.",
    )
    parser.add_argument(
        "--grid-step",
        type=int,
        default=100,
        help="Read/transcript length step of the KDE grid.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads to run KDE and GridSearchCV.",
    )
    parser.add_argument(
        "--model-lengths",
        default=False,
        action='store_true',
        help="Model read lengths instead of truncation lengths",
    )
    parser.add_argument(
        "--model-separate",
        default=False,
        action='store_true',
        help="Model read lengths instead of truncation lengths",
    )
    parser.add_argument(
        "--model-3D",
        default=False,
        action='store_true',
        help="Model read lengths instead of truncation lengths",
    )
    class ListPrinter(argparse.Action):
        def __call__(self, parser, namespace, values, option_string):
            txt = "\n".join([getattr(k, "dest") for k in parser._actions])
            print(txt)
            parser.exit()

    parser.add_argument("--list", nargs=0, action=ListPrinter)
    parser.add_argument(
        "--end-ratio",
        type=float,
        default=-1,
        help="Ratio of truncation from the end of the molecule. "
        + "If set to -1, TKSM will build a distribution of the end ratios from the PAF file."
        + "If set between 0 and 1, TKSM will use this value as the end ratio for all molecules.",
    )
    args = parser.parse_args()
    if args.end_ratio != -1:
        assert 0 <= args.end_ratio <= 1
    return args


def score_samples_runner(xy):
    return xy, kde_vals.score_samples(xy)



def get_truncation_lens_paired_with_transcript_lens_53(paf):
    truncation_lengths3 = []
    truncation_lengths5 = []
    transcript_lens = []
    for line in tqdm(open(paf, "r")):
        if "tp:A:P" not in line:
            continue
        line = line.rstrip("\n").split("\t")
        strand = line[4]
        tlen = int(line[6])
        start = int(line[7])
        end = int(line[8])

        truncation_lengths5.append(start)
        truncation_lengths3.append(tlen-end)
        transcript_lens.append(tlen)
    return truncation_lengths5, truncation_lengths3, transcript_lens

def get_truncation_lens_paired_with_transcript_lens_X4(paf):
    truncation_lengths3 = [[],[]]
    truncation_lengths5 = [[],[]]
    transcript_lens = [[],[]]
    for line in tqdm(open(paf, "r")):
        if "tp:A:P" not in line:
            continue
        line = line.rstrip("\n").split("\t")
        strand = line[4]
        tlen = int(line[6])
        start = int(line[7])
        end = int(line[8])
        idx = 0 if strand=="+" else 1
        truncation_lengths5[idx].append(start)
        truncation_lengths3[idx].append(tlen-end)
        transcript_lens[idx].append(tlen)
    return truncation_lengths5, truncation_lengths3, transcript_lens

def get_truncation_lens_paired_with_transcript_lens(paf):
    truncation_lengths = list()
    transcript_lens = list()
    end_rati = list()
    for line in tqdm(open(paf, "r")):
        if "tp:A:P" not in line:
            continue
        line = line.rstrip("\n").split("\t")
        strand = line[4]
        tlen = int(line[6])
        start = int(line[7])
        end = int(line[8])
        qstart = int(line[2])
        qend = int(line[3])
        qlen = int(line[1])

        #truncation_length = qstart + qlen - qend # head_truncation (start) + tail_truncation  (tlen-end)
        truncation_length = start + tlen - end # head_truncation (start) + tail_truncation  (tlen-end)
        transcript_lens.append(tlen)
        truncation_lengths.append(truncation_length)
        if truncation_length != 0:
            end_truncation = qlen - qend if strand == "+" else qstart

            end_rati.append(end_truncation/truncation_length)

    return truncation_lengths, transcript_lens, end_rati

def get_alignment_lens(paf):
    tlens = list()
    alens = list()
    end_ratios = list()
    for line in tqdm(open(paf, "r")):
        if "tp:A:P" not in line:
            continue
        line = line.rstrip("\n").split("\t")
        strand = line[4]
        tlen = int(line[6])
        start = int(line[7])
        end = int(line[8])
        alen = end - start
        tlens.append(tlen)
        alens.append(alen)
        trunc_len = tlen - alen
        if trunc_len == 0:
            continue
        if strand == "+":
            end_trunc = tlen - end
            end_ratios.append(end_trunc / trunc_len)
        else:
            end_trunc = start
            end_ratios.append(end_trunc / trunc_len)
    return tlens, alens, end_ratios


def sort_pp(X, Y, p):
    orderY = Y.argsort(kind="stable")
    XX = X[orderY]
    YY = Y[orderY]
    pp = p[orderY]

    orderX = XX.argsort(kind="stable")
    XX = XX[orderX]
    YY = YY[orderX]
    pp = pp[orderX]

    return XX, YY, pp

def CV_KDE_bandwidth(len_values,args ):
    bandwidths = list()
    print(
        "Non-positive bandwidth selected selected: recomputing bandwidth with GridSearchCV"
    )
    for _ in range(3):
        grid_finder = GridSearchCV(
            KernelDensity(),
            {"bandwidth": np.arange(50, 1000, 100)},
            n_jobs=args.threads,
            cv=3,
            verbose=3,
        )
        grid_finder.fit(
            len_values[np.random.randint(len_values.shape[0], size=100_000), :]
        )
        bandwidths.append(grid_finder.best_params_["bandwidth"])
        print(bandwidths[-1])
    print(bandwidths)
    return np.median(bandwidths)

def ComputeKDELikelihoods(len_values, args):
    kd = KernelDensity(bandwidth=args.bandwidth)
    global kde_vals
    kde_vals = kd.fit(len_values)

    print("Computing KDE with {} threads".format(args.threads))

    X_idxs = np.arange(args.grid_start, args.grid_end + 1, args.grid_step)
    Y_idxs = np.arange(args.grid_start, args.grid_end + 1, args.grid_step)

    XY = list()
    P = list()
    linearized_grid = list()
    for i, x in enumerate(X_idxs[:-1]):
        for j, y in enumerate(Y_idxs[:-1]):
            linearized_grid.append(
                [
                    (
                        (x + X_idxs[i + 1]) // 2,
                        (y + Y_idxs[j + 1]) // 2,
                    )
                ]
            )
    if args.threads > 1:
        p = Pool(args.threads)
        mapper = functools.partial(p.imap_unordered, chunksize=10)
    else:
        mapper = map
    for xy, log_likelihood in tqdm(
        mapper(score_samples_runner, linearized_grid), total=len(linearized_grid)
    ):
        XY.append(xy)
        P.append(log_likelihood)
    if args.threads > 1:
        p.close()
    XY = np.array(XY).reshape(len(linearized_grid), 2)
    X = np.array(XY)[:, 0]
    Y = np.array(XY)[:, 1]
    P = np.array(P)
    X, Y, P = sort_pp(X, Y, P)
    P = np.exp(P).reshape(X_idxs.shape[0] - 1, Y_idxs.shape[0] - 1)

    return X_idxs, Y_idxs, P

def PrintEndRatios(end_ratios, args):
    if args.end_ratio != -1:
        end_ratios = [args.end_ratio] * len(end_ratios)
    with open(f"{args.output}.sider.tsv", "w+") as outfile:
        for w, v in zip(*np.histogram(end_ratios, bins=np.arange(0, 1.01, 0.01))):
            outfile.write(f"{w:d}\t{v:.5f}\n")


def main():
    args = parse_args()

    print("Reading {}".format(args.input))
    
    if args.model_lengths:
        print("Modelling read lengths", file=sys.stderr)
        tlens, alens, end_ratios = get_alignment_lens(args.input)
        len_values = np.vstack([tlens, alens]).T
        args.bandwidth = args.bandwidth if args.bandwidth > 0 else CV_KDE_bandwidth(len_values, args)
        PrintEndRatios(end_ratios, args)
        X_idxs, Y_idxs, P = ComputeKDELikelihoods(len_values, args)

        print("Writing output...")
        np.save(f"{args.output}.X_idxs.npy", X_idxs)
        np.save(f"{args.output}.Y_idxs.npy", Y_idxs)
        np.save(f"{args.output}.grid.npy", P)

    elif args.model_separate:
        print("Modelling every truncation type")
        trc5, trc3, tlens = get_truncation_lens_paired_with_transcript_lens_X4(args.input) 

        args.bandwidth = args.bandwidth if args.bandwidth > 0 else CV_KDE_bandwidth(len_values, args)
        
        for strand, trunc_lens, transcript_lens in zip(["+", "-"], trc5, tlens):
            len_values = np.vstack([transcript_lens, trunc_lens]).T
            print(f"Writing output {strand}-5'...")
            X_idxs, Y_idxs, P = ComputeKDELikelihoods(len_values, args)
            np.save(f"{args.output}.grid{strand}5.npy", P)
        for strand, trunc_lens, transcript_lens in zip(["+", "-"], trc3, tlens):
            len_values = np.vstack([transcript_lens, trunc_lens]).T
            print(f"Writing output {strand}-3'...")
            X_idxs, Y_idxs, P = ComputeKDELikelihoods(len_values, args)
            np.save(f"{args.output}.grid{strand}3.npy", P)
        
        np.save(f"{args.output}.X_idxs.npy", X_idxs)
        np.save(f"{args.output}.Y_idxs.npy", Y_idxs)
    elif args.model_3D:
        print("Modelling 3D truncation type")
        trc5, trc3, tlens = get_truncation_lens_paired_with_transcript_lens_53(args.input) 
        bins, edges = np.histogramdd((trc5, trc3, tlens), bins= args.grid_step, density=True,
                              range=((args.grid_start, args.grid_end),
                                     (args.grid_start, args.grid_end),
                                     (args.grid_start, args.grid_end)) )
        np.save(f"{args.output}.grid.npy", bins)
        np.save(f"{args.output}.edges.npy", edges)
    else:
        print("Modelling truncation lengths",file=sys.stderr)
        tlens, alens, end_ratios = get_truncation_lens_paired_with_transcript_lens(args.input)
        len_values = np.vstack([tlens, alens]).T
        args.bandwidth = args.bandwidth if args.bandwidth > 0 else CV_KDE_bandwidth(len_values, args)
        PrintEndRatios(end_ratios, args)
        X_idxs, Y_idxs, P = ComputeKDELikelihoods(len_values, args)
        print("Writing output...")
        np.save(f"{args.output}.X_idxs.npy", X_idxs)
        np.save(f"{args.output}.Y_idxs.npy", Y_idxs)
        np.save(f"{args.output}.grid.npy", np.transpose(P))

"""
    print(f"Using bandwidth = {args.bandwidth}")
    
    X_idxs, Y_idxs, P = ComputeKDELikelihoods(len_values, args)

    print("Writing output...")
    np.save(f"{args.output}.X_idxs.npy", X_idxs)
    np.save(f"{args.output}.Y_idxs.npy", Y_idxs)
    if args.model_lengths:
        np.save(f"{args.output}.grid.npy", P)
    else:
        np.save(f"{args.output}.grid.npy", np.transpose(P))
"""

if __name__ == "__main__":
    main()
