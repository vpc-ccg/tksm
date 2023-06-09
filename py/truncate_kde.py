import argparse
import functools
from multiprocessing import Pool

from tqdm import tqdm
import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV

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
        help="Path prefix for the output files: <>.X_idxs.npy, <>.Y_idxs.npy, <>.grid.npy",
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
        default=100,
        help="Read/transcript length start of the KDE grid.",
    )
    parser.add_argument(
        "--grid-end",
        type=int,
        default=5000,
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
    args = parser.parse_args()
    return args


def score_samples_runner(xy):
    return xy, kde_vals.score_samples(xy)


def get_alignment_lens(paf):
    tlens = list()
    alens = list()
    for line in tqdm(open(paf, "r")):
        if "tp:A:P" not in line:
            continue
        line = line.rstrip("\n").split("\t")
        tlen = int(line[6])
        alen = int(line[3]) - int(line[2])
        tlens.append(tlen)
        alens.append(alen)
    return tlens, alens


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


def main():
    args = parse_args()

    print("Reading {}".format(args.input))

    tlens, alens = get_alignment_lens(args.input)
    len_values = np.vstack([tlens, alens]).T
    if args.bandwidth <= 0:
        print(
            "Non-positive bandwidth selected selected: recomputing bandwidth with GridSearchCV"
        )
        bandwidths = list()
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
        args.bandwidth = np.median(bandwidths)
    print(f"Using bandwidth = {args.bandwidth}")
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

    print("Writing output...")
    np.save(f"{args.output}.X_idxs.npy", X_idxs)
    np.save(f"{args.output}.Y_idxs.npy", Y_idxs)
    np.save(f"{args.output}.grid.npy", P)


if __name__ == "__main__":
    main()
