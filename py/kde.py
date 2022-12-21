import pickle
import argparse
import sys
import time

from tqdm import tqdm
import numpy as np
from sklearn.neighbors import KernelDensity
from matplotlib import pyplot as plt

from multiprocessing import Pool

def run_batch_sk(prc_arg):
    X = prc_arg
    p = kde_vals.score_samples(X)
    return X, p

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def read_paf( paf):

    lens = []
    tlens = []
    with open(paf, 'r') as hand:
        for line in hand:
            if "tp:A:P" not in line:
                continue
            fields = line.rstrip().split("\t")
            start = int(fields[2])
            end = int(fields[3])
            lens.append(end-start)
            tlens.append(int(fields[6]))
    return lens, tlens


def sort_pp(X, Y, p):
    orderY = Y.argsort(kind='stable')
    XX = X[orderY]
    YY = Y[orderY]
    pp = p[orderY]
    orderX = XX.argsort(kind='stable')

    XX = XX[orderX]
    YY = YY[orderX]
    pp = pp[orderX]

    return XX,YY,pp
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="KDE computation module")
    parser.add_argument("-i", "--input", type=str, help="Path to paf file of the real dataset", required=True)
    parser.add_argument("-b", "--bandwidth", help="Bandwidth value")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to run KDE")


    parser.add_argument("--fro", type=int, default=100, help="KDE grid min")
    parser.add_argument("--to", type=int, default=100, help="KDE grid max")
    parser.add_argument("--spaces", type=int, default=5, help="KDE grid spaces between adjacent points")
    

    parser.add_argument("--grid", type=str, help="Path to grid.npy output")
    parser.add_argument("--labelsx", type=str, help="Path to labelsx.npy output")
    parser.add_argument("--labelsy", type=str, help="Path to labelsy.npy output")



    args = parser.parse_args(sys.argv[1:])
    
    print("Reading {}".format(args.input), file=sys.stderr)
    m1, m2 = read_paf(args.input)
    values = np.vstack([m1, m2])

    kd = KernelDensity(bandwidth=float(args.bandwidth))

    kde_vals = kd.fit(values.T)

    J = args.threads


    print("Computing KDE with {} threads".format(J), file=sys.stderr)
    X_base = np.arange( args.fro, args.to, args.spaces)
    Y_base = np.arange( args.fro, args.to, args.spaces)

    XY_r = []
    pp_r = []
    proc_args = []
    start = time.time()
    L = [(x,y) for x in X_base for y in Y_base]
    for x in chunks(L, len(L)//J):
        proc_args.append(x)
    with Pool(min(J,len(proc_args))) as p:
        for XY, p in tqdm(p.imap_unordered(run_batch_sk, proc_args), total=len(proc_args)):
            XY_r.extend(XY)
            pp_r.extend(p)
    XX = np.array(XY_r)[:,0]
    YY = np.array(XY_r)[:,1]
    pp = np.array(pp_r)
    end = time.time()

    print("Pickling the computed values", file=sys.stderr)
    pickle.dump((XX,YY,pp), open("temp.pickle","wb")) 
    print("Formatting the computed probabilities", file=sys.stderr)

    XXo,YYo,ppo = sort_pp(XX,YY,pp)

    np.save(args.grid, np.exp(ppo).reshape(X_base.shape[0],Y_base.shape[0]))

    np.save(args.labelsx, X_base)
    np.save(args.labelsy, Y_base)

    print(end - start)
