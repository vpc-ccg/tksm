
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV

import pickle


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
            if len(lens)>1000:
                break
    return lens, tlens




if __name__ == "__main__":
    m1, m2 = read_paf("/groups/hachgrp/projects/dev-simulation/code/RNAInFuser/experiments/sim-from-nanopore/MCF7-sgnex.paf")
    values = np.vstack([m1, m2])


    params = {"bandwidth": np.arange(50, 120, 10)}
    grid_finer = GridSearchCV(KernelDensity(), params, n_jobs=16,cv=3,verbose=3)
    grid_finer.fit(values.T)

    pickle.dump(grid_finer, open('grid.pickle','wb'))
    print(grid_finer.best_params_)
