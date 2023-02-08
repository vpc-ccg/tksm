
import sys
import gzip
import pickle
import pathlib 
import os
import random
import time
import json

from joblib import Parallel, delayed
from multiprocessing import Pool
from tqdm import tqdm
from sklearn.neighbors import KernelDensity
from sklearn.preprocessing import normalize
import numpy as np
from matplotlib import pyplot as plt


swirl = [(-1)**i * i for i in range(1,101)]


def run_batch_sk(prc_arg):
    kde_vals = prc_arg[0]
    X = prc_arg[1]
    p = kde_vals.score_samples(X)

    return X, normalize(p, axis=0, norm="l1")

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def sort_pp(X, Y, p):
    orderY = Y.argsort(kind='stable')
    XX = X[orderY]
    YY = Y[orderY]
    pp = p[orderY]
    orderX = XX.argsort(kind='stable')

    XX = XX[orderX]
    YY = YY[orderX]
    pp = pp[orderX]

    return XX, YY, pp

class CustomDist:
    def __init__(self, pdf, indices, ):
        self.pdf = pdf
        self.indices = indices
        
        sum_pdf = np.sum(self.pdf)
        self.cdf = np.zeros(len(self.pdf))
        for i, v in enumerate(self.pdf):
            self.cdf[i] = v / sum_pdf + self.cdf[i-1]
        
    def __call__(self):
        val = random.uniform(0,1)
        pos = np.searchsorted(self.cdf, val)
        return self.indices[pos]

class Custom2Dist:
    @classmethod
    def from_file(cls, path):
        lx = np.load(f'{path}/labelsx.npy')
        ly = np.load(f'{path}/labelsy.npy')
        grid = np.load(f'{path}/grid.npy').T
        return cls(lx,ly,grid)
    @classmethod
    def from_dict(cls, dc):
        lx = np.array(dc["lx"])
        ly = np.array(dc["ly"])
        grid = np.array(dc["grid"])
        return cls(lx,ly,grid)
    def __init__(self, lx, ly, grid):
        self.lx = lx
        self.ly = ly
        self.grid = grid
        self.distributions = []
        for i, y in enumerate(self.ly):
            disto = CustomDist(self.grid[i,:], self.lx)
            self.distributions.append(disto)
    def serialize(self):
        return {"lx" : self.lx.tolist(), 
                "ly" : self.ly.tolist(),
                "grid" : self.grid.tolist()}
    def save(self, file):
        json.dump( self.serialize(), file)
    @classmethod
    def load(cls, file):
        return cls.from_dict(json.load(file))
    def __call__(self, y):
        pos = np.searchsorted(self.ly,y)
        if pos < len(self.ly) - 1:
            if np.abs(self.ly[pos] - y) > np.abs(self.ly[pos+1] - y):
                pos+=1
        if pos >= len(self.ly):
            mult = pos / self.ly[-1]
            pos = len(self.ly) - 1
        else:
            mult = 1
        return int(self.distributions[pos]() * mult)
    def __getitem__(self, y): #returns the subprobability dist without calling
        pos = np.searchsorted(self.ly,y)
        if pos < len(self.ly) - 1:
            if np.abs(self.ly[pos] - y) > np.abs(self.ly[pos+1] - y):
                pos+=1
        if pos >= len(self.ly):
            mult = y / self.ly[-1]
            pos = len(self.ly) - 1
        else:
            mult = 1
        return lambda : int(mult * self.distributions[pos]())

class Mock_noise_generator:
    def __init__(self):
        pass
    def sample(self, x):
        return 0
    def noise_seq(self, i):
        return ""

class KDE_noise_generator:
    @classmethod
    def from_data(cls, mapped, unmapped, lx, ly, transition_matrix, bw, threads=64):
        kde_vals = KernelDensity(bandwidth=bw).fit(list(zip(mapped,unmapped)))
        
        full_grid = np.array( [[x, y] for x in lx for y in ly]).reshape(lx.shape[0],ly.shape[0],2)
        results=Parallel(n_jobs=threads)(delayed(kde_vals.score_samples)(row) for row in full_grid)

        grid = np.exp(np.stack(results))
        ratio = np.sum(np.array(unmapped) > 0) / len(unmapped)

        return cls(lx, ly, grid, transition_matrix, ratio)

    def __init__(self, lx, ly, grid, transition_matrix, ratio, bases=list("AGTC")):
        self.ly = np.array(ly)
        self.lx = np.array(lx)
        self.grid = np.array(grid)
        self.ratio = ratio
        self.begin_W = np.array(transition_matrix[0])
        self.transition_matrix = np.array(transition_matrix[1])
        
        self.disto = Custom2Dist(self.lx, self.ly, self.grid)
        self.BASES = bases
    

    def sample(self, x):
        return self.disto(x)
    

    def noise_seq(self, i):
        if random.uniform(0,1) > self.ratio:
            return ""
        x = self.sample(i)
        if( x == 0):
            return ""
        if (isinstance(x, list) and len(x) == 1):
            x = x[0]
        index = random.choices(np.arange(0,4,dtype=np.int32))[0]
        gen_seq = ["A"]*x
        for i in range(x):
            w = self.transition_matrix[index,:]
            index = random.choices(np.arange(0,4,dtype=np.int32),w)[0]
            gen_seq[i] = self.BASES[index]
        return "".join(gen_seq)
    
    def save(self,file):
        dc = self.disto.serialize()
        dc["bases"] = self.BASES
        dc["ratio"] = self.ratio
        dc["trans"] = self.transition_matrix.tolist()
        dc["begin"] = self.begin_W.tolist()
        json.dump(dc, file)
    @classmethod
    def load(cls, model_type_or_filename):
        this_script_dir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))

        if model_type_or_filename == 'nanopore':

            dc = json.load(gzip.open(str(this_script_dir / 'tail_models' / 'nanopore.json.gz'),'rt'))
        elif model_type_or_filename.lower() in ["none", "pacbio"]:
            return Mock_noise_generator()       
        else:

            dc = json.load(gzip.open(model_type_or_filename,'rt'))
        return cls(dc["lx"], dc["ly"], dc["grid"], (dc["begin"], dc["trans"]), dc["ratio"], dc["bases"])

class KDE_noise_sampler:
    
    def get_noise_lengths(self, N):
        _samp = self.kde_vals.sample(N)
        _samp = sorted(_samp, key=lambda x: x[0])
        X, Y = zip(*_samp)
        return X, Y
    
    def get_bin(self,x):
        return (self.bins[self.ind[x]-1],self.X[x])

    @classmethod
    def from_data(cls, mapped, unmapped, b, W, bin_len=10, bw=50, max_to_skip = 100):
        kde_vals = KernelDensity(bandwidth=bw).fit(list(zip(mapped,unmapped)))
        ratio = np.sum(np.array(unmapped) > 0) / len(unmapped)
        return cls(kde_vals, b, W, bin_len, ratio)
    
    
    def __init__(self, kde_vals, b, W, bin_len=10, none=False, ratio = 0.1):
        self.kde_vals = kde_vals
        self.values = {}
        self.bin_len=bin_len
        self.W = W
        self.begin_W =  b.reshape(-1)
        self.resample()
        self.none =none
        self.ratio = ratio
        self.BASES = list("AGTC")
    @classmethod
    def load_from_file(cls, path):
        opener = open
        open_mode = "r"
        if path.endswith(".gz"):
            opener = gzip.open
            open_mode = "rb"
        

        beg, tm, kde_model = pickle.load(opener(path, open_mode))
        return beg, tm, kde_model


    @classmethod
    def frompickle(cls, model_type_or_filename, bin_len=10):
        this_script_dir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))

        if model_type_or_filename == 'nanopore':
            none = False
            beg, transition_matrix, kde_model, ratio = cls.load_from_file(str(this_script_dir / 'tail_models' / 'nanopore.pickle.gz'))
        elif model_type_or_filename == 'pacbio': #Pacbio doesnt have I think
#            beg, transition_matrix, kde_model = cls.load_from_file(str(this_script_dir / 'tail_models' / 'pacbio.pickle.gz'))
            none = True
            transition_matrix =[]
            kde_model = []
            beg = []
            ratio = 0
        elif model_type_or_filename == 'none':
            none = True
            transition_matrix =[]
            kde_model = []
            beg = []
            ratio = 0
        else:
            none = False
            beg, transition_matrix, kde_model, ratio = cls.load_from_file(model_type_or_filename, output)
        return cls(kde_model, beg, transition_matrix, bin_len, none, ratio)
    
    def resample(self, N=10000):
        self.X, self.Y = self.get_noise_lengths(N)
        self.bins = np.arange(0,np.max(self.X),self.bin_len, dtype=np.int32)
        self.ind = np.digitize(self.X,bins=self.bins)
        self.values.clear()
        for i,y in enumerate(self.Y):
            b, x = self.get_bin(i)
            if b not in self.values:
                self.values[b] = []
            self.values[b].append(int(np.round(y)))

    def sample(self, x):
        if self.none:
            return 0
        ind = np.digitize(x,bins=self.bins) - 1
        ind = self.bins[ind]
        if not isinstance(ind,list):
            ind=[ind]
        ret = []
        for i in ind:
            for j in range(50):#try resampling 50 times
                if i not in self.values:
                    self.resample()
                else:
                    break
            if j == 49:
                enum = 0
                while i not in self.values:
                    i+= swirl[enum] * (self.bin_len)
                    enum+=1
                    while i in self.values and len(self.values[i]) == 0:
                        i+= swirl[enum] * (self.bin_len)
                        enum+=1
            l = self.values[i]
            while len(l) == 0:
                self.resample()
                if i not in self.values:
                    continue
                l = self.values[i]
            ret.append(l.pop())
        return ret
    def noise_seq(self, i):
        if random.uniform(0,1) > self.ratio:
            return ""
        x = self.sample(i)
        if( x == 0):
            return ""
        if (isinstance(x,list) and len(x) == 1):
            x = x[0]
        index = random.choices(np.arange(0,4,dtype=np.int32),self.begin_W)[0]
        gen_seq = ["A"]*x
        for i in range(x):
            w = self.W[index,:]
            index = random.choices(np.arange(0,4,dtype=np.int32),w)[0]
            gen_seq[i] = self.BASES[index]
        return "".join(gen_seq)
#Loads primary alignments from paf file skips reads with multiple alignments
def load_paf(paf_path):
    reads = {}
    banned = set()
    with open(paf_path, 'r') as hand:
        for line in hand:
            if "tp:A:P" not in line:
                continue
            fields = line.split("\t")
            rid = fields[0]
            if rid in banned:
                continue
            if rid in reads:
                banned.add(rid)
                del reads[rid]
                continue
            if int(fields[1]) - int(fields[3]) > int(fields[3]):
                continue
            reads[rid] = (int(fields[3]),int(fields[1]),int(fields[7])) #q-end, q-length, t-start
    return reads, banned

def load_seq(fastq_path, reads, banned, first_N=1000000):
    opener = open
    open_mode = "r"
    if fastq_path.endswith(".gz"):
        opener = gzip.open
        open_mode = "rb"

    seqs = {}
    c = 0

    with opener(fastq_path, open_mode) as hand:
        line =  hand.readline()
        while line:

            fields = line.split()
            rid = fields[0][1:]
            if (isinstance(rid,bytes)):
                rid = rid.decode("utf-8")
            seq =  hand.readline()
            if (isinstance(seq,bytes)):
                seq = seq.decode("utf-8")
            if rid not in banned and rid in reads:
                noise_interval = reads[rid]
                seqs[rid] = seq[noise_interval[0]:]
                c+=1
                if c > first_N:
                    break
            line =  hand.readline()
            line =  hand.readline()
            line =  hand.readline()
    return seqs
            

def make_tail_noise_model(args, output=sys.stderr):


    reads, banned = load_paf(args.alignment)
    if args.reads:
        seqs = load_seq(args.reads, reads, banned)
    else:
        seqs = {}
    
    mapped_vs_unmapped = [(v[0],v[1] - v[0]) for k,v in reads.items()]
    mapped, unmapped = list(zip(*mapped_vs_unmapped[:args.sample_size]))
    BASES = list("AGTC")
    base_index = { a : i for i, a in enumerate(BASES)}
    transition_matrix = np.ones((len(BASES),len(BASES)))
    beg = np.zeros(len(BASES))
    for rid, seq in seqs.items():
        if(seq == "\n"):
            continue
        beg[base_index[seq[0]]] +=1
        for i,c in enumerate(seq.rstrip()):
            if i > 0:
                transition_matrix[base_index[seq[i-1]],base_index[c]] +=1

    transition_matrix = transition_matrix / np.sum(transition_matrix, axis=1)
    fro, to, spaces = args.range
    #Convert KDE to a 2D probability map
    X_base = np.arange( fro, to, spaces)
    Y_base = np.arange( fro, to, spaces)

    kdm = KDE_noise_generator.from_data(mapped, unmapped, X_base, Y_base, (beg,transition_matrix), args.bw, args.threads) 

    kdm.save(gzip.open(args.out, 'wt', encoding='ascii'))


if __name__ == "__main__":
    #Testing the tail noise sampler
    sampler = KDE_noise_sampler.frompickle("nanopore")
    
    mock_rl_dist_sample = np.random.lognormal(np.log(1200), 0.45, 1000)

    noise_lens = np.array([ 0 if random.uniform(0.0,1.0) < 0.9 else len(sampler.noise_seq(x)) for x in tqdm(mock_rl_dist_sample)])
#    noise_lens = np.zeros(len(mock_rl_dist_sample))
    fig,axs = plt.subplots(4, figsize=(8,20))
    axs[0].hist(mock_rl_dist_sample,bins=np.arange(0,8000,50))
    axs[1].hist(noise_lens + mock_rl_dist_sample,bins=np.arange(0,8000,50))
    axs[2].hist(noise_lens,bins=np.arange(0,8000,50))
    axs[3].scatter(mock_rl_dist_sample, noise_lens)

    axs[0].set_title("Random sample")
    axs[1].set_title("Generated Noise + sample")
    axs[2].set_title("Generated Noise")


    plt.savefig("dist.test.png")
