import sys
import gzip
from multiprocessing import Pool

import re
from itertools import groupby
import functools

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from snakemake.utils import min_version
from tqdm import tqdm

min_version("6.0")

if len(config)==0:
    configfile: 'config.yaml'
config['outpath'] = f'{config["outpath"]}/benchmark'

module RI_smk:
    snakefile: "Snakefile"
    config: config


outpath   = config['outpath']
preproc_d = f'{outpath}/preprocess'
RI_d      = f'{outpath}/RI'
nanosim_d = f'{outpath}/NS'
plots_d   = f'{outpath}/plots'

    
use rule * from RI_smk as RI_*

def experiment_prefix(exprmnt):
    prefix = list()
    for module in config["experiments"][exprmnt]["pipeline"]:
        prefix.append(list(module)[0])
    return '.'.join(prefix)

def longest_polys(seq, s, e, step, match_score=1, mismatch_score=-2, char='A'):
    if e-s == 0:
        return
    if seq[s] == char:
        scores = [match_score]
    else:
        scores = [0]
    for m in (match_score if c == char else mismatch_score for c in seq[s+step:e:step]):
        scores.append(max(0, scores[-1]+m))
    for k, g in groupby(enumerate(scores), lambda x: x[1] > 0):
        if not k:
            continue
        i, S = list(zip(*g))
        max_s, max_i = max(zip(S, i))
        l = max_i+1-i[0]
        yield i[0], l, seq[s:e:step][i[0]:i[0]+l].count(char)/l

def run_longest_polys(seq):
    L = 0
    for c in ['A','T']:
        for i,l,p in longest_polys(seq, 0, len(seq), 1):
            if p > .8:
                L = max(L,l)
    return L

rule all:
    input:
        # [
        #     f'{RI_d}/{exprmnt}/{experiment_prefix(exprmnt)}.fastq'
        #     for exprmnt in config['experiments']
        # ],
        [
            f'{plots_d}/{sample}/substition_analysis.tsv'
            for sample in config['samples']
        ],
        f'{plots_d}/all_raw_lengths.png',
        f'{plots_d}/all_mapping_raw_lengths.png',
        f'{plots_d}/all_mapped_lengths.png',
        f'{plots_d}/all_error_dist.png',
        f'{plots_d}/all_polyA.png',
    default_target: True

rule plot_all_raw_lengths:
    input:
        fastqs = [config['samples'][s]['fastq'] for s in config['samples']],
    output:
        png = f'{plots_d}/all_raw_lengths.png',
    params:
        bins=50
    run:
        fig = plt.figure()
        bb = None
        lens_arrays = list()
        max_up = 0
        samples = list(config['samples'].keys())
        print('Reading FASTQ/A files')
        for fastq in tqdm(input.fastqs):
            lens = list()
            if fastq.endswith('.gz'):
                infile = gzip.open(fastq, 'rt')
            else:
                infile = open(fastq, 'rt')
            for idx,line in tqdm(enumerate(infile)):
                if idx == 0:
                    mod = 4 if line[0]=='@' else 2
                if idx % mod == 1:
                    lens.append(len(line)-1)
            lens = np.array(lens)
            lens_arrays.append(lens)
            max_up = max(max_up, np.percentile(lens, 99))
        for lens,sample in zip(lens_arrays, samples):
            lens = lens[lens < max_up]
            _, bb, _ = plt.hist(
                lens,
                bins = params.bins if bb is None else bb,
                density = True,
                label = sample,
                alpha = .3,
            )

        plt.legend()
        plt.title("Length distribution (whole reads)")
        plt.savefig(output[0],dpi=300)        

rule plot_all_mapping_raw_lengths:
    input:
        pafs = [f'{preproc_d}/minimap2/{s}.cDNA.paf' for s in config['samples']],
    output:
        png = f'{plots_d}/all_mapping_raw_lengths.png',
    params:
        bins=50
    run:
        fig = plt.figure()
        bb = None
        lens_arrays = list()
        max_up = 0
        rids = set()
        samples = list(config['samples'].keys())
        print('Reading PAF files')
        for paf in tqdm(input.pafs):
            lens = list()
            for line in tqdm(open(paf)):
                    if "tp:A:P" not in line:
                        continue
                    fields = line.rstrip('\n').split('\t')
                    rid = fields[0]                    
                    if rid in rids:
                        continue
                    rids.add(rid)
                    qlen = int(fields[1])
                    lens.append(qlen)
            lens = np.array(lens)
            lens_arrays.append(lens)
            max_up = max(max_up, np.percentile(lens, 99))
        for lens,sample in zip(lens_arrays, samples):
            lens = lens[lens < max_up]
            _, bb, _ = plt.hist(
                lens,
                bins = params.bins if bb is None else bb,
                density = True,
                label = sample,
                alpha = .3,
            )
        plt.xlim(left=0, right=max_up)
        plt.legend()
        plt.title("Length distribution (whole reads from PAF)")
        plt.savefig(output[0],dpi=300)        

rule plot_all_mapped_lengths:
    input:
        pafs = [f'{preproc_d}/minimap2/{s}.cDNA.paf' for s in config['samples']],
    output:
        png = f'{plots_d}/all_mapped_lengths.png',
    params:
        bins=50
    run:
        fig = plt.figure()
        bb = None
        lens_arrays = list()
        max_up = 0
        rids = set()
        samples = list(config['samples'].keys())
        print('Reading PAF files')
        for paf in tqdm(input.pafs):
            lens = list()
            for line in tqdm(open(paf)):
                    if "tp:A:P" not in line:
                        continue
                    fields = line.rstrip('\n').split('\t')
                    rid = fields[0]                    
                    if rid in rids:
                        continue
                    rids.add(rid)
                    start = int(fields[2])
                    end = int(fields[3])
                    lens.append(end-start)
            lens = np.array(lens)
            lens_arrays.append(lens)
            max_up = max(max_up, np.percentile(lens, 99))
        for lens,sample in zip(lens_arrays, samples):
            lens = lens[lens < max_up]
            _, bb, _ = plt.hist(
                lens,
                bins = params.bins if bb is None else bb,
                density = True,
                label = sample,
                alpha = .3,
            )
        plt.xlim(left=0, right=max_up)
        plt.legend()
        plt.title("Length distribution (mapping part of reads from PAF)")
        plt.savefig(output[0],dpi=300)        


rule substition_analysis:
    input:
        paf = f'{preproc_d}/minimap2/{{sample}}.cDNA.paf',
    output:
        tsv = f'{plots_d}/{{sample}}/substition_analysis.tsv',
    run:
        cigar_re = re.compile(r'(\d+)([M|I|D|N|S|H|P|=|X]{1})')
        with open(input.paf, 'r') as hand, open(output.tsv, 'w') as whand:
            for line in tqdm(hand):
                if "tp:A:P" not in line:
                    continue
                fields = line.rstrip('\n').split('\t')
                for f in fields[11:]:
                    if "cg" in f:
                        cigars = [(int(x),y) for x,y in cigar_re.findall(f[5:])]
                        break
                else:
                    continue
                match_c = np.sum([x for x,y in cigars if y == '='])
                subs_c = np.sum([x for x,y in cigars if y == 'X'])
                insert_c = np.sum([x for x,y in cigars if y == 'I'])
                del_c = np.sum([x for x,y in cigars if y == 'D'])
                print(match_c, subs_c, insert_c, del_c, int(fields[3])-int(fields[2]), sep="\t", file=whand)


rule plot_substition_analysis:
    input:
        tsvs = [f'{plots_d}/{s}/substition_analysis.tsv' for s in config['samples'].keys()] ,
    output:
        f'{plots_d}/all_error_dist.png'
    params:
        bins=50
    run:
        samples = list(config['samples'].keys())
        fig, axs = plt.subplots(4)
        fig.tight_layout()

        bb = [None,None,None,None]
        maxi = [0,0,0,0]
        print('Reading TSV files...')
        for tsv,sample in tqdm(zip(input.tsvs, samples), total=len(input.tsvs)):
            counts = list()
            for line in tqdm(open(tsv, 'r')):
                match_c, subs_c, insert_c, del_c, length = [int(float(x)) for x in line.rstrip().split("\t")]
                l100 = 100/length
                counts.append((
                    l100*match_c,
                    l100*subs_c,
                    l100*insert_c,
                    l100*del_c,
                )) 
            for j,title in enumerate(["match counts", "substition counts", "insertion counts", "deletion counts"]):
                cc = np.array([x[j] for x in counts])
                try:
                    perc = np.percentile(cc, 99)
                except:
                    print(counts)
                    raise
                cc = cc[cc < perc]
                _, bb[j], _ = axs[j].hist(cc, bins=params.bins if bb[j] is None else bb[j], alpha=0.3, label = sample, density=True)
        for i,title in enumerate(["match counts", "substition counts", "insertion counts", "deletion counts"]):
            axs[i].legend()
            axs[i].set_title(title + " per 100 bases")
        plt.savefig(output[0], dpi=300)

rule plot_polyA:
    input:
        fastqs = [config['samples'][s]['fastq'] for s in config['samples']],
    output:
        png = f'{plots_d}/all_polyA.png',
    params:
        bins=30,
    threads:
        32
    run:
        fig = plt.figure()
        bb = None
        max_up = 0
        samples = list(config['samples'].keys())
        print('Reading FASTQ/A files')
        polys_arrays = list()
        for fastq in tqdm(input.fastqs):
            seqs = list()
            if fastq.endswith('.gz'):
                infile = gzip.open(fastq, 'rt')
            else:
                infile = open(fastq, 'rt')
            for idx,line in tqdm(enumerate(infile)):
                if idx == 0:
                    mod = 4 if line[0]=='@' else 2
                if idx % mod != 1:
                    continue
                seqs.append(line.rstrip()) 
            if threads > 1:
                p = Pool(threads)
                mapper = functools.partial(p.imap_unordered, chunksize=100)
            else:
                mapper = map
            polys = list()
            for l in tqdm(mapper(run_longest_polys, seqs), total=len(seqs)):
                polys.append(l)
            if threads > 1:
                p.close()
            polys_arrays.append(np.array(polys))
            max_up = max(max_up, np.percentile(polys, 99))
        for polys,sample in zip(polys_arrays, samples):
            polys = polys[polys < max_up]
            _, bb, _ = plt.hist(
                polys,
                bins = params.bins if bb is None else bb,
                density = True,
                label = sample,
                alpha = .5,
            )
        plt.legend()
        plt.title("Poly(A,T) length distribution")
        plt.savefig(output[0],dpi=300)