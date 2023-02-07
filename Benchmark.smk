import sys
import gzip
from multiprocessing import Pool

import re
import itertools
import functools

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from snakemake.utils import min_version
from tqdm import tqdm

min_version('6.0')

if len(config)==0:
    configfile: 'Benchmark_config.yaml'

module RI_smk:
    snakefile: 'Snakefile'
    config: config

outpath   = config['outpath']
preproc_d = f'{outpath}/preprocess'
RI_d      = f'{outpath}/RI'
NS_d = f'{outpath}/NS'
plots_d   = f'{outpath}/plots'

use rule * from RI_smk as RI_*

use rule minimap_cdna from RI_smk as RI_minimap_cdna with:
    input:
        reads = lambda wc: get_sample_fastqs(wc.sample),
        ref = lambda wildcards: config['refs'][get_sample_ref(wildcards.sample)]['cDNA'],
use rule scTagger_lr_seg from RI_smk as RI_scTagger_lr_seg with:
    input:
        reads = lambda wc: get_sample_fastqs(wc.sample),

def longest_polys(seq, s, e, step, match_score=1, mismatch_score=-2, char='A'):
    if e-s == 0:
        return
    if seq[s] == char:
        scores = [match_score]
    else:
        scores = [0]
    for m in (match_score if c == char else mismatch_score for c in seq[s+step:e:step]):
        scores.append(max(0, scores[-1]+m))
    for k, g in itertools.groupby(enumerate(scores), lambda x: x[1] > 0):
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

def get_sample_ref(name):
    if name in config['samples']:
        return config['samples'][name]['ref']
    elif name in config['RI_experiments']:
        return get_sample_ref(config['RI_experiments'][name]['sample'])
    elif name in config['NS_experiments']:
        return get_sample_ref(config['NS_experiments'][name]['sample'])
    else:
        print(f'Invalid experiment/sample name! {name}')
        1/0

def get_sample_fastqs(name):
    if name in config['samples']:
        sample = name
        return config['samples'][sample]['fastq']
    elif name in config['RI_experiments']:
        exprmnt = name
        return [f'{RI_d}/{exprmnt}/{RI_smk.experiment_prefix(exprmnt)}.fastq']
    elif name in config['NS_experiments']:
        exprmnt = name
        return [f'{NS_d}/{exprmnt}/simulation_aligned_reads.fasta']
    else:
        print(f'Invalid experiment/sample name! {name}')
        1/0

def NS_exprmnt_sample(exprmnt):
    return config['NS_experiments'][exprmnt]['sample']

rule all:
    input:
        [
            f'{plots_d}/{r}/{".".join(p["data"])}.png'
            for p in config['plots']
                for r in p['rules']
        ]
    default_target: True

rule raw_lengths:
    input:
        fastqs = lambda wc: itertools.chain.from_iterable([
            get_sample_fastqs(s) for s in wc.samples.split('.')
        ]),
    output:
        png = f'{plots_d}/raw_lengths/{{samples}}.png',
    params:
        bins=50
    run:
        fig = plt.figure()
        bb = None
        lens_arrays = list()
        max_up = 0
        samples = wildcards.samples.split('.')
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
        plt.title('Length distribution (whole reads)')
        plt.savefig(output[0],dpi=300)        

rule mapped_raw_lengths:
    input:
        pafs = lambda wc: [f'{preproc_d}/minimap2/{s}.cDNA.paf' for s in wc.samples.split('.')],
    output:
        png = f'{plots_d}/mapped_raw_lengths/{{samples}}.png',        
    params:
        bins=50
    run:
        fig = plt.figure()
        bb = None
        lens_arrays = list()
        max_up = 0
        rids = set()
        samples = wildcards.samples.split('.')
        print('Reading PAF files')
        for paf in tqdm(input.pafs):
            lens = list()
            for line in tqdm(open(paf)):
                    if 'tp:A:P' not in line:
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
        plt.title('Length distribution (whole reads from PAF)')
        plt.savefig(output[0],dpi=300)        

rule mapped_lengths:
    input:
        pafs = lambda wc: [f'{preproc_d}/minimap2/{s}.cDNA.paf' for s in wc.samples.split('.')],
    output:
        png = f'{plots_d}/mapped_lengths/{{samples}}.png',        
    params:
        bins=50
    run:
        fig = plt.figure()
        bb = None
        lens_arrays = list()
        max_up = 0
        rids = set()
        samples = wildcards.samples.split('.')
        print('Reading PAF files')
        for paf in tqdm(input.pafs):
            lens = list()
            for line in tqdm(open(paf)):
                    if 'tp:A:P' not in line:
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
        plt.title('Length distribution (mapping part of reads from PAF)')
        plt.savefig(output[0],dpi=300)        


rule substition_stats:
    input:
        paf = f'{preproc_d}/minimap2/{{sample}}.cDNA.paf',
    output:
        tsv = f'{plots_d}/substition_stats/{{sample}}.substition_stats.tsv',
    run:
        cigar_re = re.compile(r'(\d+)([M|I|D|N|S|H|P|=|X]{1})')
        with open(input.paf, 'r') as hand, open(output.tsv, 'w') as whand:
            for line in tqdm(hand):
                if 'tp:A:P' not in line:
                    continue
                fields = line.rstrip('\n').split('\t')
                for f in fields[11:]:
                    if 'cg' in f:
                        cigars = [(int(x),y) for x,y in cigar_re.findall(f[5:])]
                        break
                else:
                    continue
                match_c = np.sum([x for x,y in cigars if y == '='])
                subs_c = np.sum([x for x,y in cigars if y == 'X'])
                insert_c = np.sum([x for x,y in cigars if y == 'I'])
                del_c = np.sum([x for x,y in cigars if y == 'D'])
                print(match_c, subs_c, insert_c, del_c, int(fields[3])-int(fields[2]), sep='\t', file=whand)


rule substition:
    input:
        tsvs = lambda wc: [
            f'{plots_d}/substition_stats/{s}.substition_stats.tsv' 
            for s in wc.samples.split('.')
        ],
    output:
        f'{plots_d}/substition/{{samples}}.png'
    params:
        bins=50
    run:
        samples = wildcards.samples.split('.')
        fig, axs = plt.subplots(4)
        fig.tight_layout()

        bb = [None,None,None,None]
        maxi = [0,0,0,0]
        print('Reading TSV files...')
        for tsv,sample in tqdm(zip(input.tsvs, samples), total=len(input.tsvs)):
            counts = list()
            for line in tqdm(open(tsv, 'r')):
                match_c, subs_c, insert_c, del_c, length = [
                    int(float(x)) 
                    for x in line.rstrip().split('\t')
                ]
                l100 = 100/length
                counts.append((
                    l100*match_c,
                    l100*subs_c,
                    l100*insert_c,
                    l100*del_c,
                )) 
            for j,title in enumerate([
                'match counts',
                'substition counts',
                'insertion counts',
                'deletion counts'
            ]):
                cc = np.array([x[j] for x in counts])
                try:
                    perc = np.percentile(cc, 99)
                except:
                    print(counts)
                    raise
                cc = cc[cc < perc]
                _, bb[j], _ = axs[j].hist(
                    cc,
                    bins    = params.bins if bb[j] is None else bb[j],
                    alpha   = 0.3,
                    label   = sample,
                    density = True,
                )
        for i,title in enumerate([
            'match counts',
            'substition counts',
            'insertion counts',
            'deletion counts'
        ]):
            axs[i].legend()
            axs[i].set_title(f'{title} per 100 bases')
        plt.savefig(output[0], dpi=300)

rule polyA:
    input:
        fastqs = lambda wc: itertools.chain.from_iterable([
            get_sample_fastqs(s) for s in wc.samples.split('.')
        ]),
    output:
        png = f'{plots_d}/polyA/{{samples}}.png',
    params:
        bins=30,
    threads:
        32
    run:
        fig = plt.figure()
        bb = None
        max_up = 0
        samples = wildcards.samples.split('.')
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
        plt.title('Poly(A,T) length distribution')
        plt.savefig(output[0],dpi=300)

rule expression:
    input:
        script = config['exec']['transcript_abundance'],
        paf = f'{preproc_d}/minimap2/{{sample}}.cDNA.paf',
    output:
        tsv = f'{plots_d}/expression_stats/{{sample}}.expression_stats.tsv',
    shell:
       'python {input.script} -p {input.paf} -o {output.tsv}' 

rule tpm_plot:
    input:
        tsvs = lambda wc: [
            f'{plots_d}/expression_stats/{s}.expression_stats.tsv'
            for s in wc.samples.split('.')
        ],
    output:
        png = f'{plots_d}/tpm_plot/{{samples}}.png',        
    run:

        frames = dict()
        real = pd.read_csv(input.tsvs[0], sep="\t", header=0)
        real_keys = set(real['target_id'])
        real.set_index('target_id', inplace=True)

        samples = wildcards.samples.split('.')
        fl = len(samples)-1
        fig,axs = plt.subplots(fl,1,sharex=True,sharey=True,figsize=(12,12*fl))
        for i, (k,v) in enumerate(zip(samples[1:], input.tsvs[1:])):
            ax = axs[i]
            frame = pd.read_csv(v, sep="\t", header=0)
            keys = real_keys.intersection(set(frame['target_id']))
            frame.set_index('target_id',inplace=True)

            X = frame.loc[list(keys)]['tpm']
            Y = real.loc[list(keys)]['tpm']

            ax.plot(X,Y, 'o', alpha=.5, label = k)

            ax.set_xlabel('Simulated read TPM')
            ax.set_ylabel('Simulation input reads TPM')
            ax.legend()
            ax.set_xscale('log')
            ax.set_yscale('log')
        plt.savefig(output[0])

rule NS_analysis:
    input:
        reads = lambda wildcards: config['samples'][wildcards.sample]['fastq'],
        dna  = lambda wildcards: config['refs'][get_sample_ref(wildcards.sample)]['DNA'],
        cdna = lambda wildcards: config['refs'][get_sample_ref(wildcards.sample)]['cDNA'],
        gtf  = lambda wildcards: config['refs'][get_sample_ref(wildcards.sample)]['GTF'],
    output:
        analysis_dir = directory(f'{NS_d}/{{sample}}/analysis'),
    params:
        output_prefix = f'{NS_d}/{{sample}}/analysis/sim',
    threads:
        32
    shell:
        'read_analysis.py transcriptome'
        ' -i {input.reads}'
        ' -rg {input.dna}'
        ' -rt {input.cdna}'
        ' --annotation {input.gtf}'
        ' -o {params.output_prefix}'
        ' -t {threads}'
        ' --no_intron_retention'

rule NS_quantify:
    input:
        reads = lambda wildcards: config['samples'][wildcards.sample]['fastq'],
    output:
        quantify_tsv = f'{NS_d}/{{sample}}/abundance/sim_transcriptome_quantification.tsv',
    params:
        dna  = lambda wildcards: config['refs'][get_sample_ref(wildcards.sample)]['DNA'],
        cdna = lambda wildcards: config['refs'][get_sample_ref(wildcards.sample)]['cDNA'],
        gtf  = lambda wildcards: config['refs'][get_sample_ref(wildcards.sample)]['GTF'],
        output_prefix = f'{NS_d}/{{sample}}/abundance/sim',
    threads:
        32
    shell:
        'read_analysis.py quantify'
        ' -i {input.reads} '
        ' -rt {params.cdna}'
        ' -o {params.output_prefix}'
        ' -t {threads}'
        ' -e trans'

rule NS_simulate:
    input:
        analysis_dir = lambda wc: f'{NS_d}/{NS_exprmnt_sample(wc.exprmnt)}/analysis',
        quantify_tsv = lambda wc: f'{NS_d}/{NS_exprmnt_sample(wc.exprmnt)}/abundance/sim_transcriptome_quantification.tsv',
        dna  = lambda wildcards: config['refs'][get_sample_ref(wildcards.exprmnt)]['DNA'],
        cdna = lambda wildcards: config['refs'][get_sample_ref(wildcards.exprmnt)]['cDNA'],
    output:
        fasta = f'{NS_d}/{{exprmnt}}/simulation_aligned_reads.fasta',
    params:
        analysis_model = lambda wc: f'{NS_d}/{NS_exprmnt_sample(wc.exprmnt)}/analysis/sim',
        out_prefix = lambda wc: f'{NS_d}/{wc.exprmnt}/simulation',
        other = lambda wc: config['NS_experiments'][wc.exprmnt]['simulation_params'],
    threads:
        32
    shell:
        'simulator.py transcriptome'
        ' -rt {input.cdna}'
        ' -rg {input.dna}'
        ' --exp {input.quantify_tsv}'
        ' -c {params.analysis_model}'
        ' -o {params.out_prefix}'
        ' -t {threads}'
        ' --no_model_ir'
        ' {params.other}'