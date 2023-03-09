import sys
import gzip
from multiprocessing import Pool

import re
import itertools
import functools
from collections import Counter

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib 
from snakemake.utils import min_version
from tqdm import tqdm
from sklearn.metrics import r2_score,mean_squared_error

min_version('6.0')

if len(config)==0:
    configfile: 'Benchmark_config.yaml'

module TS_smk:
    snakefile: 'Snakefile'
    config: config

outpath   = config['outpath']
preproc_d = f'{outpath}/preprocess'
TS_d      = f'{outpath}/TS'
NS_d = f'{outpath}/NS'
time_d = f'{outpath}/time'
plots_d   = f'{outpath}/plots'

use rule * from TS_smk as TS_*

use rule minimap_cdna from TS_smk as TS_minimap_cdna with:
    input:
        reads = lambda wc: get_sample_fastqs(wc.sample),
        ref = lambda wildcards: config['refs'][get_sample_ref(wildcards.sample)]['cDNA'],
use rule scTagger_lr_seg from TS_smk as TS_scTagger_lr_seg with:
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
            # if p > .8:
            L = max(L,l)
    return L

def get_sample_ref(name):
    if name in config['samples']:
        return config['samples'][name]['ref']
    elif name in config['TS_experiments']:
        return get_sample_ref(config['TS_experiments'][name]['sample'])
    elif name in config['NS_experiments']:
        return get_sample_ref(config['NS_experiments'][name]['sample'])
    else:
        print(f'Invalid experiment/sample name! {name}')
        1/0

def get_sample_fastqs(name):
    if name in config['samples']:
        sample = name
        return config['samples'][sample]['fastq']
    elif name in config['TS_experiments']:
        exprmnt = name
        return [f'{TS_d}/{exprmnt}/{TS_smk.experiment_prefix(exprmnt)}.fastq']
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
        for fastq in tqdm(input.fastqs, total=len(input.fastqs), desc='[raw_lengths] Reading FASTQ/A'):
            lens = list()
            if fastq.endswith('.gz'):
                infile = gzip.open(fastq, 'rt')
            else:
                infile = open(fastq, 'rt')
            for idx,line in tqdm(enumerate(infile), desc=f'[raw_lengths] Processing {fastq}'):
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
            print(f'{sample} length mean: {np.mean(lens):.2f}')
            print(f'{sample} length std: {np.std(lens):.2f}')
            print(f'{sample} length median: {np.median(lens):.2f}')
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
        for paf in tqdm(input.pafs, total=len(input.pafs), desc='[mapped_raw_lengths] Reading PAF'):
            lens = list()
            for line in tqdm(open(paf), desc=f'[mapped_raw_lengths] Processing {paf}'):
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
            print(f'{sample} length (PAF) mean: {np.mean(lens):.2f}')
            print(f'{sample} length (PAF) std: {np.std(lens):.2f}')
            print(f'{sample} length (PAF) median: {np.median(lens):.2f}')
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
        for paf in tqdm(input.pafs, total=len(input.pafs), desc='[mapped_lengths] Reading PAF'):
            lens = list()
            for line in tqdm(open(paf), desc=f'[mapped_lengths] Processing {paf}'):
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
            print(f'{sample} mapping length (PAF) mean: {np.mean(lens):.2f}')
            print(f'{sample} mapping length (PAF) std: {np.std(lens):.2f}')
            print(f'{sample} mapping length (PAF) median: {np.median(lens):.2f}')
        plt.xlim(left=0, right=max_up)
        plt.legend()
        plt.title('Length distribution (mapping part of reads from PAF)')
        plt.savefig(output[0],dpi=300)        


rule substitution_stats:
    input:
        paf = f'{preproc_d}/minimap2/{{sample}}.cDNA.paf',
    output:
        tsv = f'{plots_d}/substitution_stats/{{sample}}.substitution_stats.tsv',
    run:
        cigar_re = re.compile(r'(\d+)([M|I|D|N|S|H|P|=|X]{1})')
        with open(input.paf, 'r') as hand, open(output.tsv, 'w') as whand:
            for line in tqdm(hand, desc=f'[substitution_stats] Processing {input.paf}'):
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


rule substitution:
    input:
        tsvs = lambda wc: [
            f'{plots_d}/substitution_stats/{s}.substitution_stats.tsv' 
            for s in wc.samples.split('.')
        ],
    output:
        f'{plots_d}/substitution/{{samples}}.png'
    params:
        bins=50
    run:
        samples = wildcards.samples.split('.')
        fig, axs = plt.subplots(2,2, figsize=(12,4))
        axs = axs.flatten()
        bb = [None,None,None,None]
        maxi = [0,0,0,0]
        print('Reading TSV files...')
        for tsv,sample in tqdm(zip(input.tsvs, samples), total=len(input.tsvs), desc='[substitution] Reading TSV files'):
            counts = list()
            for line in tqdm(open(tsv, 'r'), desc=f'[substitution] Processing {tsv}'):
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
                'substitution counts',
                'insertion counts',
                'deletion counts'
            ]):
                cc = np.array([x[j] for x in counts])
                _, bb[j], _ = axs[j].hist(
                    cc,
                    bins    = params.bins if bb[j] is None else bb[j],
                    alpha   = 0.3,
                    label   = f'{sample}: (µ={np.mean(cc):.2f}, σ={np.std(cc):.2f})',
                    density = True,
                )
                axs[j].legend()
        for i,title in enumerate([
            'match counts',
            'substitution counts',
            'insertion counts',
            'deletion counts'
        ]):
            axs[i].set_title(f'{title} per 100 bases')
        fig.tight_layout()
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
        for fastq in tqdm(input.fastqs, total=len(input.fastqs), desc='[polyA] Reading FASTQ/A files'):
            seqs = list()
            if fastq.endswith('.gz'):
                infile = gzip.open(fastq, 'rt')
            else:
                infile = open(fastq, 'rt')
            for idx,line in tqdm(enumerate(infile), desc=f'[polyA] Processing {fastq}'):
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
            for l in tqdm(mapper(run_longest_polys, seqs), total=len(seqs), desc=f'[polyA] Finding longest tails in {fastq}'):
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
                label = f'{sample} (µ={np.mean(polys):.2f}, σ={np.std(polys):.2f})',
                alpha = .5,
            )
        plt.legend()
        plt.title('Poly(A,T) length distribution')
        plt.savefig(output[0],dpi=300)

rule abundance:
    input:
        binary = config['exec']['tksm'],
        paf = f'{preproc_d}/minimap2/{{sample}}.cDNA.paf',
    output:
        tsv = f'{plots_d}/expression_stats/{{sample}}.expression_stats.tsv',
    shell:
        '{input.binary} abundance'
        ' -p {input.paf}'
        ' -o {output.tsv}' 

rule LIQA_refgene:
    input:
        gtf = '{prefix}.gtf',
    output:
        refgen = '{prefix}.gtf.refgene',
    shell:
        'python3 extern/LIQA/liqa_src/liqa.py'
        ' -task refgene'
        ' -format gtf'
        ' -ref {input.gtf}'
        ' -out {output.refgen}'

rule minimap_dna:
    input:
        reads = lambda wc: get_sample_fastqs(wc.sample),
        ref = lambda wildcards: config['refs'][get_sample_ref(wildcards.sample)]['DNA'],
    output:
        bam = f'{preproc_d}/minimap2/{{sample}}.DNA.bam',
        bai = f'{preproc_d}/minimap2/{{sample}}.DNA.bam.bai',
    benchmark:
        f'{time_d}/{{sample}}/minimap2_cdna.benchmark'
    threads:
        32
    shell:
        'minimap2'
        ' -t {threads}'
        ' -x splice'
        ' -a'
        ' {input.ref}'
        ' {input.reads}'
        ' | samtools sort '
        ' -T {output.bam}.tmp'
        ' -@ {threads}'
        ' -o {output.bam}'
        ' - '
        ' && samtools index {output.bam}'


rule LIQA_quantify:
    input:
        refgen = lambda wildcards: f"{config['refs'][get_sample_ref(wildcards.sample)]['GTF']}.refgene",
        bam = f'{preproc_d}/minimap2/{{sample}}.DNA.bam',
    output:
        tsv = f'{plots_d}/expression_stats/{{sample}}.LIQA.tsv',
    threads:
        16
    shell:
        'python3 extern/LIQA/liqa_src/liqa.py'
        ' -task quantify'
        ' -refgene {input.refgen}'
        ' -bam {input.bam}'
        ' -out {output.tsv}'
        ' -max_distance 20'
        ' -f_weight 1'
        ' -threads {threads} > /dev/null'

rule paf_quantify:
    input:
        paf = f'{preproc_d}/minimap2/{{sample}}.cDNA.paf',
    output:
        tsv = f'{plots_d}/expression_stats/{{sample}}.PAF.tsv',
    run:
        tid_to_count = Counter()
        for line in tqdm(open(input.paf), desc=f'[paf_quantify] Processing {input.paf}'):
            if 'tp:A:P' not in line:
                continue
            fields = line.rstrip('\n').split('\t')
            tid = fields[5]
            tid_to_count[tid] += 1
        with open(output.tsv, 'w') as f:
            f.write('transcript_id\tcount\n')
            for tid,count in tid_to_count.items():
                f.write(f'{tid}\t{count}\n')

rule nanosim_quantify:
    input:
        paf = f'{preproc_d}/minimap2/{{sample}}.cDNA.paf',
    output:
        tsv = f'{plots_d}/expression_stats/{{sample}}.nanosim.tsv',
    shell:
        'nanopore_transcript_abundance.py'
        ' -i {input.paf}'
        ' > {output.tsv}'

rule tpm_plot_nanosim:
    input:
        tsvs = lambda wc: [
            f'{NS_d}/{s}/abundance/sim_transcriptome_quantification.tsv'
            for s in wc.samples.split('.')
        ],
    output:
        png = f'{plots_d}/tpm_plot_nanosim/{{samples}}.png',        
    run:
        my_get_tpm = functools.partial(get_tpm, key_col=0, val_col=2, header=True)
        X_tpm = my_get_tpm(input.tsvs[0])
        Y_tpms = list(map(my_get_tpm, input.tsvs[1:]))
        samples = wildcards.samples.split('.')
        plot_tpm_func(X_tpm, Y_tpms, samples, output.png, title='TPM using Nanosim abundance estimates')


rule tpm_plot_paf:
    input:
        tsvs = lambda wc: [
            f'{plots_d}/expression_stats/{s}.PAF.tsv'
            for s in wc.samples.split('.')
        ],
    output:
        png = f'{plots_d}/tpm_plot_paf/{{samples}}.png',        
    run:
        my_get_tpm = functools.partial(get_tpm, key_col=0, val_col=1, header=True)
        X_tpm = my_get_tpm(input.tsvs[0])
        Y_tpms = list(map(my_get_tpm, input.tsvs[1:]))
        samples = wildcards.samples.split('.')
        plot_tpm_func(X_tpm, Y_tpms, samples, output.png, title='TPM using minimap2 primary alignments')

def get_tpm(tsv, min_tpm=0.001, min_val=0.0, key_col=0, val_col=1, header=True):
    tpm = Counter()
    for line_num,line in enumerate(open(tsv, 'r')):
        if header and line_num == 0:
            continue
        line = line.rstrip().split('\t')
        key = line[key_col]
        val = float(line[val_col])
        if val>min_val:
            tpm[key] = val
    thruput = sum(tpm.values())
    for k,v in tpm.items():
        tpm[k] = v / thruput * 1e6
    tpm = Counter({k:v for k,v in tpm.items() if v>min_tpm})
    thruput = sum(tpm.values())
    for k,v in tpm.items():
        tpm[k] = v / thruput * 1e6
    return tpm

def plot_tpm_func(X_tpm, Y_tpms, samples, outpath, title):
    plt.rc('font', size=22)

    plot_count = len(samples)-1
    fig,axs = plt.subplots(
        plot_count,
        3,
        sharex=True,
        sharey=True,
        figsize=(10*3,10*plot_count),
    )
    fig.suptitle(title, fontsize=30)
    axs[0,0].set_title('Transcripts in sample AND in input', fontsize=28)
    axs[0,1].set_title('Transcripts in sample OR in input', fontsize=28)
    axs[0,2].set_title('Transcripts in input (regardless of sample)', fontsize=28)

    for sample,ax in zip(samples[1:],axs[:,0]):
        ax.set_ylabel(f'TPM in {sample}', fontsize=28)
    for ax in axs[-1]:
        ax.set_xlabel('TPM in input', fontsize=28)
    X_tids = set(X_tpm.keys())
    for idx, (sample, Y_tpm) in enumerate(zip(samples[1:], Y_tpms)):
        Y_tids = set(Y_tpm.keys())

        ax = axs[idx,0]
        for ax, select_tids in zip(axs[idx], [
            X_tids & Y_tids,
            X_tids | Y_tids,
            X_tids,
        ]):
            X = np.array([X_tpm[tid] for tid in select_tids])
            Y = np.array([Y_tpm[tid] for tid in select_tids])
            ax.plot(X, Y, 'o', alpha=.5, label = sample)
            ax.text(
                .15,
                .75,
                f"N = {len(select_tids)}\n" + \
                f"Sample: {sample}\n" + \
                f"r-squared = {r2_score(X,Y):.3f}\n" + \
                f"RMSE = {mean_squared_error(X,Y, squared=False):.1f}",
                transform=ax.transAxes,
            )
            ax.set_xscale('log')
            ax.set_yscale('log')
    fig.tight_layout()            
    plt.savefig(outpath)

rule tpm_plot_liqa:
    input:
        tsvs = lambda wc: [
            f'{plots_d}/expression_stats/{s}.LIQA.tsv'
            for s in wc.samples.split('.')
        ],
    output:
        png = f'{plots_d}/tpm_plot_liqa/{{samples}}.png',        
    run:
        my_get_tpm = functools.partial(get_tpm, min_val=0.99, key_col=1, val_col=2, header=True)
        X_tpm = my_get_tpm(input.tsvs[0])
        Y_tpms = list(map(my_get_tpm, input.tsvs[1:]))
        samples = wildcards.samples.split('.')
        plot_tpm_func(X_tpm, Y_tpms, samples, output.png, 'TPM using LIQA')

rule tpm_plot_tksm:
    input:
        tsvs = lambda wc: [
            f'{plots_d}/expression_stats/{s}.expression_stats.tsv'
            for s in wc.samples.split('.')
        ],
    output:
        png = f'{plots_d}/tpm_plot_tksm/{{samples}}.png',        
    run:
        my_get_tpm = functools.partial(get_tpm, key_col=0, val_col=1, header=True)
        X_tpm = my_get_tpm(input.tsvs[0])
        Y_tpms = list(map(my_get_tpm, input.tsvs[1:]))
        samples = wildcards.samples.split('.')
        plot_tpm_func(X_tpm, Y_tpms, samples, output.png, title='TPM using TKSM abundance estimates')

rule NS_analysis:
    input:
        reads = lambda wildcards: config['samples'][wildcards.sample]['fastq'],
        dna  = lambda wildcards: config['refs'][get_sample_ref(wildcards.sample)]['DNA'],
        cdna = lambda wildcards: config['refs'][get_sample_ref(wildcards.sample)]['cDNA'],
        gtf  = lambda wildcards: config['refs'][get_sample_ref(wildcards.sample)]['GTF'],
    output:
        analysis_dir = directory(f'{NS_d}/{{sample}}/analysis'),
    benchmark:
        f'{time_d}/{{sample}}/NS_analysis.benchmark',
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
        reads = lambda wildcards: get_sample_fastqs(wildcards.sample),
    output:
        quantify_tsv = f'{NS_d}/{{sample}}/abundance/sim_transcriptome_quantification.tsv',
    benchmark:
        f'{time_d}/{{sample}}/NS_quantify.benchmark',
    params:
        cdna = lambda wildcards: config['refs'][get_sample_ref(wildcards.sample)]['cDNA'],
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
    benchmark:
        f'{time_d}/{{exprmnt}}/NS_simulate.benchmark',
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