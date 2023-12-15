import sys
import gzip
from multiprocessing import Pool
import pickle
import glob

import re
import itertools
import functools
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
from snakemake.utils import min_version
from tqdm import tqdm
from sklearn.metrics import r2_score, mean_squared_error

min_version("6.0")

if len(config) == 0:

    configfile: "Benchmark_config.yaml"


module TS_smk:
    snakefile:
        "Snakefile"
    config:
        config


outpath = config["outpath"]
preproc_d = f"{outpath}/preprocess"
TS_d = f"{outpath}/TS"
NS_d = f"{outpath}/NS"
plots_d = f"{outpath}/plots"


### Import TKSM Snakemake rules
use rule * from TS_smk as TS_*


use rule minimap_cdna from TS_smk as TS_minimap_cdna with:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "cDNA"),
        time=ancient(TS_smk.time_tsv),


use rule scTagger_lr_seg from TS_smk as TS_scTagger_lr_seg with:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        time=ancient(TS_smk.time_tsv),


def get_sample_ref_names(sample):
    try:
        # Try should work for TSKM and real samples
        return TS_smk.get_sample_ref_names(sample)
    except ValueError:
        if sample in config["NS_experiments"]:
            return TS_smk.get_sample_ref_names(
                config["NS_experiments"][sample]["sample"]
            )
        raise ValueError(f"Invalid sample name! {sample}")


def get_sample_ref(sample, ref_type):
    ref_names = get_sample_ref_names(sample)
    ref_name = ":".join(ref_names)
    if ref_type in ["DNA", "cDNA"]:
        file_type = "fasta"
    elif ref_type == "GTF":
        file_type = "gtf"
    else:
        raise ValueError(f"Invalid reference type! {ref_type}")
    return f"{preproc_d}/refs/{ref_name}.{ref_type}.{file_type}"


def get_sample_fastqs(name):
    try:
        return TS_smk.get_sample_fastqs(name)
    except ValueError:
        if name in config["NS_experiments"]:
            exprmnt = name
            return [f"{NS_d}/{exprmnt}/simulation_aligned_reads.fasta"]
        raise ValueError(f"Invalid experiment/sample name! {name}")


def NS_exprmnt_sample(exprmnt):
    return config["NS_experiments"][exprmnt]["sample"]


### Trans-Nanosim rules
rule NS_analysis:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        dna=lambda wc: get_sample_ref(wc.sample, "DNA"),
        cdna=lambda wc: get_sample_ref(wc.sample, "cDNA"),
        gtf=lambda wc: get_sample_ref(wc.sample, "GTF"),
        time=ancient(TS_smk.time_tsv),
    output:
        analysis_dir=directory(f"{NS_d}/{{sample}}/analysis"),
    params:
        output_prefix=f"{NS_d}/{{sample}}/analysis/sim",
    threads: 32
    conda:
        "Benchmark_env.yaml"
    shell:
        f"{TS_smk.format_gnu_time_string(process='NS_analysis', exprmnt='{wildcards.sample}', prefix='')}"
        "read_analysis.py transcriptome"
        " -i {input.reads}"
        " -rg {input.dna}"
        " -rt {input.cdna}"
        " --annotation {input.gtf}"
        " -o {params.output_prefix}"
        " -t {threads}"
        " --no_intron_retention"


rule NS_quantify:
    input:
        reads=lambda wildcards: get_sample_fastqs(wildcards.sample),
        cdna=lambda wc: get_sample_ref(wc.sample, "cDNA"),
        time=ancient(TS_smk.time_tsv),
    output:
        quantify_tsv=f"{NS_d}/{{sample}}/abundance/sim_transcriptome_quantification.tsv",
    params:
        output_prefix=f"{NS_d}/{{sample}}/abundance/sim",
    threads: 32
    resources:
        time=60 * 6 - 1,
    conda:
        "Benchmark_env.yaml"
    shell:
        f"{TS_smk.format_gnu_time_string(process='NS_quantify', exprmnt='{wildcards.sample}', prefix='')}"
        "read_analysis.py quantify"
        " -i {input.reads} "
        " -rt {input.cdna}"
        " -o {params.output_prefix}"
        " -t {threads}"
        " -e trans"


rule NS_simulate:
    input:
        analysis_dir=lambda wc: f"{NS_d}/{NS_exprmnt_sample(wc.exprmnt)}/analysis",
        quantify_tsv=lambda wc: f"{NS_d}/{NS_exprmnt_sample(wc.exprmnt)}/abundance/sim_transcriptome_quantification.tsv",
        dna=lambda wc: get_sample_ref(wc.exprmnt, "DNA"),
        cdna=lambda wc: get_sample_ref(wc.exprmnt, "cDNA"),
        time=ancient(TS_smk.time_tsv),
    output:
        fasta=f"{NS_d}/{{exprmnt}}/simulation_aligned_reads.fasta",
    params:
        analysis_model=lambda wc: f"{NS_d}/{NS_exprmnt_sample(wc.exprmnt)}/analysis/sim",
        out_prefix=lambda wc: f"{NS_d}/{wc.exprmnt}/simulation",
        other=lambda wc: config["NS_experiments"][wc.exprmnt]["simulation_params"],
    threads: 32
    conda:
        "Benchmark_env.yaml"
    shell:
        f"{TS_smk.format_gnu_time_string(process='NS_simulate', exprmnt='{wildcards.exprmnt}', prefix='')}"
        "simulator.py transcriptome"
        " -rt {input.cdna}"
        " -rg {input.dna}"
        " --exp {input.quantify_tsv}"
        " -c {params.analysis_model}"
        " -o {params.out_prefix}"
        " -t {threads}"
        " --no_model_ir"
        " {params.other}"


### Ploting helper methods
tpm_method_settings = {
    "minimap2": dict(
        get_tpm_args=dict(
            key_col=0,
            val_col=1,
            header=True,
        ),
        title="TPM using minimap2 primary alignments",
    ),
    "nanosim": dict(
        get_tpm_args=dict(
            key_col=0,
            val_col=2,
            header=True,
        ),
        title="TPM using Nanosim abundance estimates",
    ),
    "liqa": dict(
        get_tpm_args=dict(
            key_col=1,
            val_col=2,
            header=True,
        ),
        title="TPM using LIQA abundance estimates",
    ),
    "tksm": dict(
        get_tpm_args=dict(
            key_col=0,
            val_col=1,
            header=True,
        ),
        title="TPM using TKSM abundance estimates",
    ),
    "tksm_sc": dict(
        get_tpm_args=dict(
            key_col=[0, 2],
            val_col=1,
            header=True,
        ),
        title="TPM using TKSM SC abundance estimates",
    ),
}


def longest_polys(seq, s, e, step, match_score=1, mismatch_score=-2, char="A"):
    if e - s == 0:
        return
    if seq[s] == char:
        scores = [match_score]
    else:
        scores = [0]
    for m in (
        match_score if c == char else mismatch_score for c in seq[s + step : e : step]
    ):
        scores.append(max(0, scores[-1] + m))
    for k, g in itertools.groupby(enumerate(scores), lambda x: x[1] > 0):
        if not k:
            continue
        i, S = list(zip(*g))
        max_s, max_i = max(zip(S, i))
        l = max_i + 1 - i[0]
        yield i[0], l, seq[s:e:step][i[0] : i[0] + l].count(char) / l


def run_longest_polys(seq):
    L = 0
    for c in ["A", "T"]:
        for i, l, p in longest_polys(seq, 0, len(seq), 1):
            L = max(L, l)
    return L


def reverse_complement(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(complement.get(base, base) for base in reversed(seq))


def get_rule_extension(r):
    if r.startswith("tpm_plot"):
        return "png"
    if r == "gene_fusion_accuracy":
        return "csv"
    return "pdf"


rule all:
    input:
        [
            f"{plots_d}/{r}/{'.'.join(p['data'])}.{get_rule_extension(r)}"
            for p in config["plots"]
            for r in p["rules"]
        ],
    default_target: True


rule raw_lengths:
    input:
        fastqs=lambda wc: itertools.chain.from_iterable(
            [get_sample_fastqs(s) for s in wc.samples.split(".")]
        ),
    output:
        pdf=f"{plots_d}/raw_lengths/{{samples}}.pdf",
        png=f"{plots_d}/raw_lengths/{{samples}}.png",
    params:
        bins=50,
    run:
        fig = plt.figure()
        bb = None
        lens_arrays = list()
        max_up = 0
        samples = wildcards.samples.split(".")
        for fastq in tqdm(
            input.fastqs, total=len(input.fastqs), desc="[raw_lengths] Reading FASTQ/A"
        ):
            lens = list()
            if fastq.endswith(".gz"):
                infile = gzip.open(fastq, "rt")
            else:
                infile = open(fastq, "rt")
            for idx, line in tqdm(
                enumerate(infile), desc=f"[raw_lengths] Processing {fastq}"
            ):
                if idx == 0:
                    mod = 4 if line[0] == "@" else 2
                if idx % mod == 1:
                    lens.append(len(line) - 1)
            lens = np.array(lens)
            lens_arrays.append(lens)
            max_up = max(max_up, np.percentile(lens, 99))
        for lens, sample in zip(lens_arrays, samples):
            lens = lens[lens < max_up]
            _, bb, _ = plt.hist(
                lens,
                bins=params.bins if bb is None else bb,
                density=True,
                label=sample,
                histtype="step",
            )
            print(f"{sample} length mean: {np.mean(lens):.2f}")
            print(f"{sample} length std: {np.std(lens):.2f}")
            print(f"{sample} length median: {np.median(lens):.2f}")
        plt.legend()
        plt.title("Length distribution (whole reads)")
        plt.savefig(output.pdf, dpi=300)
        plt.savefig(output.png, dpi=1000)


rule mapped_raw_lengths:
    input:
        pafs=lambda wc: [
            f"{preproc_d}/minimap2/{s}.cDNA.paf" for s in wc.samples.split(".")
        ],
    output:
        pdf=f"{plots_d}/mapped_raw_lengths/{{samples}}.pdf",
        png=f"{plots_d}/mapped_raw_lengths/{{samples}}.png",
    params:
        bins=50,
    run:
        fig = plt.figure()
        bb = None
        lens_arrays = list()
        max_up = 0
        rids = set()
        samples = wildcards.samples.split(".")
        print("Reading PAF files")
        for paf in tqdm(
            input.pafs, total=len(input.pafs), desc="[mapped_raw_lengths] Reading PAF"
        ):
            lens = list()
            for line in tqdm(open(paf), desc=f"[mapped_raw_lengths] Processing {paf}"):
                if "tp:A:P" not in line:
                    continue
                fields = line.rstrip("\n").split("\t")
                rid = fields[0]
                if rid in rids:
                    continue
                rids.add(rid)
                qlen = int(fields[1])
                lens.append(qlen)
            lens = np.array(lens)
            lens_arrays.append(lens)
            max_up = max(max_up, np.percentile(lens, 99))

        for lens, sample in zip(lens_arrays, samples):
            lens = lens[lens < max_up]
            _, bb, _ = plt.hist(
                lens,
                bins=params.bins if bb is None else bb,
                density=True,
                label=sample,
                histtype="step",
            )
            print(f"{sample} length (PAF) mean: {np.mean(lens):.2f}")
            print(f"{sample} length (PAF) std: {np.std(lens):.2f}")
            print(f"{sample} length (PAF) median: {np.median(lens):.2f}")
        plt.xlim(left=0, right=max_up)
        plt.legend()
        plt.title("Length distribution (whole reads from PAF)")
        plt.savefig(output.pdf, dpi=300)
        plt.savefig(output.png, dpi=1000)


rule mapped_lengths:
    input:
        pafs=lambda wc: [
            f"{preproc_d}/minimap2/{s}.cDNA.paf" for s in wc.samples.split(".")
        ],
    output:
        pdf=f"{plots_d}/mapped_lengths/{{samples}}.pdf",
        png=f"{plots_d}/mapped_lengths/{{samples}}.png",
    params:
        bins=50,
    run:
        fig = plt.figure()
        bb = None
        lens_arrays = list()
        max_up = 0
        rids = set()
        samples = wildcards.samples.split(".")
        print("Reading PAF files")
        for paf in tqdm(
            input.pafs, total=len(input.pafs), desc="[mapped_lengths] Reading PAF"
        ):
            lens = list()
            for line in tqdm(open(paf), desc=f"[mapped_lengths] Processing {paf}"):
                if "tp:A:P" not in line:
                    continue
                fields = line.rstrip("\n").split("\t")
                rid = fields[0]
                if rid in rids:
                    continue
                rids.add(rid)
                start = int(fields[7])
                end = int(fields[8])
                lens.append(end - start)
            lens = np.array(lens)
            lens_arrays.append(lens)
        max_up = 5_000
        for lens, sample in zip(lens_arrays, samples):
            lens = lens[lens < max_up]
            _, _, _ = plt.hist(
                lens,
                bins=np.arange(0, max_up, 100),
                density=True,
                label=sample,
                histtype="step",
            )
            print(f"{sample} mapping length (PAF) mean: {np.mean(lens):.2f}")
            print(f"{sample} mapping length (PAF) std: {np.std(lens):.2f}")
            print(f"{sample} mapping length (PAF) median: {np.median(lens):.2f}")
        plt.xlim(left=0, right=max_up)
        plt.legend()
        plt.title("Length distribution (mapping part of reads from PAF)")
        plt.savefig(output.pdf, dpi=300)
        plt.savefig(output.png, dpi=1000)


rule mapped_truncated_ratio:
    input:
        pafs=lambda wc: [
            f"{preproc_d}/minimap2/{s}.cDNA.paf" for s in wc.samples.split(".")
        ],
    output:
        pdf=f"{plots_d}/mapped_truncated_ratio/{{samples}}.pdf",
        png=f"{plots_d}/mapped_truncated_ratio/{{samples}}.png",
    params:
        bins=50,
    run:
        end_fracts = list()
        samples = wildcards.samples.split(".")
        print("Reading PAF files")
        for paf in tqdm(
            input.pafs,
            total=len(input.pafs),
            desc="[mapped_truncated_ratio] Reading PAF",
        ):
            cur_end_fracs = list()
            for line in tqdm(
                open(paf), desc=f"[mapped_truncated_ratio] Processing {paf}"
            ):
                if "tp:A:P" not in line:
                    continue
                line = line.rstrip("\n").split("\t")
                strand = line[4]
                tlen = int(line[6])
                start = int(line[7])
                end = int(line[8])
                alen = end - start
                trunc_len = tlen - alen
                if trunc_len == 0:
                    continue
                if strand == "+":
                    end_trunc = tlen - end
                else:
                    end_trunc = start
                cur_end_fracs.append(end_trunc / trunc_len)
            end_fracts.append(cur_end_fracs)
        fig, axes = plt.subplots(2, 1, figsize=(12, 10))
        for sample, V in zip(samples, end_fracts):
            counts, bins = np.histogram(V, bins=np.arange(0, 1.01, 0.01), density=True)
            counts /= 100
            axes[0].stairs(
                counts,
                bins,
                label=sample,
            )
            counts = np.cumsum(counts)
            axes[1].stairs(
                counts,
                bins,
                label=sample,
            )

        for ax in axes:
            ax.set_xlim(left=0, right=1)
            ax.grid(
                linestyle="-",
                alpha=0.5,
            )
            ax.legend(loc="upper center")
        axes[0].set_ylabel("Density")
        axes[1].set_ylabel("Cumulative density")
        fig.suptitle("5' truncation fraction of total truncation length")

        fig.tight_layout()
        plt.savefig(output.pdf, dpi=300)
        plt.savefig(output.png, dpi=1000)


rule substitution_stats:
    input:
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
    output:
        tsv=f"{plots_d}/substitution_stats/{{sample}}.substitution_stats.tsv",
    run:
        cigar_re = re.compile(r"(\d+)([M|I|D|N|S|H|P|=|X]{1})")
        with open(input.paf, "r") as hand, open(output.tsv, "w") as whand:
            for line in tqdm(hand, desc=f"[substitution_stats] Processing {input.paf}"):
                if "tp:A:P" not in line:
                    continue
                fields = line.rstrip("\n").split("\t")
                for f in fields[11:]:
                    if "cg" in f:
                        cigars = [(int(x), y) for x, y in cigar_re.findall(f[5:])]
                        break
                else:
                    continue
                match_c = np.sum([x for x, y in cigars if y == "="])
                subs_c = np.sum([x for x, y in cigars if y == "X"])
                insert_c = np.sum([x for x, y in cigars if y == "I"])
                del_c = np.sum([x for x, y in cigars if y == "D"])
                print(
                    match_c,
                    subs_c,
                    insert_c,
                    del_c,
                    int(fields[3]) - int(fields[2]),
                    sep="\t",
                    file=whand,
                )


rule substitution:
    input:
        tsvs=lambda wc: [
            f"{plots_d}/substitution_stats/{s}.substitution_stats.tsv"
            for s in wc.samples.split(".")
        ],
    output:
        pdf=f"{plots_d}/substitution/{{samples}}.pdf",
        png=f"{plots_d}/substitution/{{samples}}.png",
    params:
        bins=50,
    run:
        samples = wildcards.samples.split(".")
        fig, axs = plt.subplots(2, 2, figsize=(12, 4), sharey=True)
        axs = axs.flatten()
        bb = [None, None, None, None]
        maxi = [0, 0, 0, 0]
        mini = [100, 100, 100, 100]
        titles = [
            "Match",
            "Substitution",
            "Insertion",
            "Deletion",
        ]
        print("Reading TSV files...")
        for tsv, sample in tqdm(
            zip(input.tsvs, samples),
            total=len(input.tsvs),
            desc="[substitution] Reading TSV files",
        ):
            counts = tuple(list() for _ in titles)
            for line in tqdm(open(tsv, "r"), desc=f"[substitution] Processing {tsv}"):
                stats = [int(float(x)) for x in line.rstrip().split("\t")]
                assert len(stats) == 5, line
                l100 = 100 / stats[4]
                for c, s in zip(counts, stats[:4]):
                    c.append(l100 * s)
            counts = tuple(np.array(c) for c in counts)
            for i in range(len(counts)):
                maxi[i] = max(maxi[i], np.percentile(counts[i], 99))
                mini[i] = min(mini[i], np.percentile(counts[i], 1))
            for i, title in enumerate(titles):
                c = counts[i]
                _, bb[i], _ = axs[i].hist(
                    c,
                    bins=params.bins if bb[i] is None else bb[i],
                    alpha=0.3,
                    label=f"{sample}: (µ={np.mean(c):.2f}, σ={np.std(c):.2f})",
                    density=True,
                )
                axs[i].set_xlim(left=mini[i], right=maxi[i])
                axs[i].legend()
                axs[i].set_title(f"{title}")
        fig.tight_layout()
        plt.savefig(output.pdf, dpi=300)
        plt.savefig(output.png, dpi=1000)


rule polyA:
    input:
        fastqs=lambda wc: itertools.chain.from_iterable(
            [get_sample_fastqs(s) for s in wc.samples.split(".")]
        ),
    output:
        pdf=f"{plots_d}/polyA/{{samples}}.pdf",
        png=f"{plots_d}/polyA/{{samples}}.png",
    params:
        bins=30,
    threads: 32
    run:
        fig = plt.figure()
        bb = None
        max_up = 0
        samples = wildcards.samples.split(".")
        print("Reading FASTQ/A files")
        polys_arrays = list()
        for fastq in tqdm(
            input.fastqs, total=len(input.fastqs), desc="[polyA] Reading FASTQ/A files"
        ):
            seqs = list()
            if fastq.endswith(".gz"):
                infile = gzip.open(fastq, "rt")
            else:
                infile = open(fastq, "rt")
            for idx, line in tqdm(
                enumerate(infile), desc=f"[polyA] Processing {fastq}"
            ):
                if idx == 0:
                    mod = 4 if line[0] == "@" else 2
                if idx % mod != 1:
                    continue
                seqs.append(line.rstrip())
            if threads > 1:
                p = Pool(threads)
                mapper = functools.partial(p.imap_unordered, chunksize=100)
            else:
                mapper = map
            polys = list()
            for l in tqdm(
                mapper(run_longest_polys, seqs),
                total=len(seqs),
                desc=f"[polyA] Finding longest tails in {fastq}",
            ):
                polys.append(l)
            if threads > 1:
                p.close()
            polys_arrays.append(np.array(polys))
            max_up = max(max_up, np.percentile(polys, 99))
        for polys, sample in zip(polys_arrays, samples):
            polys = polys[polys < max_up]
            _, bb, _ = plt.hist(
                polys,
                bins=params.bins if bb is None else bb,
                density=True,
                label=f"{sample} (µ={np.mean(polys):.2f}, σ={np.std(polys):.2f})",
                histtype="step",
            )
        plt.legend()
        plt.title("Poly(A,T) length distribution")
        plt.savefig(output.pdf, dpi=300)
        plt.savefig(output.png, dpi=1000)


rule tksm_abundance:
    input:
        binary=config["exec"]["tksm"],
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
    output:
        tsv=f"{plots_d}/expression_stats/{{sample}}.tksm.tsv",
    shell:
        "{input.binary} abundance"
        " -p {input.paf}"
        " -o {output.tsv}"


rule tksm_abundance_sc:
    input:
        binary=config["exec"]["tksm"],
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
        lr_matches=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_matches.tsv.gz",
    output:
        tsv=f"{plots_d}/expression_stats/{{sample}}.tksm_sc.tsv",
    shell:
        "{input.binary} abundance"
        " -p {input.paf}"
        " -m {input.lr_matches}"
        " -o {output.tsv}"


rule LIQA_refgene:
    input:
        gtf="{prefix}.gtf",
    output:
        refgen="{prefix}.gtf.refgene",
    shell:
        "python3 extern/LIQA/liqa_src/liqa.py"
        " -task refgene"
        " -format gtf"
        " -ref {input.gtf}"
        " -out {output.refgen}"
        " -m 1"


# Genion Rules
rule self_align_cdna:
    input:
        "{sample}",
    output:
        "{sample}.selfalign",
    threads: 32
    resources:
        time=60 * 24 - 1,
        mem_mb=64 * 1024,
    shell:
        "minimap2  -X -2 -c "
        " -t {threads}"
        " {input}"
        " {input}"
        " -o {output}.paf &&"
        "cat {output}.paf |"
        " cut -f1,6 |"
        " sed 's/_/\t/g' |"
        " awk 'BEGIN{{OFS=\"\\t\";}}{{print substr($1,1,15),substr($2,1,15),substr($3,1,15),substr($4,1,15);}}' |"
        " awk '$1!=$3' | sort | uniq"
        " > {output}"


rule genion_run:
    input:
        fastq=lambda wc: get_sample_fastqs(wc.sample),
        dna_paf=lambda wc: f"{preproc_d}/minimap2/{wc.sample}.DNA.paf",
        cdna_selfalign=lambda wc: get_sample_ref(wc.sample, "cDNA") + ".selfalign",
        gtf=lambda wc: get_sample_ref(wc.sample, "GTF"),
        dups=config["refs"]["genion_dups"],
    output:
        tsv=f"{preproc_d}/genion/{{sample}}.tsv",
        fail=f"{preproc_d}/genion/{{sample}}.tsv.fail",
    params:
        min_support=3,
    conda:
        "Benchmark_env.yaml"
    shell:
        "genion"
        " -i {input.fastq}"
        " -g {input.dna_paf}"
        " -o {output.tsv}"
        " --gtf {input.gtf}"
        " -s {input.cdna_selfalign}"
        " -d {input.dups}"
        " --min-support={params.min_support}"


rule longgf_run:
    input:
        bam=lambda wc: f"{preproc_d}/minimap2/{wc.sample}.DNA.ns.bam",
        gtf=lambda wc: get_sample_ref(wc.sample, "GTF"),
    output:
        tsv=f"{preproc_d}/longgf/{{sample}}.tsv",
    params:
        min_overlap=80,
        bin_size=4,
        min_map_len=80,
        min_support=3,
    conda:
        "Benchmark_env.yaml"
    shell:
        "LongGF"
        " {input.bam}"
        " {input.gtf}"
        " {params.min_overlap}"
        " {params.bin_size}"
        " {params.min_map_len}"
        " 0 0 {params.min_support}"
        " > {output.tsv}"


# Genion Rules End


rule minimap_dna_paf:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "DNA"),
    output:
        paf=f"{preproc_d}/minimap2/{{sample}}.DNA.paf",
    threads: 32
    shell:
        "minimap2"
        " -t {threads}"
        " -x splice"
        " -c"
        " {input.ref}"
        " {input.reads}"
        " > {output.paf}"


rule minimap_ns_dna:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "DNA"),
    output:
        bam=f"{preproc_d}/minimap2/{{sample}}.DNA.ns.bam",
    threads: 32
    shell:
        "minimap2"
        " -t {threads}"
        " -x splice"
        " -a"
        " {input.ref}"
        " {input.reads}"
        " | samtools view -hb "
        " -o {output.bam}"


rule minimap_dna:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "DNA"),
    output:
        bam=f"{preproc_d}/minimap2/{{sample}}.DNA.bam",
        bai=f"{preproc_d}/minimap2/{{sample}}.DNA.bam.bai",
    threads: 32
    shell:
        "minimap2"
        " -t {threads}"
        " -x splice"
        " -a"
        " {input.ref}"
        " {input.reads}"
        " | samtools sort "
        " -T {output.bam}.tmp"
        " -@ {threads}"
        " -o {output.bam}"
        " - "
        " && samtools index {output.bam}"


rule LIQA_quantify:
    input:
        refgen=lambda wc: f"{get_sample_ref(wc.sample, 'GTF')}.refgene",
        bam=f"{preproc_d}/minimap2/{{sample}}.DNA.bam",
    output:
        tsv=f"{plots_d}/expression_stats/{{sample}}.liqa.tsv",
    resources:
        time=60 * 6 - 1,
    threads: 16
    shell:
        "python3 extern/LIQA/liqa_src/liqa.py"
        " -task quantify"
        " -refgene <(cat {input.refgen})"
        " -bam {input.bam}"
        " -out {output.tsv}"
        " -max_distance 20"
        " -f_weight 1"
        " -threads {threads} > /dev/null"


rule minimap_quantify:
    input:
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
    output:
        tsv=f"{plots_d}/expression_stats/{{sample}}.minimap2.tsv",
    run:
        tid_to_count = Counter()
        for line in tqdm(
            open(input.paf), desc=f"[paf_quantify] Processing {input.paf}"
        ):
            if "tp:A:P" not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            tid = fields[5]
            tid_to_count[tid] += 1
        with open(output.tsv, "w") as f:
            f.write("transcript_id\tcount\n")
            for tid, count in tid_to_count.items():
                f.write(f"{tid}\t{count}\n")


rule NS_quantify_cp:
    input:
        tsv=f"{NS_d}/{{sample}}/abundance/sim_transcriptome_quantification.tsv",
    output:
        tsv=f"{plots_d}/expression_stats/{{sample}}.nanosim.tsv",
    shell:
        "cp {input.tsv} {output.tsv}"


rule tid_to_gid:
    input:
        gtf="{prefix}.gtf",
    output:
        tid_to_gid_pickle="{prefix}.gtf.tid_to_gid.pickle",
        gname_to_gid_pickle="{prefix}.gtf.gname_to_gid.pickle",
    run:
        tid_to_gid = dict()
        gname_to_gid = dict()
        for l in tqdm(open(input.gtf), desc=f"[tid_to_gid] Processing {input.gtf}"):
            if l[0] == "#":
                continue
            l = l.rstrip("\n").split("\t")
            if l[2] != "transcript":
                continue
            info = l[8]
            info = [x.strip().split(" ") for x in info.strip(";").split(";")]
            info = {x[0]: x[1].strip('"') for x in info}
            tid = info["transcript_id"]
            gid = info["gene_id"]
            gname = info.get("gene_name", gid)
            assert tid not in tid_to_gid
            tid_to_gid[tid] = gid

            gname_to_gid[gname] = gid
        pickle.dump(tid_to_gid, open(output.tid_to_gid_pickle, "wb+"))
        pickle.dump(gname_to_gid, open(output.gname_to_gid_pickle, "wb+"))


rule tpm_plot:
    input:
        tsvs=lambda wc: [
            f"{plots_d}/expression_stats/{s}.{wc.tpm_method}.tsv"
            for s in wc.samples.split(".")
        ],
        tid_to_gid_list=lambda wc: [
            f"{get_sample_ref(s, 'GTF')}.tid_to_gid.pickle"
            for s in wc.samples.split(".")
        ],
    output:
        png=f"{plots_d}/tpm_plot_{{tpm_method}}{{merge_type}}/{{samples}}.png",
    wildcard_constraints:
        tpm_method="|".join(tpm_method_settings.keys()),
        merge_type=".by_gene|",
    run:
        by_gene_title = ""
        if wildcards.merge_type == ".by_gene":
            by_gene_title = "\nmerged by gene"
            tid_to_gid = dict()
            for tid_to_gid_file in input.tid_to_gid_list:
                tid_to_gid.update(pickle.load(open(tid_to_gid_file, "rb")))
        else:
            tid_to_gid = None
        my_get_tpm_args = {
            **tpm_method_settings[wildcards.tpm_method]["get_tpm_args"],
            **dict(tid_to_gid=tid_to_gid),
        }
        my_get_tpm = functools.partial(
            get_tpm,
            **my_get_tpm_args,
        )
        X_tpm = my_get_tpm(input.tsvs[0])
        Y_tpms = list(map(my_get_tpm, input.tsvs[1:]))
        samples = wildcards.samples.split(".")
        plot_tpm_func(
            X_tpm,
            Y_tpms,
            samples,
            output,
            title=tpm_method_settings[wildcards.tpm_method]["title"] + by_gene_title,
        )


def get_tpm(
    tsv,
    min_tpm=0.01,
    min_val=0.0,
    key_col=0,
    val_col=1,
    header=True,
    sep="\t",
    tid_to_gid=None,
):
    tpm = Counter()
    if isinstance(key_col, int):
        key_col = [key_col]
    for line_num, line in tqdm(
        enumerate(open(tsv, "r")), desc=f"[get_tpm] Processing {tsv}"
    ):
        if header and line_num == 0:
            continue
        line = line.rstrip().split(sep)
        if tid_to_gid is not None:
            tid = line[key_col[0]]
            tid = tid.split(".")[0]
            gid = tid_to_gid.get(tid, tid)
            key = tuple([gid] + [line[kc] for kc in key_col[1:]])
        else:
            key = tuple([line[kc] for kc in key_col])
        val = float(line[val_col])
        if val > min_val:
            tpm[key] = val
    thruput = sum(tpm.values())
    for k, v in tpm.items():
        tpm[k] = v / thruput * 1e6
    tpm = Counter({k: v for k, v in tpm.items() if v > min_tpm})
    thruput = sum(tpm.values())
    for k, v in tpm.items():
        tpm[k] = v / thruput * 1e6
    return tpm


def plot_tpm_func(X_tpm, Y_tpms, samples, outpaths, title, flip=True):
    plt.rc("font", size=22)

    plot_count = len(samples) - 1
    unions = [
        (lambda xy: xy[0], "Transcripts in input (regardless of sample)"),
        # (lambda xy: xy[0] & xy[1], "Transcripts in sample AND in input"),
        # (lambda xy: xy[0] | xy[1], "Transcripts in sample OR in input"),
    ]
    if flip == False:
        fig, axs = plt.subplots(
            plot_count,
            len(unions),
            sharex=True,
            sharey=True,
            figsize=(
                10 * len(unions),
                10 * plot_count,
            ),
            squeeze=False,
        )
    else:
        fig, axs = plt.subplots(
            len(unions),
            plot_count,
            sharex=True,
            sharey=True,
            figsize=(
                10 * plot_count,
                10 * len(unions),
            ),
            squeeze=False,
        )
    fig.suptitle(title, fontsize=30)

    if flip == False:
        for (_, title), ax in zip(unions, axs[0]):
            ax.set_title(title, fontsize=28)
        for sample, ax in zip(samples[1:], axs[:, 0]):
            ax.set_ylabel(f"TPM in {sample}", fontsize=28)
        for ax in axs[-1]:
            ax.set_xlabel(f"TPM in input ({samples[0]})", fontsize=28)
    else:
        pass
        for c_idx, sample in enumerate(samples[1:]):
            for ax in axs[:, c_idx]:
                ax.set_ylabel(f"TPM in {sample}", fontsize=28)
                if c_idx == 0 and len(unions) > 1:
                    title = unions[c_idx][1]
                    ax.set_ylabel(f"{title}\n\nTPM in {sample}", fontsize=28)

        for ax in axs[-1, :]:
            ax.set_xlabel(f"TPM in input ({samples[0]})", fontsize=28)

    X_tids = set(X_tpm.keys())
    for idx, (sample, Y_tpm) in enumerate(zip(samples[1:], Y_tpms)):
        Y_tids = set(Y_tpm.keys())
        for ax, (select_tids_func, _) in zip(
            axs[idx] if not flip else axs[:, idx],
            unions,
        ):
            select_tids = select_tids_func((X_tids, Y_tids))
            X = np.array([X_tpm[tid] for tid in select_tids])
            Y = np.array([Y_tpm[tid] for tid in select_tids])
            ax.plot(X, Y, "o", alpha=0.5, label=sample)
            t = ax.text(
                0.15,
                0.75,
                f"N = {len(select_tids)}\n"
                + f"Sample: {sample}\n"
                + f"r-squared = {r2_score(X,Y):.3f}\n"
                + f"RMSE = {mean_squared_error(X,Y, squared=False):.1f}",
                transform=ax.transAxes,
            )
            t.set_bbox(dict(facecolor="red", alpha=0.5, edgecolor="red"))
            ax.set_xscale("log")
            ax.set_yscale("log")
    fig.tight_layout()
    for outpath in outpaths:
        plt.savefig(outpath)


rule lr_cell_stats:
    input:
        script=config["exec"]["lr_cell_stats"],
        barcodes=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.bc_whitelist.tsv.gz",
        reads=lambda wc: get_sample_fastqs(wc.sample),
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
        lr_matches=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_matches.tsv.gz",
        lr_br=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz",
    output:
        barcodes_p=f"{plots_d}/lr_cell_stats/{{sample}}.lr_cell_stats.barcodes.pickle",
        reads_p=f"{plots_d}/lr_cell_stats/{{sample}}.lr_cell_stats.reads.pickle",
        seqs_p=f"{plots_d}/lr_cell_stats/{{sample}}.lr_cell_stats.seqs.pickle",
    params:
        outpath=f"{plots_d}/lr_cell_stats/{{sample}}.lr_cell_stats",
    shell:
        "python3 {input.script}"
        " -b {input.barcodes}"
        " -r {input.reads}"
        " -p {input.paf}"
        " -m {input.lr_matches}"
        " -l {input.lr_br}"
        " -o {params.outpath}"


rule lr_cell_matching:
    input:
        script=config["exec"]["lr_cell_matching"],
        barcodes=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.bc_whitelist.tsv.gz",
        reads=lambda wc: get_sample_fastqs(wc.sample),
    output:
        tsv=f"{plots_d}/lr_cell_stats/{{sample}}.lr_cell_matching.tsv",
    threads: 64
    shell:
        "python3 {input.script}"
        " -i {input.reads}"
        " -b {input.barcodes}"
        " -o {output.tsv}"
        " -t {threads}"


rule lr_cell_matching_plot:
    input:
        tsvs=lambda wc: [
            f"{plots_d}/lr_cell_stats/{s}.lr_cell_matching.tsv"
            for s in wc.samples.split(".")
        ],
    output:
        pdf=f"{plots_d}/lr_cell_matching_plot/{{samples}}.pdf",
        png=f"{plots_d}/lr_cell_matching_plot/{{samples}}.png",
    run:
        samples = wildcards.samples.split(".")
        fig = plt.figure(figsize=(8, 4))
        ax = fig.add_subplot(111)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
        width = 0.4
        counts = dict()
        for IDX, (tsv, sample) in enumerate(zip(input.tsvs, samples)):
            for idx, l in enumerate(open(tsv)):
                if idx == 0:
                    continue
                l = l.rstrip("\n").split("\t")
                k = tuple(l[:2])
                v = float(l[3].strip("%"))
                if not k in counts:
                    counts[k] = [0.0 for _ in samples]
                counts[k][IDX] = v
        x_ticks, Y = zip(*sorted(counts.items()))
        x_ticks = [f"{x[0]}\n d={x[1]}" for x in x_ticks]
        X = np.array(np.arange(len(x_ticks)), dtype=float)
        ax.set_xticks(X + (len(samples) - 1) * width / 2, x_ticks)

        for sample, sample_Y in zip(samples, zip(*Y)):
            ax.bar(X, sample_Y, label=sample, width=width)
            for x, y in zip(X, sample_Y):
                plt.text(
                    x,
                    y + 0.5 if y < 6.5 else y - 6.5,
                    f"{y:.1f}%",
                    ha="center",
                    size="small",
                    rotation=-60,
                )
            X += width
        plt.legend()
        fig.tight_layout()
        plt.savefig(output.pdf, dpi=300)
        plt.savefig(output.png, dpi=1000)


rule lr_cell_plot:
    input:
        reads_pickles=lambda wc: [
            f"{plots_d}/lr_cell_stats/{s}.lr_cell_stats.reads.pickle"
            for s in wc.samples.split(".")
        ],
    output:
        pdf=f"{plots_d}/lr_sr_adapt/{{samples}}.pdf",
        png=f"{plots_d}/lr_sr_adapt/{{samples}}.png",
    run:
        samples = wildcards.samples.split(".")
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
        fig.suptitle(
            f"Density distribution of long-reads, with a single primary alignment,\n in terms of strand and Illumina adapter mapping location"
        )
        sample_reads = list()
        for p in tqdm(
            input.reads_pickles,
            total=len(samples),
            desc="[lr_cell_plot] Loading pickles",
        ):
            sample_reads.append(pickle.load(open(p, "rb")))
        totals = [
            sum([1 for read in reads if len(read["mappings"]) == 1])
            for reads in sample_reads
        ]
        for ax, Tstrand in zip(axes, ["+", "-"]):
            ax.set_title(f"Transcriptome strand {Tstrand}")
            for reads, sample, total in tqdm(
                zip(sample_reads, samples, totals),
                total=len(samples),
                desc=f"[lr_cell_plot] Processing samples on transcriptome strand {Tstrand}",
            ):
                vals = Counter()
                for read in tqdm(reads, total=len(reads), desc=f"Processing {sample}"):
                    if len(read["mappings"]) != 1:
                        continue
                    tid, tstrand = list(read["mappings"])[0][:2]
                    if tstrand != Tstrand:
                        continue
                    vals[read["br_seg"][1]] += 1 / total
                bin_edges = np.arange(-150, 151, 5)
                keys, weights = zip(*vals.items())
                hist, _ = np.histogram(keys, bins=bin_edges, weights=weights)
                print(hist)
                MIN_VAL = 0.01
                f = 0
                for i, v in enumerate(hist):
                    if v > MIN_VAL:
                        f = i
                        break
                l = len(hist)
                for i, v in enumerate(reversed(hist)):
                    if v > MIN_VAL:
                        l = len(hist) - i
                        break
                ax.bar(
                    bin_edges[f:l],
                    hist[f:l],
                    width=5,
                    align="edge",
                    alpha=0.3,
                    label=sample,
                )
            ax.legend()
        fig.tight_layout()
        plt.savefig(output.pdf, dpi=300)
        plt.savefig(output.png, dpi=1000)


### Gene fusion plots
def get_last_mdf(exprmnt, rule_name):
    step_names = [
        list(step)[0] for step in config["TS_experiments"][exprmnt]["pipeline"]
    ]
    try:
        rule_idx = len(step_names) - 1 - step_names[::-1].index(rule_name)
        prefix = ".".join(step_names[: rule_idx + 1])
        return f"{TS_d}/{exprmnt}/{prefix}.mdf"
    except ValueError:
        if step_names[0] == "Mrg":
            sources = config["TS_experiments"][exprmnt]["pipeline"][0]["Mrg"]["sources"]
            return get_last_mdf(sources[0], rule_name)
        else:
            raise Exception(f"Rule {rule_name} not found in {exprmnt} pipeline")


rule gene_fusion_intersection:
    input:
        script=config["exec"]["gene_fusion_intersection"],
        tid_to_gid=lambda wc: f"{get_sample_ref(wc.sample, 'GTF')}.tid_to_gid.pickle",
        gname_to_gid=lambda wc: f"{get_sample_ref(wc.sample, 'GTF')}.gname_to_gid.pickle",
        tsb=lambda wc: get_last_mdf(wc.sample, "Tsb"),
        mrg=lambda wc: get_last_mdf(wc.sample, "Mrg"),
        uns=lambda wc: get_last_mdf(wc.sample, "Uns"),
        genion_pass=f"{preproc_d}/genion/{{sample}}.tsv",
        genion_fail=f"{preproc_d}/genion/{{sample}}.tsv.fail",
        longGF=f"{preproc_d}/longgf/{{sample}}.tsv",
    output:
        pickle=f"{plots_d}/gene_fusion_intersection/{{sample}}.pickle",
    conda:
        "Benchmark_env.yaml"
    shell:
        "python {input.script}"
        " --tid_to_gid {input.tid_to_gid}"
        " --gname_to_gid {input.gname_to_gid}"
        " --tsb {input.tsb}"
        " --mrg {input.mrg}"
        " --uns {input.uns}"
        " --genion_pass {input.genion_pass}"
        " --genion_fail {input.genion_fail}"
        " --longGF {input.longGF}"
        " --output {output.pickle}"


rule gene_fusion_upset_plot:
    input:
        script=config["exec"]["gene_fusion_upset_plot"],
        pickle=f"{plots_d}/gene_fusion_intersection/{{sample}}.pickle",
    output:
        pdf=f"{plots_d}/gene_fusion_upset/{{sample}}.pdf",
        png=f"{plots_d}/gene_fusion_upset/{{sample}}.png",
    params:
        outpath=f"{plots_d}/gene_fusion_upset/{{sample}}",
    conda:
        "Benchmark_env.yaml"
    shell:
        "python {input.script}"
        " --pickle {input.pickle}"
        " --sample {wildcards.sample}"
        " --outpath {params.outpath}"


rule gene_fusion_accuracy:
    input:
        pickles=lambda wc: [
            f"{plots_d}/gene_fusion_intersection/{s}.pickle"
            for s in wc.samples.split(".")
        ],
    output:
        csv=f"{plots_d}/gene_fusion_accuracy/{{samples}}.csv",
    run:
        tools = ["LongGF", "PASS:GF", "Truth"]
        delim = ","
        outfile = open(output.csv, "w+")
        record = list()
        record.append("Truncation")
        record.append("Glue rate")
        record.append("Tool")
        record.append("#")
        record.append("TP")
        record.append("FP")
        record.append("FN")
        record.append("F1")
        print(delim.join(record), file=outfile)
        for p in input.pickles:
            gene_fusions = pickle.load(open(p, "rb"))
            p = p.split("/")[-1][len("TKSM_gene_fusion") :][: -len(".pickle")]
            trc = False
            if p.endswith("_trc"):
                trc = True
                p = p[: -len("_trc")]
            rate = float(p) / 100
            tool_gf = {t: set() for t in tools}
            for k, v in gene_fusions.items():
                for t in v:
                    if not t in tool_gf:
                        continue
                    tool_gf[t].add(k)
            X = tool_gf["Truth"]
            for t in tools:
                Y = tool_gf[t]
                TP = len(Y & X)
                FP = len(Y - X)
                FN = len(X - Y)
                F1 = TP / (TP + 0.5 * (FP + FN))

                if t == "PASS:GF":
                    t = "Genion"
                record = [
                    f"{trc}",
                    f"{rate:.0%}",
                    t,
                    f"{len(Y):d}",
                    f"{TP:d}",
                    f"{FP:d}",
                    f"{FN:d}",
                    f"{F1:.2f}",
                ]
                print(delim.join(record), file=outfile)
