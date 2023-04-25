import sys
import gzip
from multiprocessing import Pool
import pickle

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
time_d = f"{outpath}/time"
plots_d = f"{outpath}/plots"


use rule * from TS_smk as TS_*


use rule minimap_cdna from TS_smk as TS_minimap_cdna with:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "cDNA"),


use rule scTagger_lr_seg from TS_smk as TS_scTagger_lr_seg with:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),


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
        title="TPM using TKSM single cell abundance estimates",
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


def get_sample_ref_names(sample):
    try:
        # Try should work for TS and real samples
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


rule all:
    input:
        [
            f"{plots_d}/{r}/{'.'.join(p['data'])}.svg"
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
        svg=f"{plots_d}/raw_lengths/{{samples}}.svg",
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
                alpha=0.3,
            )
            print(f"{sample} length mean: {np.mean(lens):.2f}")
            print(f"{sample} length std: {np.std(lens):.2f}")
            print(f"{sample} length median: {np.median(lens):.2f}")
        plt.legend()
        plt.title("Length distribution (whole reads)")
        plt.savefig(output[0], dpi=300)


rule mapped_raw_lengths:
    input:
        pafs=lambda wc: [
            f"{preproc_d}/minimap2/{s}.cDNA.paf" for s in wc.samples.split(".")
        ],
    output:
        svg=f"{plots_d}/mapped_raw_lengths/{{samples}}.svg",
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
                alpha=0.3,
            )
            print(f"{sample} length (PAF) mean: {np.mean(lens):.2f}")
            print(f"{sample} length (PAF) std: {np.std(lens):.2f}")
            print(f"{sample} length (PAF) median: {np.median(lens):.2f}")
        plt.xlim(left=0, right=max_up)
        plt.legend()
        plt.title("Length distribution (whole reads from PAF)")
        plt.savefig(output[0], dpi=300)


rule mapped_lengths:
    input:
        pafs=lambda wc: [
            f"{preproc_d}/minimap2/{s}.cDNA.paf" for s in wc.samples.split(".")
        ],
    output:
        svg=f"{plots_d}/mapped_lengths/{{samples}}.svg",
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
                start = int(fields[2])
                end = int(fields[3])
                lens.append(end - start)
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
                alpha=0.3,
            )
            print(f"{sample} mapping length (PAF) mean: {np.mean(lens):.2f}")
            print(f"{sample} mapping length (PAF) std: {np.std(lens):.2f}")
            print(f"{sample} mapping length (PAF) median: {np.median(lens):.2f}")
        plt.xlim(left=0, right=max_up)
        plt.legend()
        plt.title("Length distribution (mapping part of reads from PAF)")
        plt.savefig(output[0], dpi=300)


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
        f"{plots_d}/substitution/{{samples}}.svg",
    params:
        bins=50,
    run:
        samples = wildcards.samples.split(".")
        fig, axs = plt.subplots(2, 2, figsize=(12, 4))
        axs = axs.flatten()
        bb = [None, None, None, None]
        maxi = [0, 0, 0, 0]
        print("Reading TSV files...")
        for tsv, sample in tqdm(
            zip(input.tsvs, samples),
            total=len(input.tsvs),
            desc="[substitution] Reading TSV files",
        ):
            counts = list()
            for line in tqdm(open(tsv, "r"), desc=f"[substitution] Processing {tsv}"):
                match_c, subs_c, insert_c, del_c, length = [
                    int(float(x)) for x in line.rstrip().split("\t")
                ]
                l100 = 100 / length
                counts.append(
                    (
                        l100 * match_c,
                        l100 * subs_c,
                        l100 * insert_c,
                        l100 * del_c,
                    )
                )
            for j, title in enumerate(
                [
                    "match counts",
                    "substitution counts",
                    "insertion counts",
                    "deletion counts",
                ]
            ):
                cc = np.array([x[j] for x in counts])
                _, bb[j], _ = axs[j].hist(
                    cc,
                    bins=params.bins if bb[j] is None else bb[j],
                    alpha=0.3,
                    label=f"{sample}: (µ={np.mean(cc):.2f}, σ={np.std(cc):.2f})",
                    density=True,
                )
                axs[j].legend()
        for i, title in enumerate(
            [
                "match counts",
                "substitution counts",
                "insertion counts",
                "deletion counts",
            ]
        ):
            axs[i].set_title(f"{title} per 100 bases")
        fig.tight_layout()
        plt.savefig(output[0], dpi=300)


rule polyA:
    input:
        fastqs=lambda wc: itertools.chain.from_iterable(
            [get_sample_fastqs(s) for s in wc.samples.split(".")]
        ),
    output:
        svg=f"{plots_d}/polyA/{{samples}}.svg",
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
                alpha=0.5,
            )
        plt.legend()
        plt.title("Poly(A,T) length distribution")
        plt.savefig(output[0], dpi=300)


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


rule minimap_dna_paf:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "cDNA"),
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


rule minimap_dna:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "DNA"),
    output:
        bam=f"{preproc_d}/minimap2/{{sample}}.DNA.bam",
        bai=f"{preproc_d}/minimap2/{{sample}}.DNA.bam.bai",
    benchmark:
        f"{time_d}/{{sample}}/minimap2_cdna.benchmark"
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
        pickle="{prefix}.gtf.tid_to_gid.pickle",
    run:
        tid_to_gid = dict()
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
            assert tid not in tid_to_gid
            tid_to_gid[tid] = gid
        pickle.dump(tid_to_gid, open(output.pickle, "wb+"))


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
        svg=f"{plots_d}/tpm_plot_{{tpm_method}}{{merge_type}}/{{samples}}.svg",
    wildcard_constraints:
        tpm_method="|".join(tpm_method_settings.keys()),
        merge_type=".by_gene|",
    run:
        if wildcards.merge_type == ".by_gene":
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
            output.svg,
            title=tpm_method_settings[wildcards.tpm_method]["title"],
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


def plot_tpm_func(X_tpm, Y_tpms, samples, outpath, title):
    plt.rc("font", size=22)

    plot_count = len(samples) - 1
    fig, axs = plt.subplots(
        plot_count,
        3,
        sharex=True,
        sharey=True,
        figsize=(10 * 3, 10 * plot_count),
        squeeze=False,
    )
    fig.suptitle(title, fontsize=30)
    axs[0, 0].set_title("Transcripts in sample AND in input", fontsize=28)
    axs[0, 1].set_title("Transcripts in sample OR in input", fontsize=28)
    axs[0, 2].set_title("Transcripts in input (regardless of sample)", fontsize=28)

    for sample, ax in zip(samples[1:], axs[:, 0]):
        ax.set_ylabel(f"TPM in {sample}", fontsize=28)
    for ax in axs[-1]:
        ax.set_xlabel(f"TPM in input ({samples[0]})", fontsize=28)
    X_tids = set(X_tpm.keys())
    for idx, (sample, Y_tpm) in enumerate(zip(samples[1:], Y_tpms)):
        Y_tids = set(Y_tpm.keys())

        ax = axs[idx, 0]
        for ax, select_tids in zip(
            axs[idx],
            [
                X_tids & Y_tids,
                X_tids | Y_tids,
                X_tids,
            ],
        ):
            X = np.array([X_tpm[tid] for tid in select_tids])
            Y = np.array([Y_tpm[tid] for tid in select_tids])
            ax.plot(X, Y, "o", alpha=0.5, label=sample)
            ax.text(
                0.15,
                0.75,
                f"N = {len(select_tids)}\n"
                + f"Sample: {sample}\n"
                + f"r-squared = {r2_score(X,Y):.3f}\n"
                + f"RMSE = {mean_squared_error(X,Y, squared=False):.1f}",
                transform=ax.transAxes,
            )
            ax.set_xscale("log")
            ax.set_yscale("log")
    fig.tight_layout()
    plt.savefig(outpath)


rule NS_analysis:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        dna=lambda wc: get_sample_ref(wc.sample, "DNA"),
        cdna=lambda wc: get_sample_ref(wc.sample, "cDNA"),
        gtf=lambda wc: get_sample_ref(wc.sample, "GTF"),
    output:
        analysis_dir=directory(f"{NS_d}/{{sample}}/analysis"),
    benchmark:
        f"{time_d}/{{sample}}/NS_analysis.benchmark"
    params:
        output_prefix=f"{NS_d}/{{sample}}/analysis/sim",
    threads: 32
    shell:
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
    output:
        quantify_tsv=f"{NS_d}/{{sample}}/abundance/sim_transcriptome_quantification.tsv",
    benchmark:
        f"{time_d}/{{sample}}/NS_quantify.benchmark"
    params:
        output_prefix=f"{NS_d}/{{sample}}/abundance/sim",
    threads: 32
    resources:
        time=60 * 6 - 1,
    shell:
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
    output:
        fasta=f"{NS_d}/{{exprmnt}}/simulation_aligned_reads.fasta",
    benchmark:
        f"{time_d}/{{exprmnt}}/NS_simulate.benchmark"
    params:
        analysis_model=lambda wc: f"{NS_d}/{NS_exprmnt_sample(wc.exprmnt)}/analysis/sim",
        out_prefix=lambda wc: f"{NS_d}/{wc.exprmnt}/simulation",
        other=lambda wc: config["NS_experiments"][wc.exprmnt]["simulation_params"],
    threads: 32
    shell:
        "simulator.py transcriptome"
        " -rt {input.cdna}"
        " -rg {input.dna}"
        " --exp {input.quantify_tsv}"
        " -c {params.analysis_model}"
        " -o {params.out_prefix}"
        " -t {threads}"
        " --no_model_ir"
        " {params.other}"


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


rule lr_cell_plot:
    input:
        reads_pickles=lambda wc: [
            f"{plots_d}/lr_cell_stats/{s}.lr_cell_stats.reads.pickle"
            for s in wc.samples.split(".")
        ],
    output:
        f"{plots_d}/lr_sr_adapt/{{samples}}.svg",
    run:
        samples = wildcards.samples.split(".")
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        fig.suptitle(f"Read with one unique primary mapping")
        sample_reads = list()
        for p in tqdm(
            input.reads_pickles,
            total=len(samples),
            desc="[lr_cell_plot] Loading pickles",
        ):
            sample_reads.append(pickle.load(open(p, "rb")))
        for ax, Tstrand in zip(axes, ["+", "-"]):
            ax.set_title(f"Transcriptome strand {Tstrand}")
            for reads, sample in tqdm(
                zip(sample_reads, samples),
                total=len(samples),
                desc=f"[lr_cell_plot] Processing samples on transcriptome strand {Tstrand}",
            ):
                vals = list()
                for read in tqdm(reads, total=len(reads), desc=f"Processing {sample}"):
                    if len(read["mappings"]) != 1:
                        continue
                    tid, tstrand = list(read["mappings"])[0][:2]
                    if tstrand != Tstrand:
                        continue
                    vals.append(read["br_seg"][1])
                ax.hist(
                    vals,
                    np.arange(-100, 100, 5),
                    density=True,
                    alpha=0.3,
                    label=sample,
                )
            ax.legend()
        plt.legend()
        plt.savefig(output[0], dpi=300)
