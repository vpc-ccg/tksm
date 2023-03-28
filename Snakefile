import sys
import re

if len(config) == 0:

    configfile: "config.yaml"


outpath = config["outpath"]
preproc_d = f"{outpath}/preprocess"
TS_d = f"{outpath}/TS"
time_d = f"{outpath}/time"
DEBUG = True

if DEBUG:

    def pipe(X):
        return X

    def temp(X):
        return X


def exprmnt_sample(exprmnt):
    return config["TS_experiments"][exprmnt]["sample"]


def get_split_mdf(wc, SPL_K):
    pipeline = config["TS_experiments"][wc.exprmnt]["pipeline"]
    prefix = f"{wc.prefix}.Spl"
    step = get_step_pipeline(pipeline, [x.split(".") for x in prefix.split("/")])
    step = step["Spl"][SPL_K]
    suffix = ".".join([list(x.keys())[0] for x in step])
    if step == list():
        return list()
    return f"{TS_d}/{wc.exprmnt}/{prefix}/{SPL_K}.{suffix}.mdf"


def get_step_pipeline(pipeline, prefix):
    if isinstance(pipeline, dict):
        pipeline = pipeline[prefix[0][0]]
        prefix[0] = prefix[0][1:]
    pipeline_idx = len(prefix[0]) - 1
    assert len(prefix) > 0
    assert pipeline_idx >= 0

    if len(prefix) == 1:
        return pipeline[pipeline_idx]
    else:
        return get_step_pipeline(pipeline[pipeline_idx]["Spl"], prefix[1:])


def module_params(wildcards, rule_name):
    pipeline = config["TS_experiments"][wildcards.exprmnt]["pipeline"]
    prefix = [x.split(".") for x in f"{wildcards.prefix}.{rule_name}".split("/")]
    step = get_step_pipeline(pipeline, prefix)
    return step[rule_name]["params"]


def experiment_prefix(exprmnt):
    prefix = list()
    for component in config["TS_experiments"][exprmnt]["pipeline"]:
        prefix.append(list(component)[0])
    return ".".join(prefix)


def get_sample_ref_name(sample):
    if sample in config["samples"]:
        return config["samples"][sample]["ref"]
    else:
        return get_sample_ref_name(config["TS_experiments"][sample]["sample"])


def get_sample_ref(sample, ref_type):
    ref_name = get_sample_ref_name(sample)
    return config["refs"][ref_name][ref_type]


def get_sample_fastqs(name):
    if name in config["samples"]:
        sample = name
        return config["samples"][sample]["fastq"]
    elif name in config["TS_experiments"]:
        exprmnt = name
        return [f"{TS_d}/{exprmnt}/{experiment_prefix(exprmnt)}.fastq"]
    else:
        raise ValueError(f"Invalid experiment/sample name! {name}")


rule all:
    input:
        [
            f"{TS_d}/{exprmnt}/{experiment_prefix(exprmnt)}.fastq"
            for exprmnt in config["TS_experiments"]
        ],


rule make_binary:
    input:
        "src/{prefix}.cpp",
    output:
        "build/{prefix}",
    shell:
        "make {output}"


rule sequencer:
    input:
        obj=["build/obj/sequencer.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        fastas=lambda wc: get_sample_ref(wc.exprmnt, "DNA"),
        qscore_model=lambda wc: f"{preproc_d}/models/badread/{exprmnt_sample(wc.exprmnt)}.qscore.gz",
        error_model=lambda wc: f"{preproc_d}/models/badread/{exprmnt_sample(wc.exprmnt)}.error.gz",
    output:
        fastq=f"{TS_d}/{{exprmnt}}/{{prefix}}.Seq.fastq",
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.Seq.benchmark"
    threads: 32
    params:
        other=lambda wc: module_params(wc, "Seq"),
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt="|".join([re.escape(x) for x in config["TS_experiments"]]),
    shell:
        "{params.binary} sequencer"
        " -i {input.mdf}"
        " --references {input.fastas}"
        " -o {output.fastq}"
        " --threads {threads}"
        " --badread-error-model={input.error_model}"
        " --badread-qscore-model={input.qscore_model}"
        " {params.other}"


rule split_merge:
    input:
        mdf_f=lambda wc: get_split_mdf(wc, "Spl_T"),
        mdf_t=lambda wc: get_split_mdf(wc, "Spl_F"),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Spl.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.Spl.benchmark"
    wildcard_constraints:
        exprmnt="|".join([re.escape(x) for x in config["TS_experiments"]]),
    shell:
        "cat {input.mdf_f} {input.mdf_t} > {output.mdf}"


rule split:
    input:
        obj=["build/obj/filter.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf_f=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Spl/Spl_T.mdf"),
        mdf_t=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Spl/Spl_F.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.Spl.benchmark"
    params:
        other=lambda wc: module_params(wc, "Spl"),
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt="|".join([re.escape(x) for x in config["TS_experiments"]]),
    shell:
        "{params.binary} filter"
        " -i {input.mdf}"
        " -t {output.mdf_t}"
        " -f {output.mdf_f}"
        " {params.other}"


rule truncate:
    input:
        obj=["build/obj/truncate.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        x=lambda wc: f"{preproc_d}/truncate_kde/{exprmnt_sample(wc.exprmnt)}.X_idxs.npy",
        y=lambda wc: f"{preproc_d}/truncate_kde/{exprmnt_sample(wc.exprmnt)}.Y_idxs.npy",
        g=lambda wc: f"{preproc_d}/truncate_kde/{exprmnt_sample(wc.exprmnt)}.grid.npy",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Trc.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.Trc.benchmark"
    params:
        other=lambda wc: module_params(wc, "Trc"),
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt="|".join([re.escape(x) for x in config["TS_experiments"]]),
    shell:
        "{params.binary} truncate"
        " -i {input.mdf}"
        " --kde={input.g},{input.x},{input.y}"
        " -o {output.mdf}"
        " {params.other}"


rule flip:
    input:
        obj=["build/obj/strand_man.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Flp.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.Flp.benchmark"
    params:
        other=lambda wc: module_params(wc, "Flp"),
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt="|".join([re.escape(x) for x in config["TS_experiments"]]),
    shell:
        "{params.binary} flip"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.other}"


rule pcr:
    input:
        obj=["build/obj/pcr.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.PCR.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.PCR.benchmark"
    params:
        other=lambda wc: module_params(wc, "PCR"),
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt="|".join([re.escape(x) for x in config["TS_experiments"]]),
    shell:
        "{params.binary} pcr"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.other}"


rule tag:
    input:
        obj=["build/obj/tag.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Tag.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.Tag.benchmark"
    params:
        other=lambda wc: module_params(wc, "Tag"),
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt="|".join([re.escape(x) for x in config["TS_experiments"]]),
    shell:
        "{params.binary} tag"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.other}"


rule single_cell_barcoder:
    input:
        obj=["build/obj/single-cell-barcoder.o", "build/obj/tksm.o"]
        if DEBUG
        else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.SCB.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.SCB.benchmark"
    params:
        other=lambda wc: module_params(wc, "SCB"),
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt="|".join([re.escape(x) for x in config["TS_experiments"]]),
    shell:
        "{params.binary} single-cell-barcoder"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.other}"


rule polyA:
    input:
        obj=["build/obj/polyA.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.plA.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.plA.benchmark"
    params:
        other=lambda wc: module_params(wc, "plA"),
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt="|".join([re.escape(x) for x in config["TS_experiments"]]),
    shell:
        "{params.binary} polyA"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.other}"


rule splicer:
    input:
        obj=["build/obj/splicer.o", "build/obj/tksm.o"] if DEBUG else list(),
        tsv=f"{TS_d}/{{exprmnt}}/{{prefix}}.tsv",
        gtf=lambda wc: get_sample_ref(wc.exprmnt, "GTF"),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Spc.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.Spc.benchmark"
    params:
        other=lambda wc: module_params(wc, "Spc"),
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt="|".join([re.escape(x) for x in config["TS_experiments"]]),
    shell:
        "{params.binary} splicer"
        " -a {input.tsv}"
        " -g {input.gtf}"
        " -o {output.mdf}"
        " {params.other}"


rule abundance:
    input:
        obj=["build/obj/abundance.o", "build/obj/tksm.o"] if DEBUG else list(),
        paf=lambda wc: f"{preproc_d}/minimap2/{exprmnt_sample(wc.exprmnt)}.cDNA.paf",
    output:
        tsv=f"{TS_d}/{{exprmnt}}/Xpr.tsv",
    benchmark:
        f"{time_d}/{{exprmnt}}/Xpr.benchmark"
    params:
        binary=config["exec"]["tksm"],
    shell:
        "{params.binary} abundance"
        " -p {input.paf}"
        " -o {output.tsv}"


rule abundance_sc:
    input:
        obj=["build/obj/abundance.o", "build/obj/tksm.o"] if DEBUG else list(),
        paf=lambda wc: f"{preproc_d}/minimap2/{exprmnt_sample(wc.exprmnt)}.cDNA.paf",
        lr_matches=lambda wc: f"{preproc_d}/scTagger/{exprmnt_sample(wc.exprmnt)}/{exprmnt_sample(wc.exprmnt)}.lr_matches.tsv.gz",
    output:
        tsv=f"{TS_d}/{{exprmnt}}/Xpr_sc.tsv",
    benchmark:
        f"{time_d}/{{exprmnt}}/Xpr_sc.benchmark"
    params:
        binary=config["exec"]["tksm"],
    shell:
        "{params.binary} abundance"
        " -p {input.paf}"
        " -m {input.lr_matches}"
        " -o {output.tsv}"


rule truncate_kde:
    input:
        obj=["build/obj/kde.o", "build/obj/tksm.o"] if DEBUG else list(),
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
    output:
        x=f"{preproc_d}/truncate_kde/{{sample}}.X_idxs.npy",
        y=f"{preproc_d}/truncate_kde/{{sample}}.Y_idxs.npy",
        g=f"{preproc_d}/truncate_kde/{{sample}}.grid.npy",
    benchmark:
        f"{time_d}/{{sample}}/truncate_kde.benchmark"
    params:
        out_prefix=f"{preproc_d}/truncate_kde/{{sample}}",
        binary=config["exec"]["tksm"],
    threads: 32
    shell:
        "{params.binary} kde"
        " -i {input.paf}"
        " -o {params.out_prefix}"
        " --threads {threads}"


rule minimap_cdna:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "cDNA"),        
    output:
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
    benchmark:
        f"{time_d}/{{sample}}/minimap2_cdna.benchmark"
    threads: 32
    shell:
        "minimap2"
        " -t {threads}"
        " -x map-ont"
        " -c --eqx"
        " -o {output.paf}"
        " {input.ref}"
        " {input.reads}"


rule scTagger_match:
    input:
        lr_tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz",
        wl_tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.bc_whitelist.tsv.gz",
    output:
        lr_tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_matches.tsv.gz",
    benchmark:
        f"{time_d}/{{sample}}/scTagger_match.benchmark"
    threads: 32
    shell:
        "scTagger.py match_trie"
        " -lr {input.lr_tsv}"
        " -sr {input.wl_tsv}"
        " -o {output.lr_tsv}"
        " -t {threads}"


rule scTagger_extract_bc:
    input:
        tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz",
        wl=config["refs"]["10x_bc"],
    output:
        tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.bc_whitelist.tsv.gz",
    benchmark:
        f"{time_d}/{{sample}}/scTagger_extract_bc.benchmark"
    shell:
        "scTagger.py extract_sr_bc_from_lr"
        " -i {input.tsv}"
        " -wl {input.wl}"
        " -o {output.tsv}"


rule scTagger_lr_seg:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
    output:
        tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz",
    threads: 32
    benchmark:
        f"{time_d}/{{sample}}/scTagger_lr_seg.benchmark"
    shell:
        "scTagger.py extract_lr_bc"
        " -r {input.reads}"
        " -o {output.tsv}"
        " -t {threads}"


rule minimap_cdna_for_badread_models:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "cDNA"),        
    output:
        paf=f"{preproc_d}/badread/{{sample}}.badread.cDNA.paf",
    benchmark:
        f"{time_d}/{{sample}}/minimap2_cdna.benchmark"
    threads: 32
    shell:
        "minimap2"
        " -t {threads}"
        " -x map-ont"
        " -c"
        " -o {output.paf}"
        " {input.ref}"
        " {input.reads}"


rule badread_error_model:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "cDNA"),        
        paf=f"{preproc_d}/badread/{{sample}}.badread.cDNA.paf",
    output:
        model=f"{preproc_d}/models/badread/{{sample}}.error.gz",
    shell:
        "badread error_model"
        " --reads {input.reads}"
        " --reference {input.ref}"
        " --alignment {input.paf}"
        " --max_alignments 250000"
        " > {output.model}"


rule badread_qscore_model:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "cDNA"),        
        paf=f"{preproc_d}/badread/{{sample}}.badread.cDNA.paf",
    output:
        model=f"{preproc_d}/models/badread/{{sample}}.qscore.gz",
    shell:
        "badread qscore_model"
        " --reads {input.reads}"
        " --reference {input.ref}"
        " --alignment {input.paf}"
        " --max_alignments 250000"
        " > {output.model}"
