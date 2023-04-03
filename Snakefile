import sys
import re

if len(config) == 0:

    configfile: "config.yaml"


outpath = config["outpath"]
preproc_d = f"{outpath}/preprocess"
TS_d = f"{outpath}/TS"
time_d = f"{outpath}/time"
exprmnts_re = "|".join([re.escape(x) for x in config["TS_experiments"]])
DEBUG = True

if DEBUG:

    def pipe(X):
        return X

    def temp(X):
        return X


def exprmnt_final_file(exprmnt):
    prefix = [list(step)[0] for step in config["TS_experiments"][exprmnt]["pipeline"]]
    if prefix[-1] in ["Seq"]:
        prefix.append("fastq")
    elif prefix[-1] in [
        "Spc",
        "Flt",
        "PCR",
        "plA",
        "SCB",
        "Tag",
        "Flp",
        "Trc",
    ]:
        prefix.append("mdf")
    else:
        raise ValueError(f"Invalid terminal pipeline step! {prefix[-1]}")
    prefix = ".".join(prefix)
    final_file = f"{TS_d}/{exprmnt}/{prefix}"
    return final_file


def get_sample_ref_names(sample):
    # Check if sample is real
    if sample in config["samples"]:
        return [config["samples"][sample]["ref"]]
    # If not, then sample must be a TS experiment
    if sample in config["TS_experiments"]:
        step = config["TS_experiments"][sample]["pipeline"][0]
        rule_name = list(step)[0]
        step = step[rule_name]
        # If 1st step is splicer, then return its model's reference
        if rule_name == "Spc":
            return get_sample_ref_names(step["model"])
        # If 1st step is merge, then return the references of its sources
        if rule_name == "Mrg":
            ref_names = set()
            for source in step["sources"]:
                ref_names.update(get_sample_ref_names(source))
            ref_names = sorted(ref_names)
            return ref_names
        raise ValueError(f"Invalid 1st rule ({rule_name}) for sample ({sample})!")
    raise ValueError(f"Invalid sample ({sample})!")


def get_sample_refs(sample, ref_type):
    refs = list()
    for ref_name in get_sample_ref_names(sample):
        refs.append(config["refs"][ref_name][ref_type])
    return refs


def get_sample_fastqs(name):
    if name in config["samples"]:
        sample = name
        return config["samples"][sample]["fastq"]
    if name in config["TS_experiments"]:
        exprmnt = name
        fastq = exprmnt_final_file(exprmnt)
        assert fastq.endswith(".fastq")
        return [fastq]
    raise ValueError(f"Invalid experiment/sample name! {name}")


def get_step(exprmnt, prefix):
    idx = len(prefix.split(".")) - 1
    step = config["TS_experiments"][exprmnt]["pipeline"][idx]
    rule_name = list(step)[0]
    return step[rule_name]


def get_merge_mdf_input(wc):
    step = get_step(wc.exprmnt, "Mrg")
    mdfs = list()
    for source in step["sources"]:
        mdf = exprmnt_final_file(source)
        assert mdf.endswith(".mdf")
        mdfs.append(mdf)
    return mdfs


def get_sequencer_model_input(wc, model_type):
    step = get_step(wc.exprmnt, f"{wc.prefix}.Seq")
    model = step["model"]
    if model in config["samples"] or model in config["TS_experiments"]:
        return f"{preproc_d}/models/badread/{model}.{model_type}.gz"
    else:
        return list()


def get_kde_model_input(wc):
    step = get_step(wc.exprmnt, f"{wc.prefix}.Kde")
    model = step["model"]
    kde_input = [
        f"{preproc_d}/truncate_kde/{model}.{x}.npy"
        for x in ["grid", "X_idxs", "Y_idxs"]
    ]
    return kde_input


rule all:
    input:
        [exprmnt_final_file(exprmnt) for exprmnt in config["TS_experiments"]],


rule sequencer:
    input:
        obj=["build/obj/sequencer.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        fastas=lambda wc: get_sample_refs(wc.exprmnt, "DNA"),
        qscore_model=lambda wc: get_sequencer_model_input(wc, "qscore"),
        error_model=lambda wc: get_sequencer_model_input(wc, "error"),
    output:
        fastq=f"{TS_d}/{{exprmnt}}/{{prefix}}.Seq.fastq",
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.Seq.benchmark"
    threads: 32
    params:
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Seq")["params"],
        binary=config["exec"]["tksm"],
        fastas=lambda wc: ",".join(get_sample_refs(wc.exprmnt, "DNA")),
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} sequencer"
        " -i {input.mdf}"
        " --references {params.fastas}"
        " -o {output.fastq}"
        " --threads {threads}"
        " --badread-error-model={input.error_model}"
        " --badread-qscore-model={input.qscore_model}"
        " {params.other}"


rule filter:
    input:
        obj=["build/obj/filter.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Flt.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.Flt.benchmark"
    params:
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Flt")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} filter"
        " -i {input.mdf}"
        " -t {output.mdf}"
        " {params.other}"


rule truncate:
    input:
        obj=["build/obj/truncate.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        kde=get_kde_model_input,
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Trc.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.Trc.benchmark"
    params:
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Trc")["params"],
        kde=lambda wc: ",".join(get_kde_model_input(wc)),
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} truncate"
        " -i {input.mdf}"
        " --kde={params.kde}"
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
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Flp")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
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
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.PCR")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
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
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Tag")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
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
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.SCB")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
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
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.plA")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} polyA"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.other}"


### Entry rules ###
rule splicer:
    input:
        obj=["build/obj/splicer.o", "build/obj/tksm.o"] if DEBUG else list(),
        tsv=lambda wc: f"{preproc_d}/tksm_abundance/{get_step(wc.exprmnt, 'Spc')['model']}.{get_step(wc.exprmnt, 'Spc')['mode']}.tsv",
        gtfs=lambda wc: get_sample_refs(wc.exprmnt, "GTF"),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/Spc.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/Spc.benchmark"
    params:
        other=lambda wc: get_step(wc.exprmnt, f"Spc")["params"],
        gtfs=lambda wc: ",".join(get_sample_refs(wc.exprmnt, "GTF")),
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} splicer"
        " -a {input.tsv}"
        " -g {params.gtfs}"
        " -o {output.mdf}"
        " {params.other}"


rule merge:
    input:
        mdfs=get_merge_mdf_input,
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/Mrg.mdf"),
    shell:
        "cat {input.mdfs} > {output.mdf}"


### Preprocessing rules ###
rule abundance:
    input:
        obj=["build/obj/abundance.o", "build/obj/tksm.o"] if DEBUG else list(),
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
    output:
        tsv=f"{preproc_d}/tksm_abundance/{{sample}}.Xpr.tsv",
    benchmark:
        f"{time_d}/{{sample}}/Xpr.benchmark"
    params:
        binary=config["exec"]["tksm"],
    shell:
        "{params.binary} abundance"
        " -p {input.paf}"
        " -o {output.tsv}"


rule abundance_sc:
    input:
        obj=["build/obj/abundance.o", "build/obj/tksm.o"] if DEBUG else list(),
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
        lr_matches=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_matches.tsv.gz",
    output:
        tsv=f"{preproc_d}/tksm_abundance/{{sample}}.Xpr_sc.tsv",
    benchmark:
        f"{time_d}/{{sample}}/Xpr_sc.benchmark"
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
        refs=lambda wc: get_sample_refs(wc.sample, "cDNA"),
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
        " <(cat {input.refs})"
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
        refs=lambda wc: get_sample_refs(wc.sample, "cDNA"),
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
        " <(cat {input.refs})"
        " {input.reads}"


rule badread_error_model:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        refs=lambda wc: get_sample_refs(wc.sample, "cDNA"),
        paf=f"{preproc_d}/badread/{{sample}}.badread.cDNA.paf",
    output:
        model=f"{preproc_d}/models/badread/{{sample}}.error.gz",
    shell:
        "badread error_model"
        " --reads {input.reads}"
        " --reference <(cat {input.refs})"
        " --alignment {input.paf}"
        " --max_alignments 250000"
        " > {output.model}"


rule badread_qscore_model:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        refs=lambda wc: get_sample_refs(wc.sample, "cDNA"),
        paf=f"{preproc_d}/badread/{{sample}}.badread.cDNA.paf",
    output:
        model=f"{preproc_d}/models/badread/{{sample}}.qscore.gz",
    shell:
        "badread qscore_model"
        " --reads {input.reads}"
        " --reference <(cat {input.refs})"
        " --alignment {input.paf}"
        " --max_alignments 250000"
        " > {output.model}"
