import sys

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


def component_idx(prefix):
    return len(prefix.split("."))


def fastas_for_TS_sequence(wc):
    fastas = list()
    sample = exprmnt_sample(wc.exprmnt)
    fastas.append(config["refs"][get_sample_ref(sample)]["DNA"])
    for idx, component in enumerate(config["TS_experiments"][wc.exprmnt]["pipeline"]):
        component = list(component.keys())[0]
        if component in ["plA", "SCB", "UMI"]:
            prefix = ".".join(
                [
                    list(c.keys())[0]
                    for c in config["TS_experiments"][wc.exprmnt]["pipeline"][: idx + 1]
                ]
            )
            fastas.append(
                f"{TS_d}/{wc.exprmnt}/{prefix}.fasta",
            )
    return fastas


def experiment_prefix(exprmnt):
    prefix = list()
    for component in config["TS_experiments"][exprmnt]["pipeline"]:
        prefix.append(list(component)[0])
    return ".".join(prefix)


def get_sample_ref(name):
    if name in config["samples"]:
        return config["samples"][name]["ref"]
    elif name in config["TS_experiments"]:
        return get_sample_ref(config["TS_experiments"][name]["sample"])
    else:
        print(f"Invalid experiment/sample name! {name}")
        1 / 0


def get_sample_fastqs(name):
    if name in config["samples"]:
        sample = name
        return config["samples"][sample]["fastq"]
    elif name in config["TS_experiments"]:
        exprmnt = name
        return [f"{TS_d}/{exprmnt}/{experiment_prefix(exprmnt)}.fastq"]
    else:
        print(f"Invalid experiment/sample name! {name}")
        1 / 0


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


rule sequence:
    input:
        binary=config["exec"]["tksm"],
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        fastas=fastas_for_TS_sequence,
        qscore_model=lambda wc: f"{preproc_d}/models/badread/{exprmnt_sample(wc.exprmnt)}.qscore.gz",
        error_model=lambda wc: f"{preproc_d}/models/badread/{exprmnt_sample(wc.exprmnt)}.error.gz",
    output:
        fastq=f"{TS_d}/{{exprmnt}}/{{prefix}}.Seq.fastq",
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.Seq.benchmark"
    params:
        other=lambda wc: config["TS_experiments"][wc.exprmnt]["pipeline"][
            component_idx(wc.prefix)
        ]["Seq"],
    threads: 32
    shell:
        "{input.binary} sequencer"
        " -i {input.mdf}"
        " --references {input.fastas}"
        " -o {output.fastq}"
        " --threads {threads}"
        " --badread-error-model={input.error_model}"
        " --badread-qscore-model={input.qscore_model}"
        " {params.other}"


rule mirror:
    input:
        "{anything}",
    output:
        temp("{anything}.mirror"),
    shell:
        "cat {input} > {output}"


rule truncate:
    input:
        binary=config["exec"]["tksm"],
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        x=lambda wc: f"{preproc_d}/truncate_kde/{exprmnt_sample(wc.exprmnt)}.X_idxs.npy",
        y=lambda wc: f"{preproc_d}/truncate_kde/{exprmnt_sample(wc.exprmnt)}.Y_idxs.npy",
        g=lambda wc: f"{preproc_d}/truncate_kde/{exprmnt_sample(wc.exprmnt)}.grid.npy",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Trc.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.Trc.benchmark"
    params:
        lambda wc: config["TS_experiments"][wc.exprmnt]["pipeline"][
            component_idx(wc.prefix)
        ]["Trc"],
    shell:
        "{input.binary} truncate"
        " -i {input.mdf}"
        " --kde={input.g},{input.x},{input.y}"
        " -o {output.mdf}"
        " {params}"


rule pcr:
    input:
        binary=config["exec"]["tksm"],
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.PCR.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.PCR.benchmark"
    params:
        lambda wc: config["TS_experiments"][wc.exprmnt]["pipeline"][
            component_idx(wc.prefix)
        ]["PCR"],
    shell:
        "{input.binary} pcr"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params}"


rule umi:
    input:
        binary=config["exec"]["tksm"],
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.UMI.mdf"),
        fasta=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.UMI.fasta"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.UMI.benchmark"
    params:
        lambda wc: config["TS_experiments"][wc.exprmnt]["pipeline"][
            component_idx(wc.prefix)
        ]["UMI"],
    shell:
        "{input.binary} tag"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " -f {output.fasta}"
        " {params}"


rule single_cell_barcoder:
    input:
        binary=config["exec"]["tksm"],
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.SCB.mdf"),
        fasta=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.SCB.fasta"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.SCB.benchmark"
    params:
        lambda wc: config["TS_experiments"][wc.exprmnt]["pipeline"][
            component_idx(wc.prefix)
        ]["SCB"],
    shell:
        "{input.binary} single-cell-barcoder"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " -f {output.fasta}"
        " {params}"


rule polyA:
    input:
        binary=config["exec"]["tksm"],
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.plA.mdf"),
        fasta=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.plA.fasta"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.plA.benchmark"
    params:
        lambda wc: config["TS_experiments"][wc.exprmnt]["pipeline"][
            component_idx(wc.prefix)
        ]["plA"],
    shell:
        "{input.binary} polyA"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " -f {output.fasta}"
        " {params}"


rule splicer:
    input:
        binary=config["exec"]["tksm"],
        tsv=f"{TS_d}/{{exprmnt}}/{{prefix}}.tsv",
        gtf=lambda wc: config["refs"][get_sample_ref(exprmnt_sample(wc.exprmnt))]["GTF"],
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Spc.mdf"),
    benchmark:
        f"{time_d}/{{exprmnt}}/{{prefix}}.Spc.benchmark"
    params:
        lambda wc: config["TS_experiments"][wc.exprmnt]["pipeline"][
            component_idx(wc.prefix)
        ]["Spc"],
    shell:
        "{input.binary} splicer"
        " -a {input.tsv}"
        " -g {input.gtf}"
        " -o {output.mdf}"
        " {params}"


rule abundance:
    input:
        binary=config["exec"]["tksm"],
        paf=lambda wc: f"{preproc_d}/minimap2/{exprmnt_sample(wc.exprmnt)}.cDNA.paf",
    output:
        tsv=f"{TS_d}/{{exprmnt}}/Xpr.tsv",
    benchmark:
        f"{time_d}/{{exprmnt}}/Xpr.benchmark"
    shell:
        "{input.binary} abundance"
        " -p {input.paf}"
        " -o {output.tsv}"


rule abundance_sc:
    input:
        binary=config["exec"]["tksm"],
        paf=lambda wc: f"{preproc_d}/minimap2/{exprmnt_sample(wc.exprmnt)}.cDNA.paf",
        lr_matches=lambda wc: f"{preproc_d}/scTagger/{exprmnt_sample(wc.exprmnt)}/{exprmnt_sample(wc.exprmnt)}.lr_matches.tsv.gz",
    output:
        tsv=f"{TS_d}/{{exprmnt}}/Xpr_sc.tsv",
    benchmark:
        f"{time_d}/{{exprmnt}}/Xpr_sc.benchmark"
    shell:
        "{input.binary} abundance"
        " -p {input.paf}"
        " -m {input.lr_matches}"
        " -o {output.tsv}"


rule truncate_kde:
    input:
        binary=config["exec"]["tksm"],
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
    output:
        x=f"{preproc_d}/truncate_kde/{{sample}}.X_idxs.npy",
        y=f"{preproc_d}/truncate_kde/{{sample}}.Y_idxs.npy",
        g=f"{preproc_d}/truncate_kde/{{sample}}.grid.npy",
    benchmark:
        f"{time_d}/{{sample}}/truncate_kde.benchmark"
    params:
        out_prefix=f"{preproc_d}/truncate_kde/{{sample}}",
    threads: 32
    shell:
        "{input.binary} kde"
        " -i {input.paf}"
        " -o {params.out_prefix}"
        " --threads {threads}"


rule minimap_cdna:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wildcards: config["refs"][get_sample_ref(wildcards.sample)]["cDNA"],
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
        ref=lambda wildcards: config["refs"][get_sample_ref(wildcards.sample)]["cDNA"],
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
        ref=lambda wildcards: config["refs"][get_sample_ref(wildcards.sample)]["cDNA"],
        paf=f"{preproc_d}/badread/{{sample}}.badread.cDNA.paf",
    output:
        model=f"{preproc_d}/models/badread/{{sample}}.error.gz",
    shell:
        "badread error_model"
        " --reads {input.reads}"
        " --reference {input.ref}"
        " --alignment {input.paf}"
        " > {output.model}"


rule badread_qscore_model:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wildcards: config["refs"][get_sample_ref(wildcards.sample)]["cDNA"],
        paf=f"{preproc_d}/badread/{{sample}}.badread.cDNA.paf",
    output:
        model=f"{preproc_d}/models/badread/{{sample}}.qscore.gz",
    shell:
        "badread qscore_model"
        " --reads {input.reads}"
        " --reference {input.ref}"
        " --alignment {input.paf}"
        " --max_alignments 100000"
        " > {output.model}"
