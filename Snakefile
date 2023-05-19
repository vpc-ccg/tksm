import sys
import re
from collections import Counter
from datetime import datetime

if len(config) == 0:

    configfile: "config.yaml"


outpath = config["outpath"]
preproc_d = f"{outpath}/preprocess"
TS_d = f"{outpath}/TS"
time_tsv = f"{outpath}/time.tsv"
exprmnts_re = "|".join([re.escape(x) for x in config["TS_experiments"]])

DEBUG = True


def exprmnt_final_file(exprmnt):
    prefix = [list(step)[0] for step in config["TS_experiments"][exprmnt]["pipeline"]]
    if prefix[-1] in ["Seq"]:
        prefix.append("fastq")
    elif prefix[-1] in [
        "Tsb",
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
        # If 1st step is transcribe, then return its model's reference
        if rule_name == "Tsb":
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
        f"{preproc_d}/models/truncate/{model}.{x}.npy"
        for x in ["grid", "X_idxs", "Y_idxs"]
    ]
    return kde_input


def format_gnu_time_string(
    process="",
    exprmnt="{wildcards.exprmnt}",
    prefix="{wildcards.prefix}",
    threads="{threads}",
):
    if config["benchmark_time"] == False:
        return ""
    fields = list()
    fields.append((f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", "%s", ""))
    fields.append((f"{exprmnt}", "%s", ""))
    fields.append((f"{prefix}", "%s", ""))
    fields.append((f"{process}", "%s", ""))
    fields.append(("%e", "%.2f", "/60"))
    fields.append(("%U", "%.2f", "/60"))
    fields.append(("%M", "%.2f", "/(1024*1024)"))
    fields.append((f"{threads}", "%d", ""))
    fields.append(("%S", "%.2f", "/60"))
    fields.append((f"{config['enable_piping']}", "%s", ""))

    time_format = ",".join([x[0] for x in fields])
    printf_format = ",".join([x[1] for x in fields]) + "\\n"
    printf_args = ",".join([f"${i}{x[2]}" for i, x in enumerate(fields, start=1)])

    awk_cmd = (
        'awk \'BEGIN{{FS=","}} {{printf "' + printf_format + '",' + printf_args + "}}'"
    )
    return f'$(which time) -f "{time_format}" -o >({awk_cmd} >> {{input.time}}) '


rule all:
    input:
        [
            exprmnt_final_file(exprmnt)
            for exprmnt in config["TS_experiments"]
            if exprmnt_final_file(exprmnt).endswith("fastq")
        ],


if config["enable_piping"] == True:
    sys.stderr.write("Piping enabled for TKSM!\n")
    merge_source_mdf_counter = Counter()
    merge_to_numbered_sources = dict()

    for exprmnt in config["TS_experiments"]:
        step, details = tuple(config["TS_experiments"][exprmnt]["pipeline"][0].items())[
            0
        ]
        if not step == "Mrg":
            continue
        merge_to_numbered_sources[exprmnt] = list()
        for source_exprmnt in details["sources"]:
            source_mdf = exprmnt_final_file(source_exprmnt)
            merge_to_numbered_sources[exprmnt].append(
                f"{source_mdf}.{merge_source_mdf_counter[source_mdf]}"
            )
            merge_source_mdf_counter[source_mdf] += 1

    for mdf, count in merge_source_mdf_counter.items():
        rule:
            input:
                mdf=mdf,
                script="py/mdf_tee.py",
            output:
                [pipe(f"{mdf}.{number}") for number in range(count)],
            shell:
                "python {input.script} {input.mdf} {output}"

else:

    def pipe(X):
        return X


rule sequence:
    input:
        obj=["build/obj/sequence.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        fastas=lambda wc: get_sample_ref(wc.exprmnt, "DNA"),
        qscore_model=lambda wc: get_sequencer_model_input(wc, "qscore"),
        error_model=lambda wc: get_sequencer_model_input(wc, "error"),
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        fastq=f"{TS_d}/{{exprmnt}}/{{prefix}}.Seq.fastq",
    threads: 32
    params:
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Seq")["params"],
        binary=config["exec"]["tksm"],
        fastas=lambda wc: get_sample_ref(wc.exprmnt, "DNA"),
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        f"{format_gnu_time_string(process='sequence')}"
        "{params.binary} sequence"
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
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Flt.mdf"),
    params:
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Flt")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        f"{format_gnu_time_string(process='filter')}"
        "{params.binary} filter"
        " -i {input.mdf}"
        " -t {output.mdf}"
        " {params.other}"


rule truncate:
    input:
        obj=["build/obj/truncate.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        kde=get_kde_model_input,
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Trc.mdf"),
    params:
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Trc")["params"],
        kde=lambda wc: ",".join(get_kde_model_input(wc)),
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        f"{format_gnu_time_string(process='truncate')}"
        "{params.binary} truncate"
        " -i {input.mdf}"
        " --kde-model={params.kde}"
        " -o {output.mdf}"
        " {params.other}"


rule unsegment:
    input:
        obj=["build/obj/strand_man.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Uns.mdf"),
    params:
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Uns")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        f"{format_gnu_time_string(process='unsegment')}"
        "{params.binary} unsegment"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.other}"


rule shuffle:
    input:
        obj=["build/obj/strand_man.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Shf.mdf"),
    params:
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Shf")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        f"{format_gnu_time_string(process='shuffle')}"
        "{params.binary} shuffle"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.other}"


rule flip:
    input:
        obj=["build/obj/strand_man.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Flp.mdf"),
    params:
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Flp")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        f"{format_gnu_time_string(process='flip')}"
        "{params.binary} flip"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.other}"


rule pcr:
    input:
        obj=["build/obj/pcr.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.PCR.mdf"),
    params:
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.PCR")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        f"{format_gnu_time_string(process='pcr')}"
        "{params.binary} pcr"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.other}"


rule tag:
    input:
        obj=["build/obj/tag.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Tag.mdf"),
    params:
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Tag")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        f"{format_gnu_time_string(process='tag')}"
        "{params.binary} tag"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.other}"


rule single_cell_barcoder:
    input:
        obj=["build/obj/scb.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.SCB.mdf"),
    params:
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.SCB")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        f"{format_gnu_time_string(process='scb')}"
        "{params.binary} scb"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.other}"


rule polyA:
    input:
        obj=["build/obj/polyA.o", "build/obj/tksm.o"] if DEBUG else list(),
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.plA.mdf"),
    params:
        other=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.plA")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        f"{format_gnu_time_string(process='polyA')}"
        "{params.binary} polyA"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.other}"


### Entry rules ###
rule transcribe:
    input:
        obj=["build/obj/transcribe.o", "build/obj/tksm.o"] if DEBUG else list(),
        tsv=lambda wc: f"{preproc_d}/tksm_abundance/{get_step(wc.exprmnt, 'Tsb')['model']}.{get_step(wc.exprmnt, 'Tsb')['mode']}.tsv",
        gtf=lambda wc: get_sample_ref(wc.exprmnt, "GTF"),
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/Tsb.mdf"),
    params:
        other=lambda wc: get_step(wc.exprmnt, f"Tsb")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        f"{format_gnu_time_string(process='transcribe', prefix='')}"
        "{params.binary} transcribe"
        " -a {input.tsv}"
        " -g {input.gtf}"
        " -o {output.mdf}"
        " {params.other}"


if config["enable_piping"] == False:
    rule merge:
        input:
            mdfs=get_merge_mdf_input,
            time=ancient(time_tsv) if config["benchmark_time"] else list(),
        output:
            mdf=pipe(f"{TS_d}/{{exprmnt}}/Mrg.mdf"),
        shell:
            f"{format_gnu_time_string(process='merge', prefix='')}"
            "cat {input.mdfs} > {output.mdf}"
else:
    rule merge:
        input:
            script="py/mdf_cat.py",
            mdfs=lambda wc: merge_to_numbered_sources[wc.exprmnt],
            time=ancient(time_tsv) if config["benchmark_time"] else list(),
        output:
            mdf=pipe(f"{TS_d}/{{exprmnt}}/Mrg.mdf"),
        shell:
            f"{format_gnu_time_string(process='merge', prefix='')}"
            "python {input.script} {input.mdfs}  {output.mdf}"


### Preprocessing rules ###
rule abundance:
    input:
        obj=["build/obj/abundance.o", "build/obj/tksm.o"] if DEBUG else list(),
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        tsv=f"{preproc_d}/tksm_abundance/{{sample}}.Xpr.tsv",
    params:
        binary=config["exec"]["tksm"],
    shell:
        f"{format_gnu_time_string(process='abundance', exprmnt='{wildcards.sample}', prefix='')}"
        "{params.binary} abundance"
        " -p {input.paf}"
        " -o {output.tsv}"


rule abundance_sc:
    input:
        obj=["build/obj/abundance.o", "build/obj/tksm.o"] if DEBUG else list(),
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
        lr_matches=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_matches.tsv.gz",
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        tsv=f"{preproc_d}/tksm_abundance/{{sample}}.Xpr_sc.tsv",
    params:
        binary=config["exec"]["tksm"],
    shell:
        f"{format_gnu_time_string(process='abundance_sc', exprmnt='{wildcards.sample}', prefix='')}"
        "{params.binary} abundance"
        " -p {input.paf}"
        " -m {input.lr_matches}"
        " -o {output.tsv}"


rule model_truncation:
    input:
        obj=["build/obj/model_truncation.o", "build/obj/tksm.o"] if DEBUG else list(),
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        x=f"{preproc_d}/models/truncate/{{sample}}.X_idxs.npy",
        y=f"{preproc_d}/models/truncate/{{sample}}.Y_idxs.npy",
        g=f"{preproc_d}/models/truncate/{{sample}}.grid.npy",
    params:
        out_prefix=f"{preproc_d}/models/truncate/{{sample}}",
        binary=config["exec"]["tksm"],
    threads: 32
    shell:
        f"{format_gnu_time_string(process='model_truncation', exprmnt='{wildcards.sample}', prefix='')}"
        "{params.binary} model-truncation"
        " -i {input.paf}"
        " -o {params.out_prefix}"
        " --threads {threads}"


rule minimap_cdna:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "cDNA"),
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
    threads: 32
    shell:
        f"{format_gnu_time_string(process='minimap_cdna', exprmnt='{wildcards.sample}', prefix='')}"
        "minimap2"
        " -t {threads}"
        " -x map-ont"
        " -p 0.0"
        " -c --eqx"
        " -o {output.paf}"
        " {input.ref}"
        " {input.reads}"


rule scTagger_match:
    input:
        lr_tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz",
        wl_tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.bc_whitelist.tsv.gz",
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        lr_tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_matches.tsv.gz",
    threads: 32
    shell:
        f"{format_gnu_time_string(process='scTagger_match', exprmnt='{wildcards.sample}', prefix='')}"
        "scTagger.py match_trie"
        " -lr {input.lr_tsv}"
        " -sr {input.wl_tsv}"
        " -o {output.lr_tsv}"
        " -t {threads}"


rule scTagger_extract_bc:
    input:
        tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz",
        wl=config["refs"]["10x_bc"],
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.bc_whitelist.tsv.gz",
    shell:
        f"{format_gnu_time_string(process='scTagger_extract_bc', exprmnt='{wildcards.sample}', prefix='')}"
        "scTagger.py extract_sr_bc_from_lr"
        " -i {input.tsv}"
        " -wl {input.wl}"
        " -o {output.tsv}"


rule scTagger_lr_seg:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz",
    threads: 32
    shell:
        f"{format_gnu_time_string(process='scTagger_lr_seg', exprmnt='{wildcards.sample}', prefix='')}"
        "scTagger.py extract_lr_bc"
        " -r {input.reads}"
        " -o {output.tsv}"
        " -t {threads}"


rule minimap_cdna_for_badread_models:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "cDNA"),
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        paf=f"{preproc_d}/badread/{{sample}}.badread.cDNA.paf",
    threads: 32
    shell:
        f"{format_gnu_time_string(process='minimap_cdna_for_badread_models', exprmnt='{wildcards.sample}', prefix='')}"
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
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        model=f"{preproc_d}/models/badread/{{sample}}.error.gz",
    shell:
        f"{format_gnu_time_string(process='badread_error_model', exprmnt='{wildcards.sample}', prefix='')}"
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
        time=ancient(time_tsv) if config["benchmark_time"] else list(),
    output:
        model=f"{preproc_d}/models/badread/{{sample}}.qscore.gz",
    shell:
        f"{format_gnu_time_string(process='badread_qscore_model', exprmnt='{wildcards.sample}', prefix='')}"
        "badread qscore_model"
        " --reads {input.reads}"
        " --reference {input.ref}"
        " --alignment {input.paf}"
        " --max_alignments 250000"
        " > {output.model}"


rule make_time:
    output:
        time_tsv,
    shell:
        'echo "Timestamp,Experiment,Prefix,Process,Real time (min),User time (min),Memory (GB),Threads,System time (min),Piped?" > {output}'


rule cat_refs:
    input:
        refs=lambda wc: [
            config["refs"][ref_name][wc.ref_type]
            for ref_name in wc.ref_names.split(":")
        ],
    output:
        ref=f"{preproc_d}/refs/{{ref_names}}.{{ref_type}}.{{file_type}}",
    run:
        for ref in input.refs:
            if ref.endswith(".gz"):
                shell("zcat {ref} >> {output.ref}")
            else:
                shell("cat {ref} >> {output.ref}")
