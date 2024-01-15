import sys
import re
from collections import Counter
from datetime import datetime
from typing import NamedTuple

if len(config) == 0:

    configfile: "config.yaml"


outpath = config["outpath"]
preproc_d = f"{outpath}/preprocess"
TS_d = f"{outpath}/TS"
exprmnts_re = "|".join([re.escape(x) for x in config["TS_experiments"]])

DEBUG = False

for sample in config["samples"]:
    for mtype in ["Tsb", "Trc", "Seq"]:
        if mtype not in config["models"]:
            config["models"][mtype] = dict()
        if sample not in config["models"][mtype]:
            config["models"][mtype][sample] = {
                "sample": sample,
                "params": "",
            }

model_details_t = NamedTuple(
    "model_details_t",
    [
        ("name", str),
        ("sample", str),
        ("inputs", list),
        ("outputs", list),
        ("params_build", str),
        ("params_run", str),
    ],
)


def get_model_details(mtype, name):
    model_dict = config["models"][mtype][name]
    sample = model_dict["sample"]
    inputs = list()
    outputs = list()
    params_build = list()
    params_run = list()
    if "params" in model_dict:
        params_build.append(model_dict["params"])
    if mtype == "Tsb":
        assert set(model_dict.keys()) <= {"sample", "params", "cb-txt", "lr-bc"}
        # Inputs / Build params
        paf = get_sample_paf(sample, "cDNA")
        params_build.append(f"-p {paf}")
        inputs.append(paf)
        Xpr_tsv = f"{preproc_d}/tksm_abundance/{name}.Xpr.tsv"
        params_build.append(f"-o {Xpr_tsv}")
        if "cb-txt" in model_dict:
            cb_txt = get_barcode_whitelist(model_dict["cb-txt"])
            params_build.append(f"--cb-txt {cb_txt}")
            inputs.append(cb_txt)
        if "lr-bc" in model_dict:
            lr_matches_tsv = f"{preproc_d}/scTagger/{sample}/{sample}.lr_matches.tsv.gz"
            params_build.append(f"--lr-bc {lr_matches_tsv}")
            inputs.append(lr_matches_tsv)
        # Outputs / Run params
        params_run.append(f"-a {Xpr_tsv}")
        outputs.append(Xpr_tsv)
    elif mtype == "Trc":
        # Inputs / Build params
        assert set(model_dict.keys()) <= {"sample", "params"}
        paf = get_sample_paf(sample, "cDNA")
        params_build.append(f"-i {paf}")
        inputs.append(paf)
        out_prefix = f"{preproc_d}/models/truncate/{name}"
        params_build.append(f"-o {out_prefix}")
        # Outputs / Run params
        kde_files = [
            f"{out_prefix}.{x}"
            for x in ("X_idxs.npy", "Y_idxs.npy", "grid.npy", "sider.tsv")
        ]
        outputs.extend(kde_files)
        params_run.append(f"--kde-model {','.join(kde_files)}")
    elif mtype == "Seq":
        # Inputs / Build params
        assert set(model_dict.keys()) <= {"sample", "params"}
        paf = get_sample_paf(sample, "badread.cDNA")
        params_build.append(f"--alignment {paf}")
        inputs.append(paf)
        fastqs = get_sample_fastqs(sample)
        params_build.append(f"--reads {','.join(fastqs)}")
        inputs.extend(fastqs)
        ref = get_sample_ref(sample, "cDNA")
        params_build.append(f"--reference {ref}")
        inputs.append(ref)
        # Outputs / Run params
        error = f"{preproc_d}/models/badread/{name}.error.gz"
        params_run.append(f"--badread-error-model {error}")
        outputs.append(error)
        qscore = f"{preproc_d}/models/badread/{name}.qscore.gz"
        params_run.append(f"--badread-qscore-model {qscore}")
        outputs.append(qscore)
    else:
        raise ValueError(f"Invalid model type! {mtype}")
    return model_details_t(
        name=name,
        sample=sample,
        inputs=inputs,
        outputs=outputs,
        params_build=" ".join(params_build),
        params_run=" ".join(params_run),
    )


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
        "Shf",
        "Uns",
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
        # If 1st step is transcribe, then return the reference of its model's sample
        if rule_name == "Tsb":
            return get_sample_ref_names(models["Tsb", step["model"]].sample)
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


def get_sample_paf(sample, ref_type):
    return f"{preproc_d}/minimap2/{sample}.{ref_type}.paf"


def get_barcode_whitelist(name):
    return config["refs"]["barcodes"][name]


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


def get_model(wc, rule_name):
    if hasattr(wc, "prefix"):
        prefix = f"{wc.prefix}.{rule_name}"
    else:
        prefix = rule_name
    step = get_step(wc.exprmnt, prefix)
    if "model" in step:
        model_name = step["model"]
        return models[rule_name, model_name]
    else:
        return model_details_t(
            name="",
            sample="",
            inputs=list(),
            outputs=list(),
            params_build="",
            params_run="",
        )


models = dict()
for model_type in config["models"]:
    for model_name in config["models"][model_type]:
        models[model_type, model_name] = get_model_details(
            mtype=model_type,
            name=model_name,
        )


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
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        fastas=lambda wc: get_sample_ref(wc.exprmnt, "DNA"),
        model=lambda wc: get_model(wc=wc, rule_name="Seq").outputs,
    output:
        fastq=f"{TS_d}/{{exprmnt}}/{{prefix}}.Seq.fastq",
    threads: 32
    params:
        config=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Seq")["params"],
        binary=config["exec"]["tksm"],
        fastas=lambda wc: get_sample_ref(wc.exprmnt, "DNA"),
        model=lambda wc: get_model(wc=wc, rule_name="Seq").params_run,
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} sequence"
        " -i {input.mdf}"
        " --references {params.fastas}"
        " -o {output.fastq}"
        " --threads {threads}"
        " {params.config}"
        " {params.model}"


rule filter:
    input:
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Flt.mdf"),
    params:
        config=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Flt")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} filter"
        " -i {input.mdf}"
        " -t {output.mdf}"
        " {params.config}"


rule truncate:
    input:
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
        model=lambda wc: get_model(wc=wc, rule_name="Trc").outputs,
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Trc.mdf"),
    params:
        config=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Trc")["params"],
        model=lambda wc: get_model(wc=wc, rule_name="Trc").params_run,
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} truncate"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.config}"
        " {params.model}"


rule unsegment:
    input:
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Uns.mdf"),
    params:
        config=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Uns")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} unsegment"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.config}"


rule shuffle:
    input:
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Shf.mdf"),
    params:
        config=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Shf")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} shuffle"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.config}"


rule flip:
    input:
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Flp.mdf"),
    params:
        config=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Flp")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} flip"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.config}"


rule pcr:
    input:
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.PCR.mdf"),
    params:
        config=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.PCR")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} pcr"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.config}"


rule tag:
    input:
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.Tag.mdf"),
    params:
        config=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.Tag")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} tag"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.config}"


rule single_cell_barcoder:
    input:
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.SCB.mdf"),
    params:
        config=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.SCB")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} scb"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.config}"


rule polyA:
    input:
        mdf=f"{TS_d}/{{exprmnt}}/{{prefix}}.mdf",
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/{{prefix}}.plA.mdf"),
    params:
        config=lambda wc: get_step(wc.exprmnt, f"{wc.prefix}.plA")["params"],
        binary=config["exec"]["tksm"],
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} polyA"
        " -i {input.mdf}"
        " -o {output.mdf}"
        " {params.config}"


### Entry rules ###
rule transcribe:
    input:
        model=lambda wc: get_model(wc=wc, rule_name="Tsb").outputs,
        gtf=lambda wc: get_sample_ref(wc.exprmnt, "GTF"),
    output:
        mdf=pipe(f"{TS_d}/{{exprmnt}}/Tsb.mdf"),
    params:
        config=lambda wc: get_step(wc.exprmnt, f"Tsb")["params"],
        binary=config["exec"]["tksm"],
        model=lambda wc: get_model(wc=wc, rule_name="Tsb").params_run,
    wildcard_constraints:
        exprmnt=exprmnts_re,
    shell:
        "{params.binary} transcribe"
        " -g {input.gtf}"
        " -o {output.mdf}"
        " {params.model}"
        " {params.config}"


if config["enable_piping"] == False:

    rule merge:
        input:
            mdfs=get_merge_mdf_input,
        output:
            mdf=pipe(f"{TS_d}/{{exprmnt}}/Mrg.mdf"),
        shell:
            "cat {input.mdfs} > {output.mdf}"

else:

    rule merge:
        input:
            script="py/mdf_cat.py",
            mdfs=lambda wc: merge_to_numbered_sources[wc.exprmnt],
        output:
            mdf=pipe(f"{TS_d}/{{exprmnt}}/Mrg.mdf"),
        shell:
            "python {input.script} {input.mdfs}  {output.mdf}"


### Model rules ###
rule model_transcribe:
    input:
        model=lambda wc: models["Tsb", wc.model_name].inputs,
    output:
        model=f"{preproc_d}/tksm_abundance/{{model_name}}.Xpr.tsv",
    params:
        binary=config["exec"]["tksm"],
        model=lambda wc: models["Tsb", wc.model_name].params_build,
    shell:
        "{params.binary} abundance {params.model}"


rule model_truncation:
    input:
        model=lambda wc: models["Trc", wc.model_name].inputs,
    output:
        model=[
            f"{preproc_d}/models/truncate/{{model_name}}.{x}"
            for x in ("X_idxs.npy", "Y_idxs.npy", "grid.npy", "sider.tsv")
        ],
    params:
        binary=config["exec"]["tksm"],
        model=lambda wc: models["Trc", wc.model_name].params_build,
    threads: 32
    shell:
        "{params.binary} model-truncation {params.model} --threads {threads}"


rule model_sequence:
    input:
        model=lambda wc: models["Seq", wc.model_name].inputs,
    output:
        error_model=f"{preproc_d}/models/badread/{{model_name}}.error.gz",
        qscore_model=f"{preproc_d}/models/badread/{{model_name}}.qscore.gz",
    params:
        model=lambda wc: models["Seq", wc.model_name].params_build,
    shell:
        "badread qscore_model {params.model} > {output.error_model}"
        " && "
        "badread error_model {params.model} > {output.qscore_model}"


### Preprocessing rules ###
rule minimap_cdna:
    input:
        reads=lambda wc: get_sample_fastqs(wc.sample),
        ref=lambda wc: get_sample_ref(wc.sample, "cDNA"),
    output:
        paf=f"{preproc_d}/minimap2/{{sample}}.cDNA.paf",
    threads: 32
    shell:
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
    output:
        lr_tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_matches.tsv.gz",
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
        wl=lambda wc: config["refs"][config["samples"][wc.sample]["cb_wl"]],
    output:
        tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.bc_whitelist.tsv.gz",
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
    threads: 32
    shell:
        "minimap2"
        " -t {threads}"
        " -x map-ont"
        " -c"
        " -o {output.paf}"
        " {input.ref}"
        " {input.reads}"


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
