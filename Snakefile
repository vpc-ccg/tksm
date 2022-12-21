import sys
from itertools import groupby

if len(config)==0:
    configfile: "config.yaml"

outpath = config['outpath']
preproc_d = f'{outpath}/preprocess'
RI_d = f'{outpath}/RI'
nanosim_d = f'{outpath}/nanosim'

rule all:
    input:
        expand(
            f'{RI_d}/{{sample}}/{{experiment}}/sequenced.fastq',
            sample=config['samples'],
            experiment=config['experiments'],
        ),

rule RI_sequence:
    input:
        binary = config['exec']['RI_sequencer'],
        badread = config['exec']['badread'],
        splicer_mdf = f'{RI_d}/{{sample}}/{{experiment}}/splicer.mdf',
        polyA_fasta = f'{RI_d}/{{sample}}/{{experiment}}/polyA.fasta',
    output:
        fastq = f'{RI_d}/{{sample}}/{{experiment}}/sequenced.fastq',
    params:
        dna = config['refs']['DNA'],
        name = lambda wildcards: f'{wildcards.sample}-{wildcards.experiment}-RI',
    threads:
        32
    shell:           
        '{input.binary}'
        ' -m {input.splicer_mdf}'
        ' -o {output.fastq}'
        ' --badread={input.badread}'
        ' --references={params.dna},{input.polyA_fasta}'
        ' -n{params.name}'
        ' --threads={threads}'
        ' --other-br="--tail_noise_model=nanopore"'

rule RI_polyA:
    input:
        binary = config['exec']['RI_polyA'],
        mdf = f'{RI_d}/{{sample}}/{{experiment}}/truncated.mdf',
    output:
        mdf = f'{RI_d}/{{sample}}/{{experiment}}/polyA.mdf',
        fasta = f'{RI_d}/{{sample}}/{{experiment}}/polyA.fasta',
    params:
        distr = lambda wildcards: config['experiments'][wildcards.experiment]['polyA_distr'],
        min_len = lambda wildcards: config['experiments'][wildcards.experiment]['polyA_min_len'],
    shell:
        '{input.binary}'
        ' -i {input.mdf}'
        ' -o {output.mdf}'
        ' --{params.distr}'
        ' --min-length={params.min_len}'
        ' --expand-isoforms'
        ' -a {output.fasta}'

rule RI_trunc:
    input:
        binary = config['exec']['RI_truncate'],
        mdf = f'{RI_d}/{{sample}}/{{experiment}}/splicer.mdf',
        grid = config['exec']['RI_truncate_kde']['grid'],
        labels = config['exec']['RI_truncate_kde']['labels'],
    output:
        mdf = f'{RI_d}/{{sample}}/{{experiment}}/truncated.mdf',
    params:
        gtf = config['refs']['GTF'],
    shell:
        '{input.binary}'
        ' -i {input.mdf}'
        ' -o {output.mdf}'
        ' --kde={input.grid},{input.labels}'
        ' --gtf={params.gtf}'
        ' --only-single-isoform'

rule RI_splicer:
    input:
        binary = config['exec']['RI_splicer'],
        abdund_tsv = f'{preproc_d}/minimap2/{{sample}}.cDNA.abundance.tsv'
    output:
        mdf = f'{RI_d}/{{sample}}/{{experiment}}/splicer.mdf',
    params:
        gtf = config['refs']['GTF'],
        mol_count = lambda wildcards: config['experiments'][wildcards.experiment]['mol_count'],
    shell:
        '{input.binary}'
        ' -g {params.gtf}'
        ' -a {input.abdund_tsv}'
        ' --molecule-count {params.mol_count}'
        ' -o {output.mdf}'

rule compute_abundance:
    input:
        paf = f'{preproc_d}/minimap2/{{sample}}.cDNA.paf',
        script = config['exec']['ONT_trans_abund'],
    output:
        f'{preproc_d}/minimap2/{{sample}}.cDNA.abundance.tsv'
    shell:
       "python {input.script} -i {input.paf} > {output}" 

rule minimap_cdna:
    input:
        reads = lambda wildcards: config['samples'][wildcards.sample]['fastq'],
    output:
        paf = f'{preproc_d}/minimap2/{{sample}}.cDNA.paf'
    params:
        ref = config['refs']['cDNA'],
    threads:
        12
    shell:
        'minimap2 -t {threads} -x map-ont -c --eqx -p0 -o {output.paf} {params.ref} {input.reads}'

rule scTagger_match:
    input:
        lr_tsv=f'{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz',
        wl_tsv=f'{preproc_d}/scTagger/{{sample}}/{{sample}}.bc_whitelist.tsv.gz',
    output:
        lr_tsv=f'{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_matches.tsv.gz',
    params:
        script=config['exec']['scTagger'],
    threads: 32
    resources:
        mem='64G',
        time=60 * 5 - 1,
    shell:
        '{params.script} match_trie -lr {input.lr_tsv} -sr {input.wl_tsv} -o {output.lr_tsv} -t {threads}'

rule scTagger_extract_bc:
    input:
        tsv=f'{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz',
        wl=config['refs']['10x_bc'],
    output:
        tsv=f'{preproc_d}/scTagger/{{sample}}/{{sample}}.bc_whitelist.tsv.gz',
    params:
        script=config['exec']['scTagger'],
    resources:
        mem='16G',
        time=59,
    shell:
        '{params.script} extract_sr_bc_from_lr -i {input.tsv} -wl {input.wl} -o {output.tsv}'

rule scTagger_lr_seg:
    input:
        reads = lambda wildcards: config['samples'][wildcards.sample]['fastq'],
    output:
        tsv=f"{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz",
    params:
        script=config['exec']['scTagger'],
    threads:
        32
    resources:
        mem='256G',
        time=59,
    shell:
        '{params.script} extract_lr_bc -r {input.reads} -o {output.tsv} -t {threads}'