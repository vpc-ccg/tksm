import sys

if len(config)==0:
    configfile: 'config.yaml'

outpath = config['outpath']
preproc_d = f'{outpath}/preprocess'
RI_d = f'{outpath}/RI'
DEBUG=True

if DEBUG:
    def pipe(X):
        return X
    def temp(X):
        return X

def exprmnt_sample(exprmnt):
    return config['RI_experiments'][exprmnt]['sample']

def component_idx(prefix):
    return len(prefix.split('.'))

def fastas_for_RI_sequence(wc):
    fastas = list()
    fastas.append(config['refs']['DNA'])
    for idx, component in enumerate(config['RI_experiments'][wc.exprmnt]['pipeline']):
        component = list(component.keys())[0]
        if component in ['plA', 'SCB', 'UMI']:
            prefix = '.'.join(
                [
                    list(c.keys())[0]
                    for c in config['RI_experiments'][wc.exprmnt]['pipeline'][:idx+1]
                ]
            )
            fastas.append(
                f'{RI_d}/{wc.exprmnt}/{prefix}.fasta',
            )
    return fastas

def experiment_prefix(exprmnt):
    prefix = list()
    for component in config['RI_experiments'][exprmnt]['pipeline']:
        prefix.append(list(component)[0])
    return '.'.join(prefix)
    
rule all:
    input:
        [
            f'{RI_d}/{exprmnt}/{experiment_prefix(exprmnt)}.fastq'
            for exprmnt in config['RI_experiments']
        ],
       

rule make_binary:
    input:
        "src/{prefix}.cpp"
    output:
        "build/{prefix}"
    shell:
        "make {output}"

rule sequence:
    input:
        binary = config['exec']['RI_sequencer'],
        badread = config['exec']['badread'],
        mdf = f'{RI_d}/{{exprmnt}}/{{prefix}}.mdf',
        fastas = fastas_for_RI_sequence,
    output:
        fastq = f'{RI_d}/{{exprmnt}}/{{prefix}}.Seq.fastq',
    params:
        name = lambda wc: f'{exprmnt_sample(wc.exprmnt)}_{wc.exprmnt}_RI',
        tmp_dir = f'{RI_d}/{{exprmnt}}/{{prefix}}.Seq.tmps',
        fastas = lambda wc: ','.join(fastas_for_RI_sequence(wc)),
        other = lambda wc: config['RI_experiments'][wc.exprmnt]['pipeline'][component_idx(wc.prefix)]['Seq'],
    threads:
        32
    shell:
        'rm -rf {params.tmp_dir} && '
        'mkdir -p {params.tmp_dir} && '
        '{input.binary}'
        ' -i {input.mdf}'
        ' --references={params.fastas}'
        ' -o {output.fastq}'
        ' --temp "{params.tmp_dir}"'
        ' -n "{params.name}"'
        ' --badread={input.badread}'
        ' --threads={threads}'
        ' {params.other}'

rule trunc:
    input:
        binary = config['exec']['RI_truncate'],
        mdf = f'{RI_d}/{{exprmnt}}/{{prefix}}.mdf',
        x = lambda wc: f'{preproc_d}/truncate_kde/{exprmnt_sample(wc.exprmnt)}.X_idxs.npy',
        y = lambda wc: f'{preproc_d}/truncate_kde/{exprmnt_sample(wc.exprmnt)}.Y_idxs.npy',
        g = lambda wc: f'{preproc_d}/truncate_kde/{exprmnt_sample(wc.exprmnt)}.grid.npy',
        gtf = config['refs']['GTF'],
    output:
        mdf = f'{RI_d}/{{exprmnt}}/{{prefix}}.Trc.mdf',
    params:
        lambda wc: config['RI_experiments'][wc.exprmnt]['pipeline'][component_idx(wc.prefix)]['Trc'],
    shell:
        '{input.binary}'
        ' -i {input.mdf}'
        ' --kde={input.g},{input.x},{input.y}'
        ' --gtf={input.gtf}'
        ' -o {output.mdf}'
        ' {params}'

rule pcr:
    input:
        binary = config['exec']['RI_pcr'],
        mdf = f'{RI_d}/{{exprmnt}}/{{prefix}}.mdf',
    output:
        mdf = f'{RI_d}/{{exprmnt}}/{{prefix}}.PCR.mdf',
    params:
        lambda wc: config['RI_experiments'][wc.exprmnt]['pipeline'][component_idx(wc.prefix)]['PCR'],
    shell:
        '{input.binary}'
        ' -i {input.mdf}'
        ' -o {output.mdf}'
        ' {params}'

rule umi:
    input:
        binary = config['exec']['RI_umi'],
        mdf = f'{RI_d}/{{exprmnt}}/{{prefix}}.mdf',
    output:
        mdf = f'{RI_d}/{{exprmnt}}/{{prefix}}.UMI.mdf',
        fasta = f'{RI_d}/{{exprmnt}}/{{prefix}}.UMI.fasta',
    params:
        lambda wc: config['RI_experiments'][wc.exprmnt]['pipeline'][component_idx(wc.prefix)]['UMI'],
    shell:
        '{input.binary}'
        ' -i {input.mdf}'
        ' -o {output.mdf}'
        ' -f {output.fasta}'
        ' {params}'

rule single_cell:
    input:
        binary = config['exec']['RI_sc_barcoder'],
        mdf = f'{RI_d}/{{exprmnt}}/{{prefix}}.mdf',
    output:
        mdf = f'{RI_d}/{{exprmnt}}/{{prefix}}.SCB.mdf',
        fasta = f'{RI_d}/{{exprmnt}}/{{prefix}}.SCB.fasta',
    params:
        lambda wc: config['RI_experiments'][wc.exprmnt]['pipeline'][component_idx(wc.prefix)]['SCB'],
    shell:
        '{input.binary}'
        ' -i {input.mdf}'
        ' -o {output.mdf}'
        ' -f {output.fasta}'
        ' {params}'

rule polyA:
    input:
        binary = config['exec']['RI_polyA'],
        mdf = f'{RI_d}/{{exprmnt}}/{{prefix}}.mdf',
    output:
        mdf = f'{RI_d}/{{exprmnt}}/{{prefix}}.plA.mdf',
        fasta = f'{RI_d}/{{exprmnt}}/{{prefix}}.plA.fasta',
    params:
        lambda wc: config['RI_experiments'][wc.exprmnt]['pipeline'][component_idx(wc.prefix)]['plA'],
    shell:
        '{input.binary}'
        ' -i {input.mdf}'
        ' -o {output.mdf}'
        ' -f {output.fasta}'
        ' {params}'

rule splicer:
    input:
        binary = config['exec']['RI_splicer'],
        tsv = f'{RI_d}/{{exprmnt}}/{{prefix}}.tsv',
        gtf = config['refs']['GTF'],
    output:
        mdf = f'{RI_d}/{{exprmnt}}/{{prefix}}.Spc.mdf',
    params:
        lambda wc: config['RI_experiments'][wc.exprmnt]['pipeline'][component_idx(wc.prefix)]['Spc'],
    shell:
        '{input.binary}'
        ' -a {input.tsv}'
        ' -g {input.gtf}'
        ' -o {output.mdf}'
        ' {params}'

rule expression:
    input:
        script = config['exec']['transcript_abundance'],
        paf = lambda wc: f'{preproc_d}/minimap2/{exprmnt_sample(wc.exprmnt)}.cDNA.paf',
    output:
        tsv = f'{RI_d}/{{exprmnt}}/Xpr.tsv',
    shell:
       'python {input.script} -p {input.paf} -o {output.tsv}' 

rule expression_sc:
    input:
        script = config['exec']['transcript_abundance'],
        paf = lambda wc: f'{preproc_d}/minimap2/{exprmnt_sample(wc.exprmnt)}.cDNA.paf',
        lr_matches = lambda wc: f'{preproc_d}/scTagger/{exprmnt_sample(wc.exprmnt)}/{exprmnt_sample(wc.exprmnt)}.lr_matches.tsv.gz'
    output:
        tsv = f'{RI_d}/{{exprmnt}}/Xpr_sc.tsv',
    shell:
       'python {input.script} -p {input.paf} -m {input.lr_matches} -o {output.tsv}' 

rule truncate_kde:
    input:
        script = config['exec']['truncate_kde'],
        paf = f'{preproc_d}/minimap2/{{sample}}.cDNA.paf'
    output:
        x = f'{preproc_d}/truncate_kde/{{sample}}.X_idxs.npy',
        y = f'{preproc_d}/truncate_kde/{{sample}}.Y_idxs.npy',
        g = f'{preproc_d}/truncate_kde/{{sample}}.grid.npy',
    params:
        out_prefix = f'{preproc_d}/truncate_kde/{{sample}}',
    threads:
        32
    shell:
       'python {input.script} -i {input.paf} -o {params.out_prefix} --threads {threads}' 

rule minimap_cdna:
    input:
        reads = lambda wildcards: config['samples'][wildcards.sample]['fastq'],
        ref = config['refs']['cDNA'],
    output:
        paf = f'{preproc_d}/minimap2/{{sample}}.cDNA.paf'
    threads:
        32
    shell:
        'minimap2 -t {threads} -x map-ont -c --eqx -p0 -o {output.paf} {input.ref} {input.reads}'

rule scTagger_match:
    input:
        lr_tsv=f'{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz',
        wl_tsv=f'{preproc_d}/scTagger/{{sample}}/{{sample}}.bc_whitelist.tsv.gz',
    output:
        lr_tsv=f'{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_matches.tsv.gz',
    threads: 32
    resources:
        mem='64G',
        time=60 * 5 - 1,
    shell:
        'scTagger.py match_trie -lr {input.lr_tsv} -sr {input.wl_tsv} -o {output.lr_tsv} -t {threads}'

rule scTagger_extract_bc:
    input:
        tsv=f'{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz',
        wl=config['refs']['10x_bc'],
    output:
        tsv=f'{preproc_d}/scTagger/{{sample}}/{{sample}}.bc_whitelist.tsv.gz',
    resources:
        mem='16G',
        time=59,
    shell:
        'scTagger.py extract_sr_bc_from_lr -i {input.tsv} -wl {input.wl} -o {output.tsv}'

rule scTagger_lr_seg:
    input:
        reads = lambda wildcards: config['samples'][wildcards.sample]['fastq'],
    output:
        tsv=f'{preproc_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz',
    threads:
        32
    resources:
        mem='256G',
        time=59,
    shell:
        'scTagger.py extract_lr_bc -r {input.reads} -o {output.tsv} -t {threads}'
