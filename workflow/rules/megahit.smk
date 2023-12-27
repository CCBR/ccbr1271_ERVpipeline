rule megahit:
    input:
        unpack(get_input_fastqs),
    output:
        fa          = join(WORKDIR,"results","megahit","{replicate}","{replicate}.contigs.fa")
    params:
        name        = "{replicate}",
        outdir      = join(WORKDIR,"results","megahit","{replicate}"),
        peorse      = get_peorse,
        randomstr   = str(uuid.uuid4()),
    threads: getthreads("megahit")
    envmodules: TOOLS["megahit"]
    shell:"""
{SETSTR}
{TMPDIR_STR}

if [[ -d {params.outdir} ]];then rm -rf {params.outdir};fi

if [[ "{params.peorse}" == "PE" ]];then
    megahit \\
        -t{threads} \\
        -1 {input.R1} \\
        -2 {input.R2} \\
        --out-dir {params.outdir} \\
        --out-prefix {params.name} \\
        --tmp-dir $TMPDIR
else
    megahit \\
        -t{threads} \\
        --read {input.R1} \\
        --out-dir {params.outdir} \\
        --out-prefix {params.name} \\
        --tmp-dir $TMPDIR
fi

# recover 80+ percent of disk space
if [[ -d "{params.outdir}/intermediate_contigs" ]]; then rm -rf {params.outdir}/intermediate_contigs;fi
"""

rule diamond_blastx:
    input:
        fa          = join(WORKDIR,"results","megahit","{replicate}","{replicate}.contigs.fa"),
    output:
        dmnd_out    = join(WORKDIR,"results","diamond","{replicate}","{replicate}.nr.tsv")
    params:
        name        = "{replicate}",
        randomstr   = str(uuid.uuid4()),
        dmnd_nr     = join(INDEXDIR,"nr","nr.dmnd"),
    threads: getthreads("diamond_blastx")
    envmodules: TOOLS["diamond"]
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.dmnd_out})

diamond \\
    blastx \\
    -d {params.dmnd_nr} \\
    -q {input.fa} \\
    -p {threads} \\
    --outfmt 6 \\
    -o {output.dmnd_out} \\
    -k 2 \\
    --strand both
"""

rule diamond_blastx_get_metadata:
    input:
        dmnd_in     = join(WORKDIR,"results","diamond","{replicate}","{replicate}.nr.tsv"),
    output:
        dmnd_out    = join(WORKDIR,"results","diamond","{replicate}","{replicate}.nr.taxid.tsv")
    params:
        name        = "{replicate}",
        randomstr   = str(uuid.uuid4()),
        scriptsdir  = SCRIPTSDIR,
        script      = "get_taxid.py",
    threads: getthreads("diamond_blastx")
    envmodules: TOOLS["diamond"]
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.dmnd_out})

cat {input.dmnd_in} | \\
    python {params.scriptsdir}/{params.script} \\
    > {output.dmnd_out}
"""
