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

rule annotate_diamond:
    input:
        dmnd_in     = join(WORKDIR,"results","diamond","{replicate}","{replicate}.nr.tsv"),
    output:
        dmnd_out    = temp(join(WORKDIR,"results","diamond","{replicate}","{replicate}.nr.annotated.tsv.tmp"))
    params:
        name        = "{replicate}",
        randomstr   = str(uuid.uuid4()),
        script      = join(SCRIPTSDIR,"annotate_diamond.py"),
        accid2taxid = config["accid2taxid"],
        taxid2lineage = config["taxid2lineage"],
        nrtitles    = config["nrtitles"],
    threads: getthreads("annotate_diamond")
    envmodules: TOOLS["python"]
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.dmnd_out})

python {params.script} \\
	-d {input.dmnd_in} \\
	-p {params.accid2taxid} \\
	-l {params.taxid2lineage} \\
	-t {params.nrtitles} \\
	-o {output.dmnd_out}
"""

rule reannotate_diamond:
    input:
        dmnd_in     = join(WORKDIR,"results","diamond","{replicate}","{replicate}.nr.annotated.tsv.tmp"),
    output:
        dmnd_out    = join(WORKDIR,"results","diamond","{replicate}","{replicate}.nr.annotated.tsv")
    params:
        name        = "{replicate}",
        randomstr   = str(uuid.uuid4()),
        script      = join(SCRIPTSDIR,"reprocess_annotated_diamond.py"),
        taxid2lineage = config["taxid2lineage"],
    threads: getthreads("reannotate_diamond")
    envmodules: TOOLS["python"]
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.dmnd_out})

python {params.script} \\
        -a {input.dmnd_in} \\
        -l {params.taxid2lineage} \\
        -o {output.dmnd_out}
"""
