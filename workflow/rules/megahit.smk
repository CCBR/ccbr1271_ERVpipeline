rule megahit:
    input:
        unpack(get_input_fastqs),
    output:
        fa          = join(WORKDIR,"results","megahit","{replicate}","{replicate}.contigs.fa")
    params:
        name        = "{replicate}",
        peorse      = get_peorse,
        randomstr   = str(uuid.uuid4()),
    threads: getthreads("megahit")
    envmodules: TOOLS["megahit"]
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.fa})

if [[ -d $outdir ]];then rm -rf $outdir;fi

if [[ "{params.peorse}" == "PE" ]];then
    megahit \\
        -t{threads} \\
        -1 {input.R1} \\
        -2 {input.R2} \\
        --out-dir $outdir \\
        --out-prefix {params.name} \\
        --tmp-dir $TMPDIR
else
    megahit \\
        -t{threads} \\
        --read {input.R1} \\
        --out-dir $outdir \\
        --out-prefix {params.name} \\
        --tmp-dir $TMPDIR
fi

# recover 80+ percent of disk space
if [[ -d "${outdir}/intermediate_contigs" ]]; then rm -rf ${outdir}/intermediate_contigs;fi
"""
