def get_strand_param(wildcards):
# if paired and and stranded then return "-sspe" else nothing.
    if REPLICATESDF["PEorSE"][wildcards.replicate] == "PE" and REPLICATESDF["stranded"][wildcards.replicate] == "Y":
        return "-sspe"
    else:
        return ""


rule homer_class:
    input:
        bam             = join(WORKDIR,"results","STAR","{replicate}","loose","{replicate}.loose.Aligned.out.bam"),
    output:
        txt             = join(WORKDIR,"results","HOMER","{replicate}","class","tagInfo.txt"),
    params:
        name            = "{replicate}",
        peorse          = get_peorse,
        strand_param    = get_strand_param,
        randomstr       = str(uuid.uuid4()),
    envmodules: TOOLS["homer"],TOOLS["samtools"],
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.txt})

makeTagDirectory \\
    $outdir \\
    {input.bam} \\
    -keepOne {params.strand_param}
"""
