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

rule homer_locus:
    input:
        bam             = join(WORKDIR,"results","STAR","{replicate}","strict","{replicate}.strict.Aligned.out.bam"),
    output:
        txt             = join(WORKDIR,"results","HOMER","{replicate}","locus","tagInfo.txt"),
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
    {input.bam} {params.strand_param}
"""

rule annotate_class:
    input:
        expand(join(WORKDIR,"results","HOMER","{replicate}","class","tagInfo.txt"),replicate=REPLICATES),
    output:
        counts          = join(WORKDIR,"results","HOMER","counts","countsTable.class.tsv"),
    params:
        genome          = GENOME,
        randomstr       = str(uuid.uuid4()),
    envmodules: TOOLS["homer"],TOOLS["samtools"],
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.counts})

list_of_dirs=""
for i in {input};do
    d=$(dirname $i)
    list_of_dirs="$list_of_dirs $d"
done

analyzeRepeats.pl \\
    repeats \\
    {params.genome} \\
    -d $list_of_dirs \\
    > {output.counts}
"""

rule annotate_locus:
    input:
        expand(join(WORKDIR,"results","HOMER","{replicate}","locus","tagInfo.txt"),replicate=REPLICATES),
    output:
        counts          = join(WORKDIR,"results","HOMER","counts","countsTable.locus.tsv"),
    params:
        genome          = GENOME,
        randomstr       = str(uuid.uuid4()),
    envmodules: TOOLS["homer"],TOOLS["samtools"],
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.counts})

list_of_dirs=""
for i in {input};do
    d=$(dirname $i)
    list_of_dirs="$list_of_dirs $d"
done

analyzeRepeats.pl \\
    repeats \\
    {params.genome} \\
    -d $list_of_dirs \\
    > {output.counts}
"""

rule annotate_class_geve:
    input:
        expand(join(WORKDIR,"results","HOMER","{replicate}","class","tagInfo.txt"),replicate=REPLICATES),
    output:
        counts          = join(WORKDIR,"results","HOMER","counts","countsTable.class.geve.tsv"),
    params:
        genome          = GENOME,
        geve_gtf        = join(INDEXDIR,GENOME,GENOME+'.geve.gtf'),
        randomstr       = str(uuid.uuid4()),
    envmodules: TOOLS["homer"],TOOLS["samtools"],
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.counts})

list_of_dirs=""
for i in {input};do
    d=$(dirname $i)
    list_of_dirs="$list_of_dirs $d"
done

analyzeRepeats.pl \\
    {params.geve_gtf} \\
    {params.genome} \\
    -count genes -noadj \\
    -d $list_of_dirs \\
    > {output.counts}
"""

rule annotate_class_geve_w_annotation_table:
    input:
        counts              =   join(WORKDIR,"results","HOMER","counts","countsTable.class.geve.tsv"),
    output:
        annotated_counts    =   join(WORKDIR,"results","HOMER","counts","countsTable.class.geve.annotated.tsv"),
    params:
        rscript             =   join(SCRIPTSDIR,"annotate_counts.R"),
        genome              =   GENOME,
        geve_annotable      =   join(INDEXDIR,GENOME,GENOME+'.geve.annotation_table.tsv'),
        randomstr           =   str(uuid.uuid4()),
    envmodules: TOOLS["R"]
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.annotated_counts})
cd $outdir

Rscript {params.rscript} \\
    --count {input.counts} \\
    --anno {params.geve_annotable} \\
    --outfile {output.annotated_counts}
"""

rule annotate_locus_geve:
    input:
        expand(join(WORKDIR,"results","HOMER","{replicate}","locus","tagInfo.txt"),replicate=REPLICATES),
    output:
        counts          = join(WORKDIR,"results","HOMER","counts","countsTable.locus.geve.tsv"),
    params:
        genome          = GENOME,
        geve_gtf        = join(INDEXDIR,GENOME,GENOME+'.geve.gtf'),
        randomstr       = str(uuid.uuid4()),
    envmodules: TOOLS["homer"],TOOLS["samtools"],
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.counts})

list_of_dirs=""
for i in {input};do
    d=$(dirname $i)
    list_of_dirs="$list_of_dirs $d"
done

analyzeRepeats.pl \\
    {params.geve_gtf} \\
    {params.genome} \\
    -count genes -noadj \\
    -d $list_of_dirs \\
    > {output.counts}
"""

rule annotate_locus_geve_w_annotation_table:
    input:
        counts              =   join(WORKDIR,"results","HOMER","counts","countsTable.locus.geve.tsv"),
    output:
        annotated_counts    =   join(WORKDIR,"results","HOMER","counts","countsTable.locus.geve.annotated.tsv"),
    params:
        rscript             =   join(SCRIPTSDIR,"annotate_counts.R"),
        genome              =   GENOME,
        geve_annotable      =   join(INDEXDIR,GENOME,GENOME+'.geve.annotation_table.tsv'),
        randomstr           =   str(uuid.uuid4()),
    envmodules: TOOLS["R"]
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.annotated_counts})
cd $outdir

Rscript {params.rscript} \\
    --count {input.counts} \\
    --anno {params.geve_annotable} \\
    --outfile {output.annotated_counts}
"""
