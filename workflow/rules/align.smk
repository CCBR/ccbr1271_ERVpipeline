rule star_loose:
    input:
        unpack(get_input_fastqs),
    output:
        bam             = join(WORKDIR,"results","STAR","{replicate}","loose","{replicate}.loose.Aligned.out.bam"),
        unmappedR1      = join(WORKDIR,"results","STAR","{replicate}","loose","{replicate}.loose.Unmapped.out.mate1.fastq.gz"),
    params:
        name            = "{replicate}",
        peorse          = get_peorse,
        starindexdir    = join(INDEXDIR,GENOME),
        randomstr       = str(uuid.uuid4()),
    threads: getthreads("star_loose")
    envmodules: TOOLS["star"],TOOLS["samtools"],TOOLS["pigz"]
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.bam})

if [[ "{params.peorse}" == "PE" ]];then
    inputfqs=" {input.R1} {input.R2} "
else
    intputfqs=" {input.R1} "
fi

TMPOUTPREFIX="${{TMPDIR}}/{params.name}.loose.unsorted."

STAR --runMode alignReads \\
    --runThreadN {threads} \\
    --genomeDir {params.starindexdir} \\
    --outFilterMultimapNmax 1000000 \\
    --readFilesIn $inputfqs \\
    --outSAMtype BAM Unsorted \\
    --outFileNamePrefix $TMPOUTPREFIX \\
    --outReadsUnmapped Fastx \\
    --readFilesCommand zcat

ls -larth $TMPDIR

samtools sort \\
    -@{threads} \\
    -m4G -l9 \\
    --output-fmt BAM \\
    --write-index \\
    -o {output.bam} \\
    -T ${{TMPDIR}}/{params.name}.samtoolstmp \\
    ${{TMPOUTPREFIX}}Aligned.out.bam

    mv ${{TMPOUTPREFIX}}Unmapped.out.mate1 ${{TMPOUTPREFIX}}Unmapped.out.mate1.fastq
    pigz -p {threads} ${{TMPOUTPREFIX}}Unmapped.out.mate1.fastq
    cp -v ${{TMPOUTPREFIX}}Unmapped.out.mate1.fastq.gz {output.unmappedR1}

if [[ "{params.peorse}" == "PE" ]];then
    unmappedR2=$(echo {output.unmappedR1} | sed "s/mate1.fastq/mate2.fastq/g")
    mv ${{TMPOUTPREFIX}}Unmapped.out.mate2 ${{TMPOUTPREFIX}}Unmapped.out.mate2.fastq
    pigz -p {threads} ${{TMPOUTPREFIX}}Unmapped.out.mate2.fastq
    cp -v ${{TMPOUTPREFIX}}Unmapped.out.mate2.fastq.gz ${{unmappedR2}}
fi

# copy other files
for suffix in "SJ.out.tab" "Log.progress.out" "Log.out" "Log.final.out"
do
cp -v ${{TMPOUTPREFIX}}${{suffix}} $outdir
done
"""

rule star_strict:
    input:
        unpack(get_input_fastqs),
    output:
        bam             = join(WORKDIR,"results","STAR","{replicate}","strict","{replicate}.strict.Aligned.out.bam"),
        unmappedR1      = join(WORKDIR,"results","STAR","{replicate}","strict","{replicate}.strict.Unmapped.out.mate1.fastq.gz"),
    params:
        name            = "{replicate}",
        peorse          = get_peorse,
        starindexdir    = join(INDEXDIR,GENOME),
        randomstr       = str(uuid.uuid4()),
    threads: getthreads("star_strict")
    envmodules: TOOLS["star"],TOOLS["samtools"]
    shell:"""
{SETSTR}
{TMPDIR_STR}
outdir=$(dirname {output.bam})

if [[ "{params.peorse}" == "PE" ]];then
    inputfqs=" {input.R1} {input.R2} "
else
    intputfqs=" {input.R1} "
fi

TMPOUTPREFIX="${{TMPDIR}}/{params.name}.strict.unsorted."

STAR --runMode alignReads \\
    --runThreadN {threads} \\
    --genomeDir {params.starindexdir} \\
    --outFilterMultimapNmax 1 \\
    --readFilesIn $inputfqs \\
    --outSAMtype BAM Unsorted \\
    --outFileNamePrefix $TMPOUTPREFIX \\
    --outReadsUnmapped Fastx \\
    --readFilesCommand zcat

ls -larth $TMPDIR

samtools sort \\
    -@{threads} \\
    -m4G -l9 \\
    --output-fmt BAM \\
    --write-index \\
    -o {output.bam} \\
    -T ${{TMPDIR}}/{params.name}.samtoolstmp \\
    ${{TMPOUTPREFIX}}Aligned.out.bam

    mv ${{TMPOUTPREFIX}}Unmapped.out.mate1 ${{TMPOUTPREFIX}}Unmapped.out.mate1.fastq
    pigz -p {threads} ${{TMPOUTPREFIX}}Unmapped.out.mate1.fastq
    cp -v ${{TMPOUTPREFIX}}Unmapped.out.mate1.fastq.gz {output.unmappedR1}

if [[ "{params.peorse}" == "PE" ]];then
    unmappedR2=$(echo {output.unmappedR1} | sed "s/mate1.fastq/mate2.fastq/g")
    mv ${{TMPOUTPREFIX}}Unmapped.out.mate2 ${{TMPOUTPREFIX}}Unmapped.out.mate2.fastq
    pigz -p {threads} ${{TMPOUTPREFIX}}Unmapped.out.mate2.fastq
    cp -v ${{TMPOUTPREFIX}}Unmapped.out.mate2.fastq.gz ${{unmappedR2}}
fi

# copy other files
for suffix in "SJ.out.tab" "Log.progress.out" "Log.out" "Log.final.out"
do
cp -v ${{TMPOUTPREFIX}}${{suffix}} $outdir
done
"""
