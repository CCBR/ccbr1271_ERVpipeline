def get_hervquant_input(wildcards):
    """
    Return a dictionary of input files for hervquant rule
    """
    d=dict()
    if GENOME.startswith("hg"):
        d["R1"]=REPLICATESDF["R1"][wildcards.replicate]
        d["R2"]=REPLICATESDF["R2"][wildcards.replicate]
    return d

rule hervquant:
    input:
        unpack(get_hervquant_input),
    output:
        dummy                       = join(WORKDIR,"results","hervquant","{replicate}","{replicate}.dummy"),
        filtered_bam                = join(WORKDIR,"results","hervquant","{replicate}","{replicate}.Aligned.out.filtered.bam"),
    params:
        name                        = "{replicate}",
        outdir                      = join(WORKDIR,"results","hervquant","{replicate}"),
        peorse                      = get_peorse,
        randomstr                   = str(uuid.uuid4()),
        hervquantindexdir           = join(INDEXDIR,"hervquant"),
        hervquant_final_reference   = join(INDEXDIR,"hervquant","hervquant_final_reference.fa"),
    threads: getthreads("hervquant")
    envmodules: TOOLS["star"],
                TOOLS["samtools"],
                TOOLS["salmon"]
    shell:"""
{SETSTR}
{TMPDIR_STR}

touch {output.dummy}

# file_prefix="{params.outdir}/"
file_prefix="${{TMPDIR}}/{params.name}."

# align reads to reference
STAR \
 --runThreadN {threads} \
 --outFileNamePrefix $file_prefix \
 --outFilterMultimapNmax 10 \
 --outFilterMismatchNmax 7 \
 --genomeDir {params.hervquantindexdir} \
 --readFilesIn {input.R1} {input.R2} \
 --readFilesCommand zcat \
 --outTmpDir ${{TMPDIR}}/STARtmp

cp ${{file_prefix}}Log.* {params.outdir}

#Filter out all non-herv maps:
sam_file=${{file_prefix}}Aligned.out.sam
sam_file_filtered=${{file_prefix}}Aligned.out.filtered.sam
sed '/uc.*/d' $sam_file > $sam_file_filtered
samtools view -bS $sam_file_filtered > {output.filtered_bam}

# assemble reads
salmon quant \
 -t {params.hervquant_final_reference} \
 -l ISF \
 -a {output.filtered_bam} \
 -o {params.outdir} \
 -p {threads}
"""
