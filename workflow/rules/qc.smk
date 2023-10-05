def get_all_fastqs(wildcards):
    return ALLINPUTFASTQS

rule fastqc:
# """
# Run FASTQC on:
# * Raw fastqs
# * Blacklist filtered trimmed fastqs
# """
    input:
        unpack(get_all_fastqs),
    output:
        expand(join(QCDIR,"fastqc","{replicate}.R1_fastqc.zip"), replicate=REPLICATES),
    params:
        outdir=(join(QCDIR,"fastqc")),
        randomstr=str(uuid.uuid4()),
    threads: getthreads("fastqc")
    envmodules: TOOLS['fastqc']
    shell: """
{SETSTR}
{TMPDIR_STR}
fastqc {input} -t {threads} -o {params.outdir};
"""      

#########################################################