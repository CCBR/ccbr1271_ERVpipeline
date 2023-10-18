# ERV identification, annotation and quantification pipeline built for ccbr1271

### Table of Contents

- [Table of Contents](#table-of-contents)
  - [1. Disclaimers](#1-disclaimers)
  - [2. Background](#2-background)
  - [3. Flowchart](#3-flowchart)
  - [4. Technical Details](#4-technical-details)
    - [4.1 Usage](#41-usage)
    - [4.2 Resources](#42-resources)
    - [4.3 Location](#43-location)
  - [5. Version Notes](https://github.com/CCBR/RENEE/blob/main/CHANGELOG.md)

### 1. Disclaimers

> DISCLAIMERS:
>
> - Built specifically for [ccbr1271](https://abcs-amp.nih.gov/project/2487/view/) analysis
> - Tested only on [BIOWULF](https://hpc.nih.gov/)
> - Uses BIOWULF [modules](https://hpc.nih.gov/apps/modules.html)

<hr>
<p align="center">
	<a href="#erv-identification-annotation-and-quantification-pipeline-built-for-ccbr1271">Back to Top</a>
</p>
<hr>

### 2. Background

This pipeline detects and quantifies Endogenous Retroviruses using the scripts obtained from the Belkaid group.

<hr>
<p align="center">
	<a href="#erv-identification-annotation-and-quantification-pipeline-built-for-ccbr1271">Back to Top</a>
</p>
<hr>

### 3. Flowchart

![Flowchart](./docs/assets/images/ccbr1271_ERV_pipeline.png)

<hr>
<p align="center">
	<a href="#erv-identification-annotation-and-quantification-pipeline-built-for-ccbr1271">Back to Top</a>
</p>
<hr>

### 4. Technical Details

#### 4.1 Usage

```bash
⠠⠵ ./erv

##########################################################################################

Welcome to

╭━━━╮╱╱╱╱╭╮╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╭━━━╮╱╱╭╮╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╭━━━╮╱╱╱╱╱╱╭╮
┃╭━━╯╱╱╱╱┃┃╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱┃╭━╮┃╱╭╯╰╮╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱┃╭━╮┃╱╱╱╱╱╱┃┃
┃╰━━┳━╮╭━╯┣━━┳━━┳━━┳━╮╭━━┳━━┳╮╭┳━━╮┃╰━╯┣━┻╮╭╋━┳━━┳╮╭┳┳━┳╮╭┳━━╮┃╰━╯┣┳━━┳━━┫┃╭┳━╮╭━━╮
┃╭━━┫╭╮┫╭╮┃╭╮┃╭╮┃┃━┫╭╮┫┃━┫╭╮┃┃┃┃━━┫┃╭╮╭┫┃━┫┃┃╭┫╭╮┃╰╯┣┫╭┫┃┃┃━━┫┃╭━━╋┫╭╮┃┃━┫┃┣┫╭╮┫┃━┫
┃╰━━┫┃┃┃╰╯┃╰╯┃╰╯┃┃━┫┃┃┃┃━┫╰╯┃╰╯┣━━┃┃┃┃╰┫┃━┫╰┫┃┃╰╯┣╮╭┫┃┃┃╰╯┣━━┃┃┃╱╱┃┃╰╯┃┃━┫╰┫┃┃┃┃┃━┫
╰━━━┻╯╰┻━━┻━━┻━╮┣━━┻╯╰┻━━┻━━┻━━┻━━╯╰╯╰━┻━━┻━┻╯╰━━╯╰╯╰┻╯╰━━┻━━╯╰╯╱╱╰┫╭━┻━━┻━┻┻╯╰┻━━╯
╱╱╱╱╱╱╱╱╱╱╱╱╱╭━╯┃╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱┃┃
╱╱╱╱╱╱╱╱╱╱╱╱╱╰━━╯╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╰╯
... v0.0.0-dev

  * hg38          [Human]
  * mm10          [Mouse]

USAGE:
  /path/to/erv -w/--workdir=<WORKDIR> -m/--runmode=<RUNMODE>

Required Arguments:
1.  WORKDIR     : [Type: String]: Absolute or relative path to the output folder with write permissions.

2.  RUNMODE     : [Type: String] Valid options:
    * init      : initialize workdir
    * dryrun    : dry run snakemake to generate DAG
    * run       : run with slurm
    * runlocal  : run without submitting to sbatch
    ADVANCED RUNMODES (use with caution!!)
    * unlock    : unlock WORKDIR if locked by snakemake NEVER UNLOCK WORKDIR WHERE PIPELINE IS CURRENTLY RUNNING!
    * reconfig  : recreate config file in WORKDIR (debugging option) EDITS TO config.yaml WILL BE LOST!
    * recopy    : recreate tools.yaml, cluster.yaml and scriptsdir in WORKDIR (debugging option) EDITS TO these files WILL BE LOST!
    * reset     : DELETE workdir dir and re-init it (debugging option) EDITS TO ALL FILES IN WORKDIR WILL BE LOST!
    * local     : same as runlocal

Optional Arguments:

--genome|-g     : genome eg. hg38(default) or mm10
--manifest|-s   : absolute path to samples.tsv. This will be copied to output folder  (--runmode=init only)
--help|-h       : print this help

Example commands:
  ./erv -w=/my/output/folder -m=init [ -g="mm10" -s="/path/to/sample.tsv" ]
  ./erv -w=/my/output/folder -m=dryrun
  ./erv -w=/my/output/folder -m=run

##########################################################################################

VersionInfo:
  python          : 3.10
  snakemake       : 7.32.3
  pipeline_home   : /vf/users/EVset_RNAseq/Pipelines/ERVPipeline/dev
  git commit/tag  : 16c269535311bf9e80350a2de83c9c8273652acd
  pipeline_version: v0.0.0-dev

##########################################################################################
```

<hr>
<p align="center">
	<a href="#erv-identification-annotation-and-quantification-pipeline-built-for-ccbr1271">Back to Top</a>
</p>
<hr>

#### 4.2 Resources

- different versions of the pipeline are retained at `/data/EVset_RNAseq/Pipelines/ERVPipeline`
- mm10 (mouse) and hg38 (human) genomes were downloaded from UCSC directly. Only canonical chromosomes were retained, that is, chr1 through chr19 plus chrX, chrY and chrM for mouse; and chr1 through chr22 plus chrX, chrY and chrM for human.

```bash
⠠⠵ wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
⠠⠵ wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```

- STAR indices were created using version 2.7.9a and located at `/data/EVset_RNAseq/Pipelines/ERVPipeline/resources`.
- geve GTFs were provided by Stephanie. These seemed to be modified versions of files downloaded from [GEVE](http://geve.med.u-tokai.ac.jp/)
- DIAMOND nr database was creating using the following commands:

```bash
## NR fasta
# Download the entire NR fasta file
⠠⠵ wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

## build diamond index
# ref: [diamond wiki](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options)

### taxonmap file
⠠⠵ wget //ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz

### taxonnodes file
# `nodes.dmp` is part of this zip
⠠⠵ wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip

### taxonnames file
# taxonnodes zip file contains the `names.dmp`

### command

⠠⠵ diamond makedb \
 --threads 8 \
 --in nr.gz \
 --db nr.diamond_index \
 --taxonmap prot.accession2taxid.FULL.gz \
 --taxonnodes nodes.dmp \
 --taxonnames names.dmp
```

#### 4.3 Location

This pipeline is located at `/data/EVset_RNAseq/Pipelines/ERVPipeline` on [BIOWULF](https://hpc.nih.gov).

#### 4.4 Inputs

Samples are assumed to be single-end or paired-end fastqs which have already been adapter and UMI-trimmed. Intermediate fastqs from the [longRNA pipeline](https://github.com/CCBR/ccbr1271_longRNA) can be used as inputs. The sample manifest is expected to be tab-delimited with these columns headers:

- sampleName
- path_to_R1_fastq
- path_to_R2_fastq

<hr>
<p align="center">
	<a href="#erv-identification-annotation-and-quantification-pipeline-built-for-ccbr1271">Back to Top</a>
</p>
<hr>

> Please reach out to [Vishal Koparde, Ph.D.](mailto:vishal.koparde@nih.gov) from [CCBR](https://bioinformatics.ccr.cancer.gov/ccbr) for comments/questions/requests/etc.
