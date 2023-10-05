## Endogenous Retro Virus (ERV) detection and quantification pipeline

Straight forward usage:

```bash
╰─⠠⠵ ./erv

##########################################################################################

Welcome to

╭━━━╮╱╱╱╱╭╮╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╭━━━╮╱╱╭╮╱╱╱╱╱╱╭╮╱╱╭╮╱╱╱╱╱╱╱╱╭━━━╮╱╱╱╱╱╱╭╮
┃╭━━╯╱╱╱╱┃┃╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱┃╭━╮┃╱╭╯╰╮╱╱╱╱╱┃╰╮╭╯┃╱╱╱╱╱╱╱╱┃╭━╮┃╱╱╱╱╱╱┃┃
┃╰━━┳━╮╭━╯┣━━┳━━┳━━┳━╮╭━━┳━━┳╮╭┳━━╮┃╰━╯┣━┻╮╭╋━┳━━╮╰╮┃┃╭╋━┳╮╭┳━━╮┃╰━╯┣┳━━┳━━┫┃╭┳━╮╭━━╮
┃╭━━┫╭╮┫╭╮┃╭╮┃╭╮┃┃━┫╭╮┫┃━┫╭╮┃┃┃┃━━┫┃╭╮╭┫┃━┫┃┃╭┫╭╮┃╱┃╰╯┣┫╭┫┃┃┃━━┫┃╭━━╋┫╭╮┃┃━┫┃┣┫╭╮┫┃━┫
┃╰━━┫┃┃┃╰╯┃╰╯┃╰╯┃┃━┫┃┃┃┃━┫╰╯┃╰╯┣━━┃┃┃┃╰┫┃━┫╰┫┃┃╰╯┃╱╰╮╭┫┃┃┃╰╯┣━━┃┃┃╱╱┃┃╰╯┃┃━┫╰┫┃┃┃┃┃━┫
╰━━━┻╯╰┻━━┻━━┻━╮┣━━┻╯╰┻━━┻━━┻━━┻━━╯╰╯╰━┻━━┻━┻╯╰━━╯╱╱╰╯╰┻╯╰━━┻━━╯╰╯╱╱╰┫╭━┻━━┻━┻┻╯╰┻━━╯
╱╱╱╱╱╱╱╱╱╱╱╱╱╭━╯┃╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱┃┃
╱╱╱╱╱╱╱╱╱╱╱╱╱╰━━╯╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╰╯
... v0.0.1

##########################################################################################

This pipeline was built by CCBR (https://bioinformatics.ccr.cancer.gov/ccbr)
Please contact Vishal Koparde for comments/questions (vishal.koparde@nih.gov)

##########################################################################################

Here is a list of genome supported by this pipeline:

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
    * reset     : DELETE workdir dir and re-init it (debugging option) EDITS TO ALL FILES IN WORKDIR WILL BE LOST!
    * local     : same as runlocal

Optional Arguments:

--genome|-g     : genome eg. hg38(default) or mm10
--manifest|-s   : absolute path to samples.tsv. This will be copied to output folder  (--runmode=init only)
--help|-h       : print this help

Example commands:
  ./erv -w=/my/ouput/folder -m=init [ -g="mm10" -s="/path/to/sample.tsv" ]
  ./erv -w=/my/ouput/folder -m=dryrun
  ./erv -w=/my/ouput/folder -m=run

##########################################################################################

VersionInfo:
  python          : 3.10
  snakemake       : 7.32.3
  pipeline_home   : /vf/users/EVset_RNAseq/Pipelines/ERVPipeline/dev
  git commit/tag  : f5f99df9a56a60f7a12f83658291b0f0e533a1f7
  pipeline_version: v0.0.1

##########################################################################################
```

Please reach out to [Vishal Koparde](mailto:vishal.koparde@nih.gov) for comments/enhancements/etc. You can also [submit an issue](https://github.com/CCBR/ccbr1271_ERVpipeline/issues/new/choose).