ERV identification, annotation and quantification pipeline built for ccbr1271. This pipeline detects and quantifies Endogenous Retroviruses using the scripts obtained from the Belkaid group.

- This pipeline built using [Snakemake](https://snakemake.readthedocs.io/en/stable/).
- This pipeline uses modules preinstalled and available on [BIOWULF](https://hpc.nih.gov/).
- This pipeline takes adapter removed and preprocess FASTQs from [ccbr1271_longRNA](https://github.com/CCBR/ccbr1271_longRNA) pipeline as inputs.
- This pipeline has 3 distinct branches:
  - Alignment with [STAR](https://github.com/alexdobin/STAR) followed by repeat elements quantification using [HOMER](http://homer.ucsd.edu/homer/)
  - Assembly of reads using [MEGAHIT](https://github.com/voutcn/megahit) followed by contig alignments with NR database using [DIAMOND](https://github.com/bbuchfink/diamond)
  - Quantification + annotation of human ERVs using the [hervQuant](https://unclineberger.org/vincentlab/resources/) pipeline
