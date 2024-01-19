#!/usr/bin/env Rscript
library("argparse")
library("tidyverse")

get_sample_name <- function(path) {
  x <- unlist(strsplit(path,split="/"))
  if (length(x) > 2) {
    return(x[length(x)-1])
  } else {
    return(path)
  }
}

parser <- ArgumentParser()
parser$add_argument("-c", "--count", default="countsTable.class.geve.tsv")
parser$add_argument("-a", "--anno", default="mm10.geve.annotation_table.tsv")
parser$add_argument("-o", "--outfile", default="countsTable.class.geve.annotated.tsv")
args <- parser$parse_args()

# setwd("/data/EVset_RNAseq/Pipelines/ERVPipeline/mm10+viral/test1/results/HOMER/counts")
anno=read.csv(args$anno,sep = "\t",header=TRUE,
              check.names = FALSE,
              strip.white = TRUE,
              blank.lines.skip = TRUE)
# View(anno)
anno %>% select(!c("chr","start","end","strand")) -> anno
counts=read.table(args$count, sep = "\t",header = TRUE,
                  check.names = FALSE,
                  strip.white = TRUE,
                blank.lines.skip = TRUE)
# View(counts)
cn=colnames(counts)
cn[[1]]="ID"
# cn <- lapply(get_sample_name,cn)
colnames(counts)=cn

dim(counts)
dim(anno)
counts_annotated = merge(x=counts, y=anno, by="ID", all.x=TRUE)
# View(counts_annotated)
counts_annotated[is.na(counts_annotated)]="NA"

write.table(counts_annotated,file=args$outfile, sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
