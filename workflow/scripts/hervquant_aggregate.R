#!/usr/bin/env Rscript
setwd("/vf/users/EVset_RNAseq/Pipelines/ERVPipeline_runs/testRuns/hg38+viral/results/hervquant")

library("argparse")
library("xlsx")
parser <- ArgumentParser(description="Aggregate Salmon results from hervquant")

parser$add_argument("-q", "--hervquantfolder",
                    type="character",
                    help="path to hervquant folder under results output folder",
                    required=TRUE)
parser$add_argument("-l", "--lookuptable",
                    type="character",
                    help="hervquant lookup table",
                    required=TRUE)
parser$add_argument("-o", "--outfile",
                    type="character",
                    help="output TSV file",
                    required=TRUE)
parser$add_argument("-x", "--exceloutfile",
                    type="character",
                    help="output XLSX file",
                    required=TRUE)
parser$add_argument("-t", "--counttype",
		    type="character",
		    help="type of counts to aggregate on TPM or NumReads aka rawcounts",
		    default="NumReads",
		    required=FALSE)

args <- parser$parse_args()
counttype <- args$counttype
if ( counttype != "TPM" ) {
	counttype = "NumReads"
}
# args = NULL
# args$hervquantfolder = "/vf/users/EVset_RNAseq/Pipelines/ERVPipeline_runs/testRuns/hg38+viral/results/hervquant"
# args$lookuptable = "/vf/users/EVset_RNAseq/Pipelines/ERVPipeline/resources/hervquant/herv_lookup.psv"

hfiles <- sort(list.files(path=args$hervquantfolder,pattern="quant.sf",recursive=TRUE))

read_quant_file <- function(path) {
  print(path)
  pathsplit <- unlist(strsplit(path,split="/"))
  len_path <- length(pathsplit)
  rep_name <- pathsplit[len_path-1]
  contents <- read.csv(path,sep="\t",header=TRUE)
  subcontents <- contents[,c("Name",counttype)]
  colnames(subcontents) <- c("herv_id",rep_name)
  return(subcontents)
}

mergeddf=NULL
count = 0
for (p in hfiles) {
  count = count + 1
  if (count==1){
    mergeddf = read_quant_file(p)
  } else {
    readdf = read_quant_file(p)
    mergeddf = merge(mergeddf,readdf, by.x=c("herv_id"), by.y=c("herv_id"), all=TRUE)
  }
}

contents <- read.csv(args$lookuptable,sep="|",header=TRUE)
mergeddf2 <- merge(contents,mergeddf,by.x=c("herv_id"),by.y=c("herv_id"),all=TRUE)

write.table(mergeddf2,row.names=FALSE,col.names = TRUE,quote = FALSE,sep="\t",file=args$outfile)
write.xlsx2(mergeddf2,file=args$exceloutfile,row.names = FALSE)
