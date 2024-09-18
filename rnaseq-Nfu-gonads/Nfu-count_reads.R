## This script does the read counting from mapping results (bam files) and ends with a csv file of raw counts.
## It is NOT intended to be run fully because the bam files are not provided along with this data package.

# Original Author: Andreas Petzold, rearranged by Philipp Koch
# Date: 2014, 2024

library(AnnotationDbi)
library(GenomicFeatures)
library(GenomicAlignments)
library(edgeR)

data_path="."

############################################
# 0. create gene and metadata tables #
############################################
# read existing transcript database
# NOTE: The sqlite file is gzip compressed. It needs to be uncompressed prior to loadDb.
# If the system command does not work, the decompression needs to be done manually.
system("gunzip --force --keep nfurzeri_genebuild_v1.140704.sqlite.gz")
nfu_db = loadDb('nfurzeri_genebuild_v1.140704.sqlite')
nfu_exonsByGene = exonsBy(nfu_db, by="gene")

# get gene lengths
gene_lengths = sapply(as.list(width(nfu_exonsByGene)),sum)


###########################
# 1. create data sets     #
###########################

# meta and counts from bam
grz_runs.meta = read.table(paste(data_path,"GRZ_GONADS.meta.csv",sep="/"),sep="\t",header=TRUE)
grz_runs.meta$file = paste("/misc/vulpix/data/NFU_ANNO_20130926/rnaseq/",grz_runs.meta$run,".bam",sep="")
grz_runs.meta$age = factor(grz_runs.meta$age,levels=c("10dpf","0d","3d","3m"))
grz_runs.meta$sex = factor(grz_runs.meta$sex,levels=c("male","female"))
grz_runs.meta$run = factor(grz_runs.meta$run,levels=unique(grz_runs.meta$run))
grz_runs.meta$tissue = factor(grz_runs.meta$tissue,levels=c("embryo","trunk","testis","ovary"))

# Note: metadata must contain a column with 'file'
create_summarised_experiment_from_bam = function(exonsByGene,metadata,mode="Union",singleEnd=FALSE,ignore.strand=TRUE,fragments=TRUE){
  sampleInfo = DataFrame(metadata)
  bamLst = BamFileList(metadata$file,yieldSize=3000000)
  summarised_experiment = summarizeOverlaps(exonsByGene,bamLst,mode=mode,singleEnd=singleEnd,ignore.strand=ignore.strand,fragments=fragments)
  summarised_experimentIdx = match(colnames(summarised_experiment), basename(sampleInfo$file))
  colData(summarised_experiment) = cbind(colData(summarised_experiment),sampleInfo[summarised_experimentIdx, ])
  summarised_experiment
}

if(!exists("summarised_experiment")){
  summarised_experiment = create_summarised_experiment_from_bam(nfu_exonsByGene,grz_runs.meta)
}

# this is not saved as the subsequent script generates this object from the counts table.
#save(summarised_experiment, file = "summarised_experiment.RData")


## preparing counts

# raw counts
counts = assay(summarised_experiment)
colnames(counts) = gsub('\\.bam','',colnames(counts))
write.csv(counts, file = "counts.csv")
