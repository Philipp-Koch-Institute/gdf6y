## This script was originally run under R 3.1.3
## For this R version, no renv package was available.

# Original Author: Andreas Petzold, rearranged by Philipp Koch
# Date: 2014, 2024

library(DESeq2)
library(AnnotationDbi)

data_path="."

############################################
# 0. create gene and metadata tables #
############################################
# read existing transcript database
# NOTE: The sqlite file is gzip compressed. It needs to be uncompressed prior to loadDb.
# If the following system command does not work, the decompression needs to be done manually.
system("gunzip --force --keep nfurzeri_genebuild_v1.140704.sqlite.gz")
nfu_db = loadDb('nfurzeri_genebuild_v1.140704.sqlite')
nfu_exonsByGene = exonsBy(nfu_db, by="gene")


# read gene symbol annotation file
nfu.annot.table <- read.csv("nfu.annot.table.csv")
nfu.annot.table$X <- NULL
row.names(nfu.annot.table) <- nfu.annot.table$Nfu.Gene.ID
# replace empty entries by NA
nfu.annot.table$Nfu.Gene.symbol.nfu2015[nfu.annot.table$Nfu.Gene.symbol.nfu2015 == ""] <- NA
nfu.annot.table$Nfu.Gene.symbol.garhum[nfu.annot.table$Nfu.Gene.symbol.garhum == ""] <- NA


###########################
# 1. create data sets     #
###########################

# meta data
grz_runs.meta = read.table(paste(data_path,"GRZ_GONADS.meta.csv",sep="/"),sep="\t",header=TRUE)
grz_runs.meta$file = paste("/misc/vulpix/data/NFU_ANNO_20130926/rnaseq/",grz_runs.meta$run,".bam",sep="")
grz_runs.meta$age = factor(grz_runs.meta$age,levels=c("10dpf","0d","3d","3m"))
grz_runs.meta$sex = factor(grz_runs.meta$sex,levels=c("male","female"))
grz_runs.meta$run = factor(grz_runs.meta$run,levels=unique(grz_runs.meta$run))
grz_runs.meta$tissue = factor(grz_runs.meta$tissue,levels=c("embryo","trunk","testis","ovary"))


# download the counts from GEO
download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263626/suppl/GSE263626%5Fcounts%2Ecsv%2Egz",destfile = "GSE263626_counts.csv.gz",  method='curl')
# read the counts
counts <- read.table(file = "GSE263626_counts.csv.gz", sep = ",", header = TRUE, stringsAsFactors = FALSE)
rownames(counts) <- counts$X
counts$X <- NULL
counts <- as.matrix(counts)

# generate summarizedExperiment from counts
summarised_experiment <- SummarizedExperiment(nfu_exonsByGene, assay = list("counts" = counts), colData = DataFrame(grz_runs.meta))
colnames(summarised_experiment) <- basename(grz_runs.meta$file)



# normalized counts for all samples
ds <- DESeqDataSet(summarised_experiment,design = ~ sex)
ds <- estimateSizeFactors(ds)
data.df.normcounts <- data.frame(counts(ds, normalized = TRUE))
rm(ds)
colnames(data.df.normcounts) <- gsub(".bam","",colnames(data.df.normcounts))
data.df.normcounts <- merge(data.df.normcounts, nfu.annot.table, by = "row.names", all.x = TRUE, sort = FALSE)
data.df.normcounts <- data.df.normcounts[,c("Row.names", "Nfu.Gene.symbol.garhum",as.character(grz_runs.meta$run))]
colnames(data.df.normcounts) <- gsub("Row.names","ID",colnames(data.df.normcounts))
colnames(data.df.normcounts) <- gsub("Nfu.Gene.symbol.garhum","symbol",colnames(data.df.normcounts))
write.csv(data.df.normcounts, file="Supplementary_Data_1_counts_normalized.csv",row.names=FALSE,na="")


###########################
# 2. run DESeq2 & store   #
###########################

# dpf10_male_vs_female_DEG - results and data export

dpf10_summarised_experiment = summarised_experiment[, summarised_experiment$age == "10dpf"]
dpf10_deseq = DESeq(DESeqDataSet(dpf10_summarised_experiment,design = ~ sex))
dpf10_deseq_results = as.data.frame(results(dpf10_deseq))
dpf10_deseq_results <- merge(dpf10_deseq_results, nfu.annot.table, by = "row.names", all.x = TRUE, sort = FALSE)
dpf10_deseq_results <- dpf10_deseq_results[,c("Row.names", "Nfu.Gene.symbol.garhum", "baseMean","log2FoldChange","pvalue","padj")]
colnames(dpf10_deseq_results) <- c("gene", "symbol", "DESeq_baseMean","DESeq_log2FC","DESeq_pvalue","DESeq_FDR")

write.table(dpf10_deseq_results[order(dpf10_deseq_results$DESeq_FDR),], file="Supplementary_Data_1_dpf10_female_vs_male_DEGs.csv",sep="\t",quote=FALSE,row.names=FALSE,na="")


# d0_male_vs_female_DEG - results and data export

d0_summarised_experiment = summarised_experiment[, summarised_experiment$age == "0d"]
d0_deseq = DESeq(DESeqDataSet(d0_summarised_experiment,design = ~ sex))
d0_deseq_results = as.data.frame(results(d0_deseq))
d0_deseq_results <- merge(d0_deseq_results, nfu.annot.table, by = "row.names", all.x = TRUE, sort = FALSE)
d0_deseq_results <- d0_deseq_results[,c("Row.names", "Nfu.Gene.symbol.garhum", "baseMean","log2FoldChange","pvalue","padj")]
colnames(d0_deseq_results) <- c("gene", "symbol", "DESeq_baseMean","DESeq_log2FC","DESeq_pvalue","DESeq_FDR")

write.table(d0_deseq_results[order(d0_deseq_results$DESeq_FDR),], file="Supplementary_Data_1_d0_female_vs_male_DEGs.csv",sep="\t",quote=FALSE,row.names=FALSE,na="")



# d3_male_vs_female_DEG - results and data export

d3_summarised_experiment = summarised_experiment[, summarised_experiment$age == "3d"]
d3_deseq = DESeq(DESeqDataSet(d3_summarised_experiment,design = ~ sex))
d3_deseq_results = as.data.frame(results(d3_deseq))
d3_deseq_results <- merge(d3_deseq_results, nfu.annot.table, by = "row.names", all.x = TRUE, sort = FALSE)
d3_deseq_results <- d3_deseq_results[,c("Row.names", "Nfu.Gene.symbol.garhum", "baseMean","log2FoldChange","pvalue","padj")]
colnames(d3_deseq_results) <- c("gene", "symbol", "DESeq_baseMean","DESeq_log2FC","DESeq_pvalue","DESeq_FDR")

write.table(d3_deseq_results[order(d3_deseq_results$DESeq_FDR),], file="Supplementary_Data_1_d3_female_vs_male_DEGs.csv",sep="\t",quote=FALSE,row.names=FALSE,na="")
