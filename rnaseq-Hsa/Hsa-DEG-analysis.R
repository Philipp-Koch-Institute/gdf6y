## This script was originally run under R 4.1.3
## See README.md for instructions on how to re-construct the original environment.

# Note: the following lines were only used by the authors to initialize renv. 
# They should be left deactivated (= as comments)
#renv::init()
#renv::dependencies()
#renv::snapshot()


# load required packages
library(DESeq2)
library(openxlsx)
library(AnnotationHub)
library(AnnotationDbi)


### Annotation loading ###
ah <- AnnotationHub()
ah.query <- query(ah, c("EnsDb", "Homo sapiens", ""))
EnsDb <- ah.query[[length(ah.query)]]

# print some information on the annotation
print(metadata(EnsDb))


### Data loading ###

# The download of the counts works only when the GEO entry finally was published.
# download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233992/suppl/GSE233992%5Fcounts%2Ecsv%2Egz",destfile = "GSE233992_counts.csv.gz")
# for the moment, read-in the provided csv file.
data.df.counts <- read.table(file = "GSE233992_counts.csv.gz", sep = ",", header = TRUE)
row.names(data.df.counts) <- data.df.counts$X
data.df.counts$X <- NULL

# read the metadata.csv and convert every column to factor
metadata <- read.table(file = "metadata.csv", header = TRUE, sep = ",")

# read the coldata.txt and convert every column to factor
coldata <- read.table(file = "coldata.txt", header = TRUE, stringsAsFactors = TRUE)


### Data preparation ###

# normalized counts for all samples
ds <- DESeqDataSetFromMatrix(data.df.counts,coldata, design = formula("~ treatment"))
ds <- estimateSizeFactors(ds)
data.df.normcounts <- data.frame(counts(ds, normalized = TRUE))
rm(ds)


### DEG analysis (2 parts) ###

### ### ### 1st analysis: "treatment_gdf6Y_vs_EV" ### ### ###
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# select groups for this analysis remove un-used levels from factors in one go
coldata.selected.1 <- droplevels(coldata[coldata$treatment == "gdf6Y" | coldata$treatment == "EV",])

# filter counts table based on selected groups
data.df.counts.1 <- data.df.counts[,as.character(coldata.selected.1$sample)]

# build the data structure 
dds.1 <- DESeqDataSetFromMatrix(data.df.counts.1,coldata.selected.1, design = formula("~ treatment"))

# remove "nearly" empty lines
dds.1 <- dds.1[rowSums(counts(dds.1)) > 1,]

# run DEseq2 
dds.1 <- DESeq(dds.1, betaPrior = FALSE)


# extract results & shrink l2fc
res.unshrunk.1 <- results(dds.1, contrast = c("treatment", "gdf6Y", "EV")) 
coefficient.1 <- "treatment_gdf6Y_vs_EV"
res.1 <- lfcShrink(dds = dds.1, coef = coefficient.1, res = res.unshrunk.1)

# sum up the results
summary(res.1, alpha = 0.05)

# map the gene symbol and biotype according to the Ensembl ID
res.1$symbol <- mapIds(EnsDb, keys = row.names(res.1), column = "SYMBOL", keytype = "GENEID", multiVals = "first")
res.1$biotype <- mapIds(EnsDb, keys = row.names(res.1), column = "GENEBIOTYPE", keytype = "GENEID")


# order results table by padj (lowest on top)
res.ordered.1 <- res.1[order(res.1$padj),]

# calculate mean of normalized counts per treatment group
baseMeanPerLvl.1 <- sapply( c("gdf6Y", "EV"), function(lvl) rowMeans(counts(dds.1,normalized = TRUE)[,dds.1[["treatment"]] == lvl]))

# merge tables
res.ordered.df.1 <- merge(as.data.frame(res.ordered.1), baseMeanPerLvl.1, by = "row.names", sort = FALSE)

# change some column names
colnames(res.ordered.df.1)[which(colnames(res.ordered.df.1) == "Row.names")] <- "ID"
colnames(res.ordered.df.1)[which(colnames(res.ordered.df.1) == "gdf6Y")] <- "gdf6Y.Mean"
colnames(res.ordered.df.1)[which(colnames(res.ordered.df.1) == "EV")] <- "EV.Mean"

# reorder columns
res.ordered.df.1 <- res.ordered.df.1[,c("ID", "baseMean", "gdf6Y.Mean", "EV.Mean", "log2FoldChange", "lfcSE", "pvalue", "padj", "symbol", "biotype")]


### ### ### 2nd analysis: "treatment_gdf6Y_vs_del9" ### ### ###
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# select groups for this analysis remove un-used levels from factors in one go
coldata.selected.2 <- droplevels(coldata[coldata$treatment == "gdf6Y" | coldata$treatment == "del9",])

# filter counts table based on selected groups
data.df.counts.2 <- data.df.counts[,as.character(coldata.selected.2$sample)]

# build the data structure 
dds.2 <- DESeqDataSetFromMatrix(data.df.counts.2,coldata.selected.2, design = formula("~ treatment"))

# remove "nearly" empty lines
dds.2 <- dds.2[rowSums(counts(dds.2)) > 1,]

# run DEseq2 
dds.2 <- DESeq(dds.2, betaPrior = FALSE)


# extract results & shrink l2fc
res.unshrunk.2 <- results(dds.2, contrast = c("treatment", "gdf6Y", "del9")) 
coefficient.2 <- "treatment_gdf6Y_vs_del9"
res.2 <- lfcShrink(dds = dds.2, coef = coefficient.2, res = res.unshrunk.2)

# sum up the results
summary(res.2, alpha = 0.05)

# map the gene symbol and biotype according to the Ensembl ID
res.2$symbol <- mapIds(EnsDb, keys = row.names(res.2), column = "SYMBOL", keytype = "GENEID", multiVals = "first")
res.2$biotype <- mapIds(EnsDb, keys = row.names(res.2), column = "GENEBIOTYPE", keytype = "GENEID")


# order results table by padj (lowest on top)
res.ordered.2 <- res.2[order(res.2$padj),]

# calculate mean of normalized counts per treatment group
baseMeanPerLvl.2 <- sapply( c("gdf6Y", "del9"), function(lvl) rowMeans(counts(dds.2,normalized = TRUE)[,dds.2[["treatment"]] == lvl]))

# merge tables
res.ordered.df.2 <- merge(as.data.frame(res.ordered.2), baseMeanPerLvl.2, by = "row.names", sort = FALSE)

# change some column names
colnames(res.ordered.df.2)[which(colnames(res.ordered.df.2) == "Row.names")] <- "ID"
colnames(res.ordered.df.2)[which(colnames(res.ordered.df.2) == "gdf6Y")] <- "gdf6Y.Mean"
colnames(res.ordered.df.2)[which(colnames(res.ordered.df.2) == "del9")] <- "del9.Mean"

# reorder columns
res.ordered.df.2 <- res.ordered.df.2[,c("ID", "baseMean", "gdf6Y.Mean", "del9.Mean", "log2FoldChange", "lfcSE", "pvalue", "padj", "symbol", "biotype")]



### Data export ###

# add and adjust columns of normalized counts
data.df.normcounts$ID <- rownames(data.df.normcounts)
data.df.normcounts$symbol <- mapIds(EnsDb, keys = row.names(data.df.normcounts), column = "SYMBOL", keytype = "GENEID", multiVals = "first")
data.df.normcounts <- data.df.normcounts[,c("ID","symbol",as.character(coldata$sample))]

# save results objects for further usage - un-comment if desired
#save(list = c("res.ordered.1"), file = paste0("DESeq2_",coefficient.1,"_res.ordered.RData"))
#save(list = c("res.ordered.2"), file = paste0("DESeq2_",coefficient.2,"_res.ordered.RData"))

# export all four sheets to one Excel file
write.xlsx(list("A) Meta data" = metadata,
                "B) Normalized counts" = data.df.normcounts,
                "C) gdf6Y vs EV" = res.ordered.df.1,
                "D) gdf6Y vs del9" = res.ordered.df.2), file = paste0("Supplementary_Data_3_ARichter2024.xlsx"), rowNames = FALSE)
