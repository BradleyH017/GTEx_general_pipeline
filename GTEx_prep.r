# Bradley May 2022
# General pipeline for fine mapping using GTEx - any gene in a single script
# This script will be called by ...
# Made generalisable from the original fine mapping in the colon Transverse
# Bradley September 2021
# Read and prep of GTEx data, downloaded by Bioinformatic Resources
# GTEx_prep.r
# Must be done in version 4.1.0
# Apart from the expression matrix retrieved from Bioinformatics Resources, all other data was downloaded from https://gtexportal.org/home/datasets

# Generating the conda environment (from the command line)
#conda create --prefix /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/GTEx/conda_env/gtex_analysis_conda R=4.1 bioconductor-edger r-rnomni r-dplyr  bioconductor-pcatools
#source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/GTEx/conda_env/gtex_analysis_conda


# Set up
library(CePa)
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/GTEx")


# (First time) Load the expression data from BioinformaticsResources
all_expr <- read.gct("/exports/igmm/eddie/BioinformaticsResources/gtex/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
# Load the sample meta data
meta <- read.table("data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep = "\t", row.names = 1, header = T, quote = "", stringsAsFactors = FALSE)

# In this version of the script, want to read in the different options for the site of which we should be looking for eQTLs
library(optparse)
option_list = list(make_option(c("-s", "--site"), action = "store", default = NA, type ="character", help="cluster number"))
opt = parse_args(OptionParser(option_list=option_list))
s = opt$s;
print(paste("The option for site to study is", s, sep = " "))
gene = opt$g
print(paste("The option for gene to study is", gene, sep = " "))

# Now subset for this
meta$SMTSD <- gsub("\\ ", "", meta$SMTSD)
ct_meta <- meta[meta$SMTSD == s, ]
# Subset the expression
rownames(ct_meta) <- gsub("\\-", "\\.", rownames(ct_meta))
keep <- intersect(rownames(ct_meta), colnames(all_expr))
expr <- all_expr[,colnames(all_expr) %in% keep]
ct_meta <- ct_meta[rownames(ct_meta) %in% keep,]

# Read in the Subject phenotypes, so can generate meta with age etc
meta <- read.table("data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", sep = "\t", row.names = 1, header = T, quote = "", stringsAsFactors = FALSE)
# Make patient a variable on the expr to rename the samples so that they match the phenotype meta
expr <- data.frame(t(expr))
temp <- unlist(strsplit(rownames(expr), "\\."))[c(T, T, F, F, F)]
tot <- nrow(ct_meta)
patient <- rep("", tot)
for(p in 1:tot){
	patient[p] <- paste(temp[(2*p)-1], temp[(2*p)], sep = ".")
}
expr$patient <- patient
rownames(expr) <- expr$patient
expr <- expr[,-which(colnames(expr) == "patient")]
# All are unique
# make sure they match and reorder
rownames(meta) <- gsub("\\-", "\\.", rownames(meta))
meta <- meta[match(rownames(expr), rownames(meta)),]

# Check whether expression is given per transcript or per gene
test <- colnames(expr)
test2 <- unlist(strsplit(test, "\\."))[c(T,F)]
length(test2)
length(unique(test2))
# 56200 genes, 56156 are non-unique
dups <- levels(factor(test2[duplicated(test2)]))

# Find the names of the duplicates in the meta
gene_dups <- grep(paste(dups,collapse="|"), test, value = T)

# All of the duplicated genes have a 'PAR_Y' gene version
# This stands for pseudo-autosomal region of the Y chromosome
# The counts for these genes are all 0 across every sample
# These can therefore be removed and gene names generalised to remove the '.x'
expr <- expr[,-grep('PAR_Y', colnames(expr))]
colnames(expr) <- unlist(strsplit(colnames(expr), "\\."))[c(T,F)]

# Alter the meta. Not interested in reason of death, need to treat age like a continuous variable, need to binarise sex
meta <- meta[,-which(colnames(meta) == "DTHHRDY")]
meta$AGE <- gsub("20-29", "25", meta$AGE)
meta$AGE <- gsub("30-39", "35", meta$AGE)
meta$AGE <- gsub("40-49", "45", meta$AGE)
meta$AGE <- gsub("50-59", "55", meta$AGE)
meta$AGE <- gsub("60-69", "65", meta$AGE)
meta$AGE <- gsub("70-79", "75", meta$AGE)
meta$SEX1 <- rep(0, nrow(meta))
meta$SEX2 <- rep(0, nrow(meta))
for(p in 1:nrow(meta)){
	if(meta$SEX[p] == 1){
		meta$SEX1[p] = 1
	} else {
		meta$SEX2[p] = 1
	}
}
meta <- meta[,-which(colnames(meta) == "SEX")]

# Derive the batch variable and add this to the meta
batch_meta <- ct_meta[, which(colnames(ct_meta) == "SMGEBTCH")]
batch_df <- data.frame(patient = rownames(ct_meta), BATCH = batch_meta)
temp <- unlist(strsplit(batch_df$patient, "\\."))[c(T, T, F, F, F)]
tot <- nrow(batch_df)
patient <- rep("", tot)
for(p in 1:tot){
	patient[p] <- paste(temp[(2*p)-1], temp[(2*p)], sep = ".")
}
batch_df$patient <- patient
meta$patient <- rownames(meta)

# Remove samples missing batch
batch_df <- batch_df[-which(batch_df$BATCH == ""),]

# If more than one batch, take the first
for(r in 1:nrow(batch_df)){
	if(length(unlist(strsplit(batch_df$BATCH[r], "\\,"))) == 2){
		batch_df$BATCH[r] <- unlist(strsplit(batch_df$BATCH, "\\,"))[c(T,F)]
	}
	if(length(unlist(strsplit(batch_df$BATCH[r], "\\,"))) == 3){
		batch_df$BATCH[r] <- unlist(strsplit(batch_df$BATCH, "\\,"))[c(T, F, F)]
	}
	if(length(unlist(strsplit(batch_df$BATCH[r], "\\,"))) == 4){
		batch_df$BATCH[r] <- unlist(strsplit(batch_df$BATCH, "\\,"))[c(T, F, F, F)]
	}
}
batch_bin <- data.frame(model.matrix(~0 + BATCH, data = batch_df))
batch_bin$patient <- batch_df$patient

# Now merged with the meta
meta <- meta[meta$patient %in% batch_bin$patient, ]
meta <- merge(meta, batch_bin, by = "patient")

# Now subset the expression to remove the missing sample
rownames(meta) <- meta$patient
meta <- meta[,-which(colnames(meta) == "patient")]
expr <- expr[rownames(expr) %in% rownames(meta)]
# Remove any differences
expr <- expr[-which(rownames(expr) %in% setdiff(rownames(expr), rownames(meta))),]
# Left with 405 samples


# Now remove samples that are missing WGS genotype calls - Leaves us with 367 samples
WGS_samps <- colnames(read.delim("data/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/header_v8_WGS.vcf"))
expr <- expr[rownames(expr) %in% WGS_samps,]
meta <- meta[rownames(meta) %in% WGS_samps,]

# Save the expression and meta
write.csv(expr, paste("data", s, gene, "expr_counts.csv", sep = "/"))
write.csv(meta, paste("data", s, gene, "meta_ready.csv", sep = "/"))


#~~~~~~~~~ 1. Prepping the counts for the analysis
#~~~~~~  Remove the lowly expressed genes from the raw expression - Done this as have done for GTEx
library('edgeR')
expr <- data.frame(t(expr))
tpm <- cpm(expr, method = "TMM", log = F)
# Remove genes with less than 0.1 counts in at least 20% of samples
samp20pc <-  floor(ncol(tpm)*.20)
tmat <- apply(tpm, 1, function(x) { sum(x > 0.1)})
keepGenesT <- names(which(tmat >= samp20pc))
print(paste("Number of genes kept with > 0.1 TPM in at least 20% of samples = ", length(keepGenesT)))
# Also want to keep only genes with at lease 5 reads in 20% of samples
cmat <- apply(expr, 1, function(x) { sum(x >= 6)})
keepGenesC <- names(which(cmat >= samp20pc))
print(paste("Number of genes kept with >= 6 reads in at least 20% of samples = ", length(keepGenesC)))
keepGenes <- intersect(keepGenesT, keepGenesC)
print(paste("Number of genes kept meeting both conditions = ", length(keepGenes)))
# now subset the expression for this
expr <- expr[rownames(expr)  %in% keepGenes,]

#~~~~~~ TMM (non-log) and then INT the expression data
dge<-DGEList(
  counts=expr,
  samples=meta
)
dge <- calcNormFactors(dge, method = "TMM")
d.nonlog.norm <- cpm(dge, log=F, normalized.lib.sizes = T)
library('RNOmni')
# Direct inverse normal transformation
INT <- apply(d.nonlog.norm, 1, RankNorm)

# Save the INT expression values ready for PEER
write.csv(INT, paste("data", s, gene, "expr_TMM_INT.csv", sep = "/"))
