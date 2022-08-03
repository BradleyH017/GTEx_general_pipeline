# Bradley May 2022
# General pipeline for fine mapping in GTEx - any gene, any site
#Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
myPaths <- .libPaths()
myPaths <- c(myPaths, "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
myPaths <- c(myPaths[2], myPaths[1])
.libPaths(myPaths)
setwd("GTEx")
library('MatrixEQTL', lib="/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")


# IN this instance, doing this for multiple sites at the same time
library(optparse)
option_list = list(make_option(c("-s", "--site"), action = "store", default = NA, type ="character", help="cluster number"),
										make_option(c("-g", "--gene"), action = "store", default = NA, type ="character", help="gene name"),
										make_option(c("-c", "--chromosome"), action = "store", default = NA, type ="character", help="chromosome"))
opt = parse_args(OptionParser(option_list=option_list))
s = opt$s;
print(paste("The option for site to study is", s, sep = " "))
gene = opt$g
print(paste("The option for gene to study is", gene, sep = " "))

# Rename the s variable if using any with brackets
s<-gsub("\\(", "", s)
s<-gsub("\\)", "", s)


# Define path
pathOut <- paste("results", gsub("\\ ", "", s), gene, "PEER", sep = "/")

# Don't want to run again and again, so stop here if already ran
if (file.exists(paste("results", gsub("\\ ", "", s), gene, "MatrixeQTL", sep = "/"))) {stop("Already done, so stop here")
  } else { #continue the script
print("Continuing script")
  }

#~~~~~ Load data
#1. Expression
expression <- read.csv(paste(pathOut, "Peer_residuals_age_gender_batch.csv", sep = "/"), row.names=1)
# Need genes as rows, samples as columns
expression <- data.frame(t(expression))
colnames(expression) <- gsub("\\.", "\\-", colnames(expression))


#2. Genotypes
library(vcfR)
genos <- read.vcfR(paste("data/", s, "/", gene, "/GTEX_v8_WGS_", gene, "_header.vcf", sep = ""))

# Subset for samps in our study
samps <- colnames(expression)
genos <- genos[,c("FORMAT", samps)]

# Subset for variants passing filtering
sub <- which(colnames(genos@fix) == "FILTER")
genos <- genos[genos@fix[,sub] == "PASS", ]

# Calculate maf
maf_genos <- as.data.frame(maf(genos))

# Remove those with missing data and/or maf < 0.01
keep <- as.data.frame(maf_genos[maf_genos[,2] == 0,])
keep <- keep[keep$Frequency > 0.01,]
keep <- rownames(keep)
sub <- which(colnames(genos@fix) == "ID")
genos <- genos[genos@fix[,sub] %in% keep,]

# Redefine pathOut
pathOut <- paste("results", gsub("\\ ", "", s), gene, "MatrixeQTL", sep = "/")
if(file.exists(pathOut) == F){
	dir.create(pathOut)
}

# Now extract the genotype dataframe and convert
gt <- genos@gt[,-1]
rownames(gt) <- genos@fix[,which(colnames(genos@fix) == "ID")]
# Remove the other stuff
simp <- function(x){
	# Take first element
	temp <- unlist(strsplit(x, "\\:"))[c(T,F,F,F,F)]
	return(temp)
}
gt <- apply(gt[,grep("GTEX",colnames(gt))], 2, FUN=simp)
#remotes::install_github("etnite/bwardr")
library(bwardr)
testgenos <- gt2num(as.matrix(gt))
testgenos <- testgenos[[1]]
genos <- testgenos
rm(testgenos)


# Loading the snps - reformat
file=paste("data/", s, "/", gene, "/", gene, "_snps.txt", sep = "")
conv <- read.delim(file, sep="\t", header=F)

colnames(conv) <- c("rsid", "chrpos")
conv <- conv[complete.cases(conv),]
conv <- conv[conv$chrpos != "", ]
conv$pos <- unlist(strsplit(conv$chrpos, ":"))[c(F,T)]
conv$rsid <- paste0("rs", conv$rsid)
# Make sure these are still within the max/min values
maxmin <- as.numeric(readLines(paste("data/", s, "/", gene, "/", gene, "_gene_range.txt", sep = "")))
conv <- conv[conv$pos > min(maxmin) & conv$pos < max(maxmin),]

# test the intersection with genos and make the snpid
genpos <- unlist(strsplit(rownames(genos), "\\_"))[c(F,T,F,F,F)]
length(intersect(genpos, conv$pos))

# Optional here: Do we want to subset for only the variants that have rsids?
use_rs <- F
if(use_rs  == T){
	genos <- as.data.frame(genos_genostemp)
	genos$POS <- unlist(strsplit(rownames(genos), "\\_"))[c(F,T,F,F,F)]
	snpid <- snpid[snpid$POS %in% genos$POS,]
	genos <- genos[match(snpid$POS, genos$POS),]
	genos$ID <- snpid$ID
	rownames(genos) <- genos$ID
	genos <- genos[,-c(which(colnames(genos) == "ID"), which(colnames(genos) == "POS"))]
}


# Save the final list of variant IDs. Somehow, these were added a b381 on the end, not b38, so need to change this for subsetting
vars <- data.frame(variant=rownames(genos))
vars$POS <- unlist(strsplit(vars[,1], "\\_"))[c(F,T,F,F,F)]
vars$variant <- gsub("b381", "b38", vars$variant)
# Save a final list of variant positions
vars$rsid <- rep("", nrow(vars))
for(snp in 1:nrow(vars)){
	if(vars$POS[snp] %in% conv$pos){
		vars$rsid[snp] <- conv[conv$pos == conv$pos[snp],]$rsid
	}
}
write.table(data.frame(vars=vars$variant), paste(pathOut, "final_variants.txt", sep = "/"), sep="\t", row.names=F, quote=F)
write.table(data.frame(vars=vars$POS), paste(pathOut, "final_var_positions.txt", sep = "/"), sep="\t", row.names=F, quote=F)

# Final number of variants before analysis
print(paste("The final number of variable variants in the analysis is", nrow(genos)))

# 4. Download the gene coordinates - Using GRCh38, this is the default for the useEnsembl function
library(BiocManager)
library(biomaRt, lib ="/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
ensembl = useEnsembl(biomart="ensembl",
#	GRCh=38,
	dataset="hsapiens_gene_ensembl",
	host = "http://www.ensembl.org")
listMarts(ensembl)$version[1] # Using ensemble version 104
gene_coords <- getBM(attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "gene_biotype"),  uniqueRows=T, mart=ensembl)
gene_coords$length <- gene_coords$end_position - gene_coords$start_position
head(gene_coords)
# These match


#~~~~~~ Analysis
#Set up the parameters such as selected linear model
useModel = modelLINEAR

#Set the p-value threshold and set up the output file
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

pvOutputThreshold_cis = 1e-2;
pvOutputThreshold_tra = 1e-2;
errorCovariance = numeric();



# Reformat the snps - This are genotype data for each sample in the form of a Sliced Data object
#genos$snp <- rownames(genos)
genos <- as.matrix(genos)
#rownames(genotypes) <- NULL
snps = SlicedData$new()
snps$CreateFromMatrix(genos)
show(snps)



# Reformat the expression - SlicedData object with gene expression information. Should have columns matching those of snps.
genes.use <- intersect(rownames(expression), gene_coords$ensembl_gene_id)
expression <- expression[rownames(expression) %in% genes.use, ]
# Order to match snps
expression <- expression[,match(colnames(genos), colnames(expression))]
all(colnames(expression) == colnames(genos))
expression <- as.matrix(expression)
gene = SlicedData$new()
gene$CreateFromMatrix(expression);
# Check
show(gene)
# normally distributing gene values
for( sl in 1:length(gene) ) {
  mat = gene[[sl]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(gene)+1));
  gene[[sl]] = mat;
}
rm(sl, mat);
show(gene)


# Regerate "gene", have overwitten it here. Needs new variable name too,
gene.use = opt$g

# Prep the covariates - Need to re-load the binarised, numeric matrix
meta_bin <- read.csv(paste("data", gsub("\\ ", "", s), gene.use, "meta_ready.csv",sep = "/"), row.names = 1)
# Rename the rows
meta_bin <- as.matrix(meta_bin)
meta_bin <- t(meta_bin)
#rownames(covariates) <- NULL
cvrt = SlicedData$new()
cvrt$CreateFromMatrix(meta_bin);
show(cvrt)


# Prep the snp coordinates - data.frame with information about SNP locations, with 3 columns - SNP name, chromosome, and position.
chr=opt$c
if(use_rs  == T){
	snpspos <- snpid
	snpspos$chr <- rep(11, nrow(snpspos))
	colnames(snpspos) <- c("rsid", "position", "chr")
	snpspos <- snpspos[,match(c("rsid", "chr", "position"), colnames(snpspos))]
	snpspos[,3] <- as.numeric(snpspos[,3])
} else {
	snpspos <- data.frame(ID=rownames(genos),
		chr=unlist(strsplit(rownames(genos), "\\_"))[c(T,F,F,F,F)],
		position=unlist(strsplit(rownames(genos), "\\_"))[c(F,T,F,F,F)])
	snpspos$chr <- gsub("chr", "", snpspos$chr)
	if(chr != "X" & chr != "Y"){
		snpspos$chr <- as.numeric(snpspos$chr)
	}
	snpspos$position <- as.numeric(snpspos$position)
}



# Prep the gene coordinates - data.frame with information about transcript locations, with 4 columns - the name, chromosome, and positions of the left and right ends.
gene_coords <- gene_coords[gene_coords$ensembl_gene_id %in% genes.use, ]
#gene_coords <- gene_coords[order(gene_coords$ensembl_gene_id),]
gene_coords <- gene_coords[,colnames(gene_coords) %in% c("ensembl_gene_id", "chromosome_name", "start_position", "end_position")]
colnames(gene_coords) <- c("gene_name", "chr", "left", "right")
head(gene_coords)



# Set the threshold for local gene-snp interactions - Same as NM
cisDist = 2e6;
gene_coords <- as.data.frame(gene_coords)
gene_coords[,3] <- as.numeric(gene_coords[,3])
gene_coords[,4] <- as.numeric(gene_coords[,4])


# Run the analysis for histogram of pvals only
meh = Matrix_eQTL_main(
snps = snps,
gene = gene,
#cvrt = cvrt,
output_file_name = output_file_name_cis,
pvOutputThreshold = pvOutputThreshold_tra,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos,
genepos = gene_coords,
cisDist = cisDist,
pvalue.hist = 50,
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);


unlink(output_file_name_tra);
unlink(output_file_name_cis);

pdf(file=paste(pathOut, "Histogram_pvals.pdf", sep = "/"))
plot(meh);
dev.off()

#Run the analysis (can do this for no additonal covariates)
me = Matrix_eQTL_main(
snps = snps,
gene = gene,
#cvrt = cvrt,
output_file_name = output_file_name_cis,
pvOutputThreshold = pvOutputThreshold_tra,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos,
genepos = gene_coords,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

pdf(file=paste(pathOut, "QQ_local_and_distant_pvals.pdf", sep = "/"))
plot(me)
dev.off()


##Results
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
#show(me$cis$eqtls)


# Reformat before saving
me$cis$eqtls$marker <- paste(unlist(strsplit(me$cis$eqtls$snps, "\\_"))[c(T,F,F,F,F)],
							unlist(strsplit(me$cis$eqtls$snps, "\\_"))[c(F,T,F,F,F)],
							sep = ":")
# Add rsid for those that have it
me$cis$eqtls$rsid <- rep("", nrow(me$cis$eqtls))
for(snp in 1:nrow(me$cis$eqtls)){
	temp <- unlist(strsplit(me$cis$eqtls$marker[snp], "\\:"))[c(F,T)]
	if(temp %in% conv$pos){
		me$cis$eqtls$rsid[snp] <- paste(conv[conv$pos == temp,]$rsid, collapse = ",")
	} else {
		me$cis$eqtls$rsid[snp] <- me$cis$eqtls$marker[snp]
 		}
}
# Add ref/alt cols
me$cis$eqtls$REF <- unlist(strsplit(me$cis$eqtls$snps, "\\_"))[c(F,F,T,F,F)]
me$cis$eqtls$ALT <- unlist(strsplit(me$cis$eqtls$snps, "\\_"))[c(F,F,F,T,F)]
# Add the proper marker id
me$cis$eqtls$markerID <- paste(me$cis$eqtls$marker, paste(me$cis$eqtls$REF, me$cis$eqtls$ALT, sep = "/"), sep = "_")
me$cis$eqtls$markerID <- gsub("chr", "", me$cis$eqtls$markerID)

# add chromosome
me$cis$eqtls$chr <- unlist(strsplit(me$cis$eqtl$markerID, "\\:"))[c(T,F)]

# Append maf
maf_genos$snps <- rownames(maf_genos)
me$cis$eqtls$snps <- gsub("b381", "b38", me$cis$eqtls$snps)
me$cis$eqtls <- merge(me$cis$eqtls, maf_genos, by = "snps")
me$cis$eqtls <- me$cis$eqtls[order(me$cis$eqtls$FDR),]


# Save all together
write.csv(me$cis$eqtls, paste(pathOut, "cis_results.csv", sep = "/"))

## Extract BMP4
#if(region == "BMP4"){
#	int <- me$cis$eqtls[me$cis$eqtls$gene == "ENSG00000125378",]
#	int$POS <- unlist(strsplit(int$marker, "\\:"))[c(F,T)]
#	head(int)
#	int <- int[order(int$POS), ]
#	write.table(int, paste(pathOut, "BMP4_cis_eQTLs.txt", sep = "/"), sep = "\t", quote=F, row.names=F)
#}

# Save trans
write.csv(me$trans$eqtls, paste(pathOut, "trans_results.csv", sep = "/"))
