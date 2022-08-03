##### Bradley July 2022
##### Conditional analysis for Shroom2 following general GTEx pipeline

## First need to extract the list of variants detected in the analysis.
cd /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/GTEx/data/Colon-Transverse/SHROOM2
module load igmm/apps/bcftools/1.9
bgzip GTEX_v8_WGS_SHROOM2_header.vcf
tabix GTEX_v8_WGS_SHROOM2_header.vcf.gz
tail -n +2 /gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/GTEx/results/Colon-Transverse/SHROOM2/MatrixeQTL/final_variants.txt > vars_used_expr.txt
tabix GTEX_v8_WGS_SHROOM2_header.vcf.gz chrX:8,787,755-10,948,474 --print-header > GTEX_v8_WGS__expr_header_SHROOM2.vcf

# Then make a bed file from the variants we have and are using. Note: this also makes a .bim and .bam file for these
module load igmm/apps/plink/2.03
datadir=/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/GTEx/data/Colon-Transverse/SHROOM2
resdir=/gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/GTEx/results/Colon-Transverse/SHROOM2/ConditionalAnalysis
cd $resdir
plink2 --vcf $datadir/GTEX_v8_WGS__expr_header_SHROOM2.vcf --extract $datadir/vars_used_expr.txt --make-bed

# Make the .ma file
module load igmm/apps/R/4.0.2
R
setwd("/gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/GTEx/results/Colon-Transverse/SHROOM2/ConditionalAnalysis")
genes <- c("SHROOM2")
# Load
res <- read.csv("../MatrixeQTL/cis_results.csv")
# Calculate the SE: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/faq.html#se
res$se <- res$beta / res$statistic
# Extract the bits we want
res <- res[, colnames(res) %in% c("snps", "beta", "FDR", "REF", "ALT", "Frequency", "se")]
# Rename
colnames(res) <- c("SNP", "p", "b", "A1", "A2", "freq", "se")
# Reorder
res <- res[, match(c("SNP", "A1", "A2", "freq", "b", "se", "p"), colnames(res))]
res$N <- rep(367, nrow(res))
res$SNP <- gsub("b381", "b38", res$SNP)
#Save
write.table(res, "SHROOM2_ciseQTLs.ma", quote=F, sep = "\t", row.names=F)

# Now run conditional on each interseting snp
int_snps <- c("chrX:9783434", "chrX:9784504", "chrX:9786352", "chrX:9791765", "chrX:9789973", "chrX:9790796", "chrX:9795858", "chrX:9797979")
int_snps <- gsub("\\:", "\\_", int_snps)
int_snps_df <- data.frame(search = int_snps, var=rep("", length(int_snps)))
for(v in seq_along(int_snps)){
  tryCatch({
    int_snps_df$var[v] = res[grep(int_snps[v], res$SNP),]$SNP
  } , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
# chrX_9789973 and chrX_9786352 missing

int_snps_df <- int_snps_df[int_snps_df$var != "", ]

# Add top var
int_save <- int_snps_df$var
int_save <- c(int_save, res$SNP[1])
write.csv(int_save, "interseting_variants.txt", row.names=F, col.names=F, quote=F)

# Remove first rows

# Rename the chromsome number in the plink files. This needs to read 23 instead of X
sed -i 's/X/23/' plink2.bim

# Annotating the sex information for each sample
R
fam <- read.delim("plink2.fam", sep = "\t", header=F)
# From https://www.cog-genomics.org/plink/1.9/formats#fam
colnames(fam) <- c("FamID", "Sample", "FFamID", "MFamID", "Sex_code", "Phenotype")
meta <- read.csv("../../../../data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", sep = "\t", row.names = 1, header = T, quote = "", stringsAsFactors = FALSE)
meta_sex <- data.frame(Sample=rownames(meta), Sex_code=meta$SEX)
meta_sex <- meta_sex[meta_sex$Sample %in% fam$Sample,]
# From GTEx data dictionary, 1=Male, 2=Female. Same as is required by gcta. So just replace this
meta_sex <- meta_sex[match(fam$Sample, meta_sex$Sample),]
all(meta_sex$Sample == fam$Sample)
fam$Sex_code <- meta_sex$Sex_code
# Overwrite fam
write.table(fam, "plink2.fam", sep = "\t", quote=F, row.names=F, col.names=F)

##### Now run the conditional analysis (back in linux)
# From http://ibg.colorado.edu/cdrom2017/evansL/GCTA_practical/references/Yang%20et%20al.%202011%20GCTA_UserManual_v1.0.pdf
# The default number of autosome chromosomes is 22, and the X-chromosome is automatically detected as chromosome 23.
module load igmm/apps/gcta/1.91.4beta
while read p; do
  echo "$p" > snp.txt
  echo "WORKING ON"
  echo "$p"
  gcta64 --bfile plink2 --autosome-num 22 --chr 23 --maf 0.01 --cojo-file SHROOM2_ciseQTLs.ma --cojo-actual-geno --cojo-cond snp.txt --out $p
done <interseting_variants.txt

## Summarising the analysis
R
snps <- readLines("interseting_variants.txt")
cond_res <- vector("list", length = length(snps))
names(cond_res) <- snps
conds_res_int <- vector("list", length=length(snps))
names(conds_res_int) <- snps
for(snp in seq_along(snps)){
  # Read in
  cond_res[[snp]] <- read.delim(paste0(snps[snp], ".cma.cojo"), sep = "\t", header=T)
  # Reduce for significant
  cond_res[[snp]] <- cond_res[[snp]][cond_res[[snp]]$pC < 0.05,]
  # Reduce to complete cases  to exclude NA
  cond_res[[snp]] <- cond_res[[snp]][complete.cases(cond_res[[snp]]),]
  # Subset for the interseting snps in new list
  conds_res_int[[snp]] <- cond_res[[snp]][cond_res[[snp]]$SNP %in% snps,]
}
sapply(cond_res, nrow)
# Each of these snps has at least 983 independent signals

sapply(conds_res_int, nrow)
# Some of the interseting snps are independent signals of one another

# Plot a matrix of these
mat <- matrix(nrow=length(snps), ncol=length(snps))
dimnames(mat) <- list(snps, snps)
for(r in 1:nrow(mat)){
  for(c in 1:ncol(mat)){
    if(snps[r] %in% conds_res_int[[c]]$SNP & snps[c] %in% conds_res_int[[r]]$SNP){
      mat[r,c] <- 1
    }
  }
}
# save
write.csv(mat, "intersection_matrix_of_independent_signals_from_interseting_snps.csv")

# Rename and save
rs <- c("rs5934683", "rs773016882", "rs5934684", "rs957490", "rs2732875", "rs5934685", "rs1175673381")
dimnames(mat) <- list(rs,rs)
write.csv(mat, "intersection_matrix_of_independent_signals_from_interseting_snps_rsid.csv")
