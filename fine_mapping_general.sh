# !/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=192:00:00
#$ -l h_vmem=50G
#$ -N FM_genes_generally
#$ -M bradley.harris@ed.ac.uk
#$ -m baes


# Bradley May 2022
# Extracting genotype calls for GTEx, any gene, any site
gName="DDX5"
sName="Colon-Transverse"
chr="17"

echo "The option for gene is" $gName
echo "The option for site is" $sName
echo "The option for chr is" $chr

# For module loading
. /etc/profile.d/modules.sh

# 1.  Search UCSC genome browser to get the genomic coordinates (from the command line). Save this as the start and end of the window as variables to use later
# Installing esearch for downloading variant databases on the command line
# Following: https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Getting_Started
#sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
#y
export PATH=${PATH}:${HOME}/edirect

# Downloading the genomic coordinates
cd /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/GTEx/data
esearch -db gene -query "$gName[Gene Name]) AND human [ORGN]" | efetch -format docsum | xtract -pattern GenomicInfoType -element -ucsc-based ChrStart ChrStop > start_stop.txt
grange_suff="_gene_range.txt"
grange_f="${gName}${grange_suff}"

# Make directory to save into
mkdir $sName/$gName
awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' start_stop.txt > $sName/$gName/$grange_f
window=1000000
min=$(cat $sName/$gName/$grange_f | sort -nk1,1 | head -1)
min_range=$(($min-$window))
max=$(cat $sName/$gName/$grange_f | sort -nk1,1 | tail -1)
max_range=$(($max+$window))
echo "Found the genomic coordinates of gene and saved"

# 2.  Obtaining the rsids for each snp in this region
suff="_snps.txt"
snp_f_name="${gName}${suff}"
# Test
# test_chr=11
# test_min=111,186,401
# test_max=111,286,401
# esearch -db snp -query "(11[Chromosome]) AND $test_min:$test_max[Base Position]" | efetch -format docsum > test.xml
esearch -db snp -query "($chr[Chromosome]) AND $min_range:$max_range[Base Position]" | efetch -format docsum | xtract -pattern DocumentSummary -sep "|" -element Id CHRPOS > $sName/$gName/$snp_f_name
echo "Found the snps in the genomic window and saved"

# Extract the interesting regions
module load igmm/apps/bcftools/1.9
cd /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/GTEx/data/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1
chr_var="chr"
chr_sub="${chr_var}${chr}"
chr_f_pref="GTEX_v8_"
chr_f_sub="_WGS.vcf"
chr_f_name="${chr_f_pref}${chr_sub}${chr_f_sub}"
bcftools view GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz --regions $chr_sub > ../$chr_f_name
echo "Extracted the chromosome from GTEx vcf"

# Compress these
cd ..
zip_suff=".gz"
chr_f_gz="${chr_f_name}${zip_suff}"
bgzip $chr_f_name
tabix $chr_f_gz -f
echo "bgzipped and tabixed"

# Extract the region of interest
ref_min=$(printf "%'.0f\n" $min_range)
ref_max=$(printf "%'.0f\n" $max_range)
ref_range="${chr_sub}:${ref_min}-${ref_max}"
vcf_midvar="WGS_"
vcf_suff=".vcf"
vcf_name="${chr_f_pref}${vcf_midvar}${gName}${vcf_suff}"
tabix $chr_f_gz $ref_range > $sName/$gName/$vcf_name
echo "Window of snps extracted"

# Append the header (lost at some point)
head_suff="_header.vcf"
head_f_name="${chr_f_pref}${vcf_midvar}${gName}${head_suff}"
cat GTEX_v8_WGS_all_header.vcf $sName/$gName/$vcf_name > $sName/$gName/$head_f_name
echo "Appened header"
echo "~~~~~~~~ DONE"

# remove intermediate file
rm $sName/$gName/$vcf_name

######### Obtaining genotypes and snps as now been completed
######### Can now do the preparation for GTEx
module load anaconda/5.3.1
source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/GTEx/conda_env/gtex_analysis_conda
cd /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/GTEx/scripts/general_pipeline
echo "Starting the GTEx prep"
Rscript GTEx_prep.r -s $sName -g $gName
echo "Done the GTEx prep"

######### Can then do the PEER analysis
source deactivate
module unload anaconda/5.3.1
module unload igmm/apps/R/4.0.2.gcc.9.4.0 
module load igmm/apps/R/3.3.0
echo "Starting the PEER analysis"
Rscript GTEx_PEER.r -s $sName -g $gName
echo "Done the PEER analysis"


######### Can then do the eQTL Analysis
module unload igmm/apps/R/3.3.0
module load igmm/apps/R/4.0.2
echo "Starting eQTL analysis"
Rscript GTEx_eQTL.r -s $sName -g $gName -c $chr
echo "Done the eQTL analysis"
