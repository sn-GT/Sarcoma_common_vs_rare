module load R/4.1.0
module load intel-compiler/2021.2.0
module load intel-mkl/2021.2.0

# Obtained genotype QC'd files for ISKS, MGRB, UKB (phase 1: discovery)

# Step1: [R] Get common SNPs and matched alleles for BGI + phase 1 (ISKS, MGRB, UKB)

bgi <- fread("/scratch/public/sn0095/BGI-Updated/qc_output.bim")

disc <- fread("/home/552/sn0095/shared/SarcomaAnalysis/Results/AF/Combined_all3_harmonizedAF.bim")
bgi$CHRPOS <- paste0(bgi$V1,":", bgi$V4)
disc$CHRPOS <- paste0(disc$V1,":", disc$V4)

m <- merge(disc, bgi, by = "CHRPOS")
# Keep only those SNPs where the two alleles match in the dataset (i.e. remove allele mismatch or flip cases)
match <- m[m$V5.x == m$V5.y & m$V6.x == m$V6.y ,]
rg <- match[,c("V1.x", "V4.x", "V4.x")]
names(rg) <- c("CHR", "START", "END")
write.table(rg, "/g/data/wq2/users/sn0095/datasets/matched_alleles_bgi2.txt", row.names=F, col.names=T, quote=F,sep = "\t")

#######################################################################################################################

# Step2: [PLINK] Make plink bed files for common snps from both datasets in the same format with variant IDs as chr:pos:A1:A2: 
# Extract common alleles from discovery set phase1, and set var IDs is common format across datasets

plink_path="/home/552/sn0095/Tools/plink2_linux_x86_64_20220814/plink2"
file1="/home/552/sn0095/shared/SarcomaAnalysis/Results/AF/Combined_all3_harmonizedAF"
out="/g/data/wq2/users/sn0095/datasets/Combined_all3_harmonizedAF_varIDs2"
ex="/g/data/wq2/users/sn0095/datasets/matched_alleles_bgi2.txt"
$plink_path --bfile $file1 --extract range $ex --set-all-var-ids @:#:\$a:\$r --new-id-max-allele-len 5 missing --make-bed --out $out --max-alleles 2
$plink_path --bfile $out --freq --out "/g/data/wq2/users/sn0095/datasets/FREQ_Combined_all3_harmonizedAF_varIDs2"

plink_path="/home/552/sn0095/Tools/plink2_linux_x86_64_20220814/plink2"
file1="/scratch/public/sn0095/BGI-Updated/qc_output"
out="/g/data/wq2/users/sn0095/datasets/bgi_matched_varIDs2"
ex="/g/data/wq2/users/sn0095/datasets/matched_alleles_bgi2.txt"
$plink_path --bfile $file1 --extract range $ex --set-all-var-ids @:#:\$a:\$r --new-id-max-allele-len 5 missing --make-bed --out $out --max-alleles 2
$plink_path --bfile $out --freq --out "/g/data/wq2/users/sn0095/datasets/FREQ_bgi_matched_varIDs2"

#######################################################################################################################

# Step 3: [R] Compare allele frequencies of the two datasets and remove outlier SNPs (data harmonization before combining datasets)

bgi <- fread("/g/data/wq2/users/sn0095/datasets/FREQ_bgi_matched_varIDs2.afreq", header = T)
comb <- fread("/g/data/wq2/users/sn0095/datasets/FREQ_Combined_all3_harmonizedAF_varIDs2.afreq", header = T)

x <- merge(bgi, comb, by = "ID")

# ensure allele match for two datasets
x1 <- x[x$ALT.x == x$ALT.y,]

# correlation of AFs
cor.test(x1$ALT_FREQS.x, x1$ALT_FREQS.y)

png("bgi_combined2.png", res= 200, width = 6, height = 4.5, units = "in")
ggplot(x1, aes(x = ALT_FREQS.x, y=ALT_FREQS.y)) + geom_point() + 
  theme_bw() + 
  xlab("BGI") + ylab("DiscoveryCohort") + geom_smooth(method = "lm", color = "grey")
dev.off()

# Fit linear model to get residuals - if resid > 3sdu, remove
m <- x1
fit <- lm(data = m, ALT_FREQS.y ~ ALT_FREQS.x)
m$resid <- resid(fit)

sd2<-3*sd(resid(fit))
m$Group<-ifelse(abs(m$resid)>sd2, "outlier_3sd", "use")

msub <- m[m$Group == "use",]
dim(msub)
cor.test(msub$ALT_FREQS.x, msub$ALT_FREQS.y)
r <- "0.99"
p <- "2e-16"
png("bgi_combined_outlier2.png",  res=200, width=6, height =4.5, units = "in")
p <- ggplot(data = m, aes(x = ALT_FREQS.x,y=ALT_FREQS.y, color=Group)) + 
geom_point()+theme_bw() + 
#ggtitle("", sub = paste0("r=",c," , P-val=", p)) + 
xlab("BGI") + ylab("Validation Cohort") + 
geom_smooth(method = "lm", se = F, col = "grey") + scale_color_manual(values = c("red", "black")) + ylim(0,0.5) + xlim(0,0.5)
print({p})
dev.off()

exc <- m[m$Group == "outlier_3sd",]
#######################################################################################################################


# Step 4: [PLINK] After identifying outlier SNPs, remove them from plink files and merge the harmonized datasets to get final combined data

plink_path="/home/552/sn0095/Tools/plink_linux_x86_64_20220402/plink"
file1="/g/data/wq2/users/sn0095/datasets/Combined_all3_harmonizedAF_varIDs2"
out="/g/data/wq2/users/sn0095/datasets/Combined_all4_merged2"
$plink_path --bfile $file1 --merge-list "/home/552/sn0095/Analysis/BGI/Scripts/list.txt" --make-bed --out $out --biallelic-only strict --allow-no-sex

# PCA
plink_path="/home/552/sn0095/Tools/plink2_linux_x86_64_20220814/plink2"
file1="/g/data/wq2/users/sn0095/datasets/Combined_all4_merged2"
out="/g/data/wq2/users/sn0095/datasets/Combined_all4_merged_PCA"
extr="/home/552/sn0095/shared/SarcomaAnalysis/Results/AF/HarmonizedSNPs_AF_LDintersection_allcohort_BGIuse.txt"
$plink_path --bfile $file1 --extract $extr --pca 10 --out $out

#######################################################################################################################

# Step 5: [R] Remove outlier indivduals based on PCA
# Plot PCA and get zoomed ids tighly clustered European ancestry individuals

x <- read.table("/g/data/wq2/users/sn0095/datasets/Combined_all4_merged_PCA.eigenvec", header = F)
sarcoma <- fread("/home/552/sn0095/shared/SarcomaAnalysis/Results/AF/SarcomaCancer_harmonizedAF.fam", header = F)
isks <- fread("/home/552/sn0095/shared/SarcomaAnalysis/Results/AF/isks_harmonizedAF.fam", header = F)
mgrb <- fread("/home/552/sn0095/shared/SarcomaAnalysis/Results/AF/mgrb_harmonizedAF.fam", header = F)
bgi <- fread("/home/552/sn0095/Analysis/BGI-Updated/bgi2.fam", header = F)


x$GROUP <- "X"
x$GROUP[x$V2 %in% sarcoma$V2] <- "UKB"
x$GROUP[x$V2 %in% isks$V2] <- "ISKS"
x$GROUP[x$V2 %in% mgrb$V2] <- "MGRB"
x$GROUP[x$V2 %in% bgi$V2] <- "BGI"

x2_sub <- x[x$V3 >= -0.019 & x$V3 <= 0.019 & x$V4 >= -0.019 & x$V4 <= 0.019,]

library(ggplot2)
png("/home/552/sn0095/Analysis/BGI/AF/pca_combined2.png", res=300, units = "in", width = 5, height = 5)
ggplot(data = x, aes(x= V3, y=V4, group = GROUP, color = GROUP)) + geom_point(alpha=0.7) + 
theme_bw() + xlab("PC1") + ylab("PC2") + 
 scale_color_manual(values = c("#62825D","#A7D2CB","#F2D388", "#C98474")) + 
 geom_hline(yintercept = -0.019, linetype = "dashed") + 
geom_hline(yintercept = 0.019, linetype = "dashed") + 
geom_vline(xintercept = -0.019, linetype = "dashed") + 
geom_vline(xintercept = 0.019, linetype = "dashed") 
dev.off()

#######################################################################################################################

# Step 6: [R] Create phenotype and covariate files for GWAS

covall <- read.table("/g/data/wq2/users/sn0095/datasets/COVALL.txt")
m <- merge(x2_sub, covall, by = "IID")
mf <- m[,c("V1", "IID", "V3", "V4", "V5", "V6", "V7", "Sex")]
names(mf) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "Sex")
write.table(mf, "/g/data/wq2/users/sn0095/datasets/COV_all4_PCSEX_FINAL2.txt", row.names=F, col.names=T, quote=F, sep = "\t")

sarcoma <- sarcoma[,c("V1", "V2", "V6")]
table(sarcoma$V6)

isks <- isks[,c("V1", "V2", "V6")]
table(isks$V6)

mgrb <- mgrb[,c("V1", "V2", "V6")]
table(mgrb$V6)

bgi <- bgi[,c("V1", "V2", "V6")]
bgi$V6 <- 2
table(bgi$V6)

phe <- rbind(sarcoma, isks, mgrb, bgi)
names(phe) <- c("FID", "IID", "PHENO")
write.table(phe, "/g/data/wq2/users/sn0095/datasets/PHENO_all4.txt", row.names=F, col.names=T, quote=F, sep = "\t")


#######################################################################################################################

# Step 7: [PLINK] GWAS using plink

plink_path="/home/552/sn0095/Tools/plink2_linux_x86_64_20220814/plink2"
file1="/g/data/wq2/users/sn0095/datasets/Combined_all4_merged2"
out="/g/data/wq2/users/sn0095/GWAS_results/GWAS_all4"
pheno_file="/g/data/wq2/users/sn0095/datasets/PHENO_all4.txt"
cov_file="/g/data/wq2/users/sn0095/datasets/COV_all4_PCSEX_FINAL2.txt"
keepid="/g/data/wq2/users/sn0095/datasets/COV_all4_PCSEX_FINAL2.txt"

$plink_path --bfile $file1 --pheno $pheno_file --pheno-name PHENO --require-pheno PHENO --logistic hide-covar cols=+a1freqcc --ci 0.95 --keep $keepid --maf 0.01  --chr 1-22 --out $out --covar $cov_file --covar-name PC1, PC2, PC3, PC4, PC5, Sex --covar-variance-standardize

#######################################################################################################################

# Step 8: [R] GWAS manhattan plot after running GWAS in plink
library(qqman)
library(data.table)

gwas_data <- fread("/g/data/wq2/users/sn0095/GWAS_results/GWAS_all4.PHENO.glm.logistic.hybrid", header = T)
names(gwas_data)[1] <- "CHR"
exc <- fread("/g/data/wq2/users/sn0095/datasets/EXCLUDE_198vars_8sd.txt", header  =T)
gwas_data <- gwas_data[!gwas_data$ID %in% exc$ID,]
dim(gwas_data)


png("/g/data/wq2/users/sn0095/GWAS_results/GWAS_all4_DISCOVERY.png", res = 400, width = 8, height = 5, units = "in")
manhattan(gwas_data,
          chr = "CHR",
          bp="POS",
          p ="P",
          snp = "ID",
          col = c("skyblue4", "grey"), suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), 
          annotatePval = 5e-08, annotateTop = T, check.overlap=T)
dev.off()


png("GWAS_all4_QQ.png", res = 400, width = 8, height = 5, units = "in")
qq(gwas_data$P)
dev.off()

#######################################################################################################################

# Step 9: PRS 

# 3-fold cross validation
# Split data into 70% GWAS, 30% PRS

# Path to PLINK software on your server
plink_path='/storage/coda1/p-ggibson3/0/snagpal3/rich_project_bio-gibson/Tools/plink2_Mar19'

# PLink format (bed/bim/fam)
data='path_to_genotype_data_in_PLINKformat'

# 1. From GWAS summary statistics, select variants with p < 1e-03
# 2. Extract SNPs from your genotype data using --extract function in PLINK
# 3. LD pruning to obtain independent SNPs: 
${plink_path} --bfile ${data} --indep-pairwise 1000kb 1 0.2 --maf 0.01 --out ${out_file}

# 4. Using the pruned SNP list, make profile file usng GWAS summary stats in text format : 
# column headers: SNP, Effect Allele, Effect Size
profile='path_to_profile.txt'

# Subset European White British Individual IDs 
wb_ids='/storage/home/hcoda1/3/snagpal3/scratch/White_British_FIDs_UKBB.txt'

# 5. PRS: run score using PLINK
${plink_march} --bfile ${data} --keep $wb_ids --score ${profile} cols=+scoresums list-variants --out ${out_file}


# Once you have PRS, plot decile/percentile of PRS against - 
#1. prevalence of sarcoma
#2. proportion of cases which are carriers vs non-carriers.
# Note carriers identified based on gene list (ACMG + sarcoma) from table S2 of Ballinger et al. Science (2023) 

#######################################################################################################################

# Step 10: DEPICT for LD clumping GWAS SNPs, near genes for variants, gene set encrihcment analysis and tissue set enrichment analysis

# Install DEPICT http://www.broadinstitute.org/mpg/depict/depict_download/bundles/DEPICT_v1_rel173.tar.gz
# Get GWAS summary statistics in the format required for DEPICT
cd ~/Dropbox\ \(GaTech\)/Postdoc/Garvan/DEPICT/src/python 
conda create -n python2 python=2.7 anaconda
source activate python2

#installed python2.7, pandas numpy, sudo pip install interval tree
#plink1.9
#java installed
#cytoscape install

# GSEA analyses on gwas summary statistics
# Add paths to gwas summary statistics and plink in the template.cfg file 
python depict.py template.cfg

#######################################################################################################################
