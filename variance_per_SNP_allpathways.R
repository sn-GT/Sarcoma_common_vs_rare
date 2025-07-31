setwd("~/Dropbox (GaTech)/Postdoc/Garvan/Analysis/Replication/BGI/GWAS_withBGI_nocohort/Final/")
library(data.table)
library(dplyr)
library(ggplot2)

# GWAS sarcoma for combined dataset
mdf <- fread("GWAS_all4.PHENO.glm.logistic.hybrid", header = TRUE)
mdf$ID2 <- paste0(mdf$`#CHROM`, ":", mdf$POS)
mdf$beta <- log(mdf$OR)

# Liability scale conversion
K <- 7e-5
t <- qnorm(1 - K)
z <- dnorm(t)
scale_factor <- K * (1 - K) / z
mdf$beta <- mdf$beta * scale_factor

# Load allele frequencies
af <- fread("combined_all4_merged2_freq.afreq", header = TRUE, sep = "\t")
af$CHR <- sapply(strsplit(af$ID, ":"), `[`, 1)
af$POS <- sapply(strsplit(af$ID, ":"), `[`, 2)
af$ID2 <- paste0(af$CHR, ":", af$POS)
af$ALT_FREQS <- as.numeric(af$ALT_FREQS)
colnames(af)[which(names(af) == "ALT")] <- "ALT_af"

# Merge allele freq with GWAS
merged <- merge(mdf, af[, c("ID2", "ALT_af", "ALT_FREQS")], by = "ID2", all.x = TRUE)
merged$AF <- ifelse(merged$A1 == merged$ALT_af, merged$ALT_FREQS, 1 - merged$ALT_FREQS)
merged$AF <- as.numeric(merged$AF)
merged$var <- 2 * merged$AF * (1 - merged$AF) * merged$beta^2

# Load gene-to-SNP mapping
genes <- read.csv('/Users/sininagpal/GaTech Dropbox/Sini Nagpal/Postdoc/Garvan/MS/SupplementaryTables/SupplementaryTable2_GWAS_topHits.csv', header = TRUE)

# list of genes in each pathway
nedd8_genes <- c("RNF6", "RNF123", "KLH20", "TRAIP", "STUB1", "KCTD6", "DZIP3", "KLH11",
                 "ZBTB16", "UFL1", "UBE2F", "UBA7", "ARIH2", "TRIM41", "FBXL16", "TRIM4", 
                 "UBE2G1", "ASB8", "ATG7", "UBE2R2", "ASB2", "KLHL3", "RNF116")

mrpl_genes <- c("MRPS18A", "MRPS14", "RPL3L", "MRPL32", "MRPL18", "GFM2", "MTIF3", "MRPL44",
                "MRPL46", "MRPL39", "MRPS2", "MRPS5", "MRPS15", "POLR1C", "NSA2")

glu_genes <- c("NETO1", "DLG2")

dnadamage_tp53_genes <- c("ATR", "EEF1E1", "EPRS1", "HIC1")

smooth_muscle_genes <- c("ESR1", "FOXO1", "MAP2K5", "SOD2", "STUB1", "MYOCD", "RGS5", "TAFA5")

neural_genes <- c("BAIAP2", "CAMKV", "CTNNA2", "ITSN1", "MRTFA", "RHOA", "APP", "CNR1",
                  "CNTN6", "CNTNAP2", "NTM", "ROBO1", "ROBO2", "SEMA5A", "ALCAM", "BMPR1B",
                  "EFNA5", "PTPRM")

# Extract pathway SNPs
get_pathway_vars <- function(pathway_genes, label, merged_df, genes_df) {
  pathway_snps <- unique(genes_df$ID2[genes_df$Gene.symbol %in% pathway_genes])
  df <- merged_df[merged_df$ID2 %in% pathway_snps, ]
  df$Pathway <- label
  return(df)
}

# Get per-pathway SNP data
mrpl_df     <- get_pathway_vars(mrpl_genes,     "MRPL",     merged, genes)
nedd8_df    <- get_pathway_vars(nedd8_genes,    "NEDD8",    merged, genes)
glu_df      <- get_pathway_vars(glu_genes,      "Glutamate",merged, genes)
tp53_df     <- get_pathway_vars(dnadamage_tp53_genes, "TP53_DNA_Damage", merged, genes)
smooth_df   <- get_pathway_vars(smooth_muscle_genes, "SmoothMuscleApoptosis", merged, genes)
neural_df   <- get_pathway_vars(neural_genes,   "Neural",   merged, genes)

combined_pathway_df <- rbind(mrpl_df, nedd8_df, glu_df, tp53_df, smooth_df, neural_df)

# AF-matched null sampling for each SNP
set.seed(42)
null_matched_df <- do.call(rbind, lapply(1:nrow(combined_pathway_df), function(i) {
  af_i <- combined_pathway_df$AF[i]
  snp_id <- combined_pathway_df$ID2[i]
  candidates <- merged %>% filter(abs(AF - af_i) <= 0.04 & !is.na(var))
  if (nrow(candidates) >= 1000) {
    sampled_vars <- replicate(1000, {
      sampled_snp <- candidates[sample(nrow(candidates), 1), ]
      sampled_snp$var
    })
    data.frame(ID2 = snp_id, Pathway = "Null", var = mean(sampled_vars, na.rm = TRUE))
  } else {
    data.frame(ID2 = snp_id, Pathway = "Null", var = NA)
  }
}))
null_matched_df <- null_matched_df[!is.na(null_matched_df$var), ]

# plotting dataframe
plot_df <- combined_pathway_df[, c("ID2", "Pathway", "var")]
colnames(plot_df)[3] <- "var"
plot_df <- rbind(plot_df, null_matched_df)


# Step 1: Sort by decreasing mean variance
mean_var <- plot_df %>%
  group_by(Pathway) %>%
  summarise(MeanVar = mean(var, na.rm = TRUE)) %>%
  arrange(desc(MeanVar))

plot_df$Pathway <- factor(plot_df$Pathway, levels = mean_var$Pathway)

# Step 2: Set label mapping
label_map <- c(
  "MRPL" = "MRPL complex",
  "NEDD8" = "NEDD8 Pathway",
  "Glutamate" = "Glutamate signalling",
  "TP53_DNA_Damage" = "DNA damage, TP53",
  "SmoothMuscleApoptosis" = "Smooth Muscle Apoptosis",
  "Neural" = "Neural Axon Guidance",
  "Null" = "Null (AF-matched)"
)

# Plot variance per SNP for each pathway 
pf <- ggplot(plot_df, aes(x = Pathway, y = var * 100, fill = Pathway)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_x_discrete(labels = label_map[levels(plot_df$Pathway)]) +
  theme_bw() +
  labs(y = "Variance per SNP (2pqβ²) (%)") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10, color = "black", angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(size = 10)
  )
pf


# summary statistics for each pathway
summary_stats <- plot_df %>%
  group_by(Pathway) %>%
  summarise(
    MeanVariance = mean(var, na.rm = TRUE) * 100,
    SDVariance = sd(var, na.rm = TRUE),
    N = n()
  )

print(summary_stats)

png("../../../../../MS/code_new/pathways_var.png", width = 5.5, height = 4, units = "in", res=200)
print({pf})
dev.off()


# Rename pathway gene sets
pathway_lists <- list(
  MRPL = mrpl_genes,
  NEDD8 = nedd8_genes,
  Glutamate = glu_genes,
  TP53_DNA_Damage = dnadamage_tp53_genes,
  SmoothMuscleApoptosis = smooth_muscle_genes,
  Neural = neural_genes
)

# Filter genes per pathway to those present in the GWAS gene list
filtered_pathway_genes <- lapply(names(pathway_lists), function(pathway) {
  gene_set <- pathway_lists[[pathway]]
  matched_genes <- unique(genes$Gene.symbol[genes$Gene.symbol %in% gene_set])
  if (length(matched_genes) > 0) {
    return(data.frame(Pathway = pathway, Gene = matched_genes))
  } else {
    return(NULL)  
  }
})

pathway_gene_table <- do.call(rbind, filtered_pathway_genes)

