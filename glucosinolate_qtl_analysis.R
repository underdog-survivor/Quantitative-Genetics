# ==============================================================================
# QTL Analysis of Glucosinolate Content in Rapeseed (Brassica napus)
# Immortalized F2 Population 
# ==============================================================================

# ----------------------------
# 1. Load Packages & Data
# ----------------------------
library(qtl)

# Load genotype/phenotype data (Use generic filenames for public repos)
mxt <- read.cross(
  format = "csvs",
  dir = "data",
  genfile = "genotypes.csv",
  phefile = "phenotypes.csv",
  estimate.map = FALSE
)

# ----------------------------
# 2. Preprocessing
# ----------------------------
mxt <- est.rf(mxt)  # Estimate recombination frequencies
mxt <- formLinkageGroups(mxt, max.rf = 0.35, reorgMarkers = TRUE)
mxt <- sim.geno(mxt)  # Impute missing genotypes

cat("Genetic map structure:\n")
print(pull.map(mxt))

# ==============================================================================
# 3. Analysis Pipeline
# ==============================================================================

# ----------------------------
# (a) Genome Scan & Threshold
# ----------------------------
# Single QTL scan
s1 <- scanone(mxt, pheno.col = "glu", method = "mr")

# Permutation test (1000 replicates)
s2 <- scanone(mxt, pheno.col = "glu", method = "mr", n.perm = 1000)

# Visualization
png("results/glu_qtl_scan.png", width = 9.5, height = 4, units = "in", res = 300)
par(mfrow = c(2, 1), mar = c(4, 4, 0.5, 0.5))

# Chromosome 1
plot(s1, chr = 1, show.marker.names = FALSE, 
     col = "blue", lwd = 2, main = "Chromosome 1 QTL Scan")
add.threshold(s1, perms = s2, alpha = 0.05, col = "red", lty = 2)

# Chromosomes 2-5
plot(s1, chr = 2:5, show.marker.names = FALSE,
     col = "blue", lwd = 2, main = "Other Chromosomes")
add.threshold(s1, perms = s2, alpha = 0.05, col = "red", lty = 2)
dev.off()

# ----------------------------
# (b) Significant QTL Peaks
# ----------------------------
lod_threshold <- summary(s2)[1]  # 5% threshold
significant_qtl <- summary(s1, threshold = lod_threshold)
cat("\nSignificant QTL peaks (LOD >", lod_threshold, "):\n")
print(significant_qtl)

# ----------------------------
# (c) Effect Estimation
# ----------------------------
# Define major QTL
major_qtl <- makeqtl(mxt, 
                     chr = c(1, 5), 
                     pos = c(3840, 0),
                     qtl.name = c("sS1702", "IGF0235b"))

# Fit full model
qtl_model <- fitqtl(mxt, pheno.col = "glu", qtl = major_qtl,
                   formula = y ~ Q1 + Q2 + Q1:Q2, get.ests = TRUE)

cat("\nAdditive, Dominance & Epistatic Effects:\n")
print(summary(qtl_model))

# ----------------------------
# (d) Epistasis Scan
# ----------------------------
# Two-dimensional genome scan
s3 <- scantwo(mxt, pheno.col = "glu", method = "mr")

# Visualization
png("results/glu_epistasis.png", width = 6, height = 5, units = "in", res = 300)
plot(s3, chr = c(1, 5), main = "Epistatic Interactions")
dev.off()

# Save results
save.image("qtl_analysis_results.RData")
sink("session_info.txt")
sessionInfo()
sink()
