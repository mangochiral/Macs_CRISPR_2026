## Set up folders

dirs <- c("bam", "refs", "results", "plots")
for (d in dirs) {
  dir.create(d, showWarnings = FALSE)
}


## Sample sheet was generated manually 'samples.csv'

##Upload sample sheet

library(data.table)

samp <- fread("samples.csv")
# trim whitespace in case there are spaces after commas
samp[, `:=`(
  Sample = trimws(Sample),
  condition = trimws(condition),
  donor = trimws(donor),
  bam = trimws(bam)
)]

stopifnot(all(file.exists(samp$bam)))

# Check balanced design (each donor should have each condition)
print(table(samp$donor, samp$condition))

## Download GTF file 

## Run the following in terminal 
#cd ~/CRISPRoff_RNAseq/CRISPRoff_RNAseq/refs
#curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz


## Run featurecounts

library(Rsubread)
library(data.table)

# point to your GTF (adjust filename if needed)
gtf <- "refs/gencode.v44.annotation.gtf"
stopifnot(file.exists(gtf))

fc <- featureCounts(
  files = samp$bam,
  annot.ext = gtf,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  allowMultiOverlap = TRUE,
  nthreads = max(1, parallel::detectCores() - 1)
)

count_mat <- fc$counts
colnames(count_mat) <- samp$Sample

# QC: assignment rates --> How many 
# fc$stat is a data.frame: rows = status categories, cols = BAM files
stat <- fc$stat

# pull assigned counts per BAM
assigned <- as.numeric(stat[stat$Status == "Assigned", -1])

# total alignments per BAM = sum across all status rows
total <- colSums(stat[, -1])

qc <- data.frame(
  Sample = samp$Sample,
  assigned = assigned,
  total = total,
  pct_assigned = round(100 * assigned / total, 2)
)

qc <- qc[order(qc$pct_assigned), ]
print(qc)

write.csv(qc, file = "results/featureCounts_QC_assignment_rates.csv", row.names = FALSE)


## Run DESEQ2

library(DESeq2)
library(apeglm)

# Make sure count matrix columns match your sample order
stopifnot(all(colnames(count_mat) == samp$Sample))

coldata <- data.frame(
  row.names = samp$Sample,
  condition = factor(samp$condition),
  donor = factor(samp$donor)
)

dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData   = coldata,
  design    = ~ donor + condition
)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Set reference condition (change if you want a different baseline)
dds$condition <- relevel(dds$condition, ref = "NT")

dds <- DESeq(dds)

# confirm levels
levels(dds$condition)


## Plot PCA 

dir.create("plots", showWarnings = FALSE)

vsd <- vst(dds, blind = FALSE)

pdf("plots/PCA_condition.pdf")
print(plotPCA(vsd, intgroup = "condition"))
dev.off()

pdf("plots/PCA_donor.pdf")
print(plotPCA(vsd, intgroup = "donor"))
dev.off()


## RUN DE, each condition vs. NT

dir.create("results", showWarnings = FALSE)

coef_names <- resultsNames(dds)
cond_coefs <- grep("^condition_", coef_names, value = TRUE)  # all condition effects vs ref (NT)

for (cf in cond_coefs) {
  res_sh <- lfcShrink(dds, coef = cf, type = "apeglm")
  
  out <- as.data.frame(res_sh)
  out$gene_id <- rownames(out)
  out <- out[order(out$padj), ]
  
  comp <- sub("^condition_", "", cf)  # e.g. "CD81_vs_NT"
  write.csv(out,
            file = sprintf("results/DESeq2_%s_shrunk.csv", comp),
            row.names = FALSE)
}

list.files("results", pattern = "^DESeq2_.*_shrunk\\.csv$")




# MA Plot CD81 vs NT
dir.create("plots", showWarnings = FALSE)

res_cd81_sh <- lfcShrink(dds, coef = "condition_CD81_vs_NT", type = "apeglm")

pdf("plots/MA_CD81_vs_NT.pdf")
plotMA(res_cd81_sh, ylim = c(-5, 5))
dev.off()

## MA plot CD45g1 vs NT



res_CD45g1_sh <- lfcShrink(dds, coef = "condition_CD45g1_vs_NT", type = "apeglm")

pdf("plots/MA_CD45g1_vs_NT.pdf")
plotMA(res_CD45g1_sh, ylim = c(-5, 5))
dev.off()
