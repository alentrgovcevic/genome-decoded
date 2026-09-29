# HLA imputation with HIBAG from the raw (not imputed) array genotypes, using the pre-fit
# Illumina Global Screening Array v2.0 European model (hg19, two-field resolution).
# Run (from the run directory): Rscript scripts/hla.R [model.RData]   → results/hla/hla_calls.tsv
# The model file lives in reference/hibag/ (download_refs.sh); config: hibag_model.
suppressMessages(library(HIBAG))

args <- commandArgs(trailingOnly = TRUE)
model_file <- file.path("reference/hibag", if (length(args)) args[1] else "InfiniumGlobal-European-HLA4-hg19.RData")
out_dir <- "results/hla"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

geno <- hlaBED2Geno("data/hla/mhc.bed", "data/hla/mhc.fam", "data/hla/mhc.bim", assembly = "hg19")
models <- get(load(model_file))

rows <- list()
for (locus in names(models)) {
  model <- hlaModelFromObj(models[[locus]])
  n_model <- length(model$snp.id)
  n_found <- sum(model$snp.position %in% geno$snp.position)
  pred <- hlaPredict(model, geno, type = "response+prob", match.type = "Position", verbose = FALSE)
  v <- pred$value
  rows[[locus]] <- data.frame(locus = locus, allele1 = v$allele1, allele2 = v$allele2, prob = v$prob,
                              matching = v$matching, model_snps = n_model, snps_found = n_found)
  hlaClose(model)
}
res <- do.call(rbind, rows)
write.table(res, file.path(out_dir, "hla_calls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE)
