#!/usr/bin/env Rscript
# Import per-sample Salmon or kallisto transcript estimates without count scaling.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript import_tximport.R salmon|kallisto design.csv tx2gene.tsv quant_dir output_dir")
}

method <- args[1]
if (!method %in% c("salmon", "kallisto")) stop("Method must be salmon or kallisto")
if (!requireNamespace("tximport", quietly = TRUE)) stop("R package tximport is required")

design <- read.csv(args[2], stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(names(design), c("sample", "group"))) stop("Expected design.csv columns: sample,group")
if (anyDuplicated(design$sample) || anyNA(design$sample)) stop("Sample names must be unique and nonmissing")
txmap <- read.delim(args[3], stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(names(txmap), c("transcript_id", "gene_id", "gene_name"))) {
  stop("Expected tx2gene.tsv columns: transcript_id, gene_id, gene_name")
}
if (anyDuplicated(txmap$transcript_id)) stop("Duplicate transcript IDs in tx2gene.tsv")

filename <- if (method == "salmon") "quant.sf" else "abundance.tsv"
files <- file.path(args[4], design$sample, filename)
names(files) <- design$sample
if (!all(file.exists(files))) stop("Missing quantification: ", paste(files[!file.exists(files)], collapse = ", "))

txi <- tximport::tximport(files, type = method, txOut = TRUE,
                         countsFromAbundance = "no", dropInfReps = TRUE)
if (!all(rownames(txi$counts) %in% txmap$transcript_id)) {
  stop("Quantification contains transcript IDs absent from the GTF")
}
gene <- tximport::summarizeToGene(txi, txmap[, c("transcript_id", "gene_id")],
                                 countsFromAbundance = "no")

dir.create(args[5], recursive = TRUE, showWarnings = FALSE)
write_matrix <- function(values, ids, path, is_tx = FALSE) {
  key <- if (is_tx) "transcript_id" else "gene_id"
  label <- if (is_tx) {
    txmap$gene_name[match(ids, txmap$transcript_id)]
  } else {
    unique(txmap[, c("gene_id", "gene_name")])$gene_name[
      match(ids, unique(txmap[, c("gene_id", "gene_name")])$gene_id)]
  }
  if (anyNA(label)) stop("Missing gene names for output IDs")
  table <- data.frame(ids, label, values, check.names = FALSE)
  names(table)[1:2] <- c(key, "gene_name")
  write.table(table, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}

write_matrix(txi$counts, rownames(txi$counts),
             file.path(args[5], paste0(method, ".transcript.expected.tsv")), TRUE)
write_matrix(gene$counts, rownames(gene$counts),
             file.path(args[5], paste0(method, ".gene.expected.tsv")))
write_matrix(txi$abundance, rownames(txi$abundance),
             file.path(args[5], paste0(method, ".transcript.tpm.tsv")), TRUE)
write_matrix(gene$abundance, rownames(gene$abundance),
             file.path(args[5], paste0(method, ".gene.tpm.tsv")))

# Separate, explicitly labelled length-scaled counts for downstream uses.
scaled <- tximport::tximport(files, type = method, tx2gene = txmap[, 1:2],
                            countsFromAbundance = "lengthScaledTPM", dropInfReps = TRUE)
write_matrix(scaled$counts, rownames(scaled$counts),
             file.path(args[5], paste0(method, ".gene.lengthScaledTPM.tsv")))
cat(method, ":", nrow(txi$counts), "transcripts,", nrow(gene$counts), "genes\n")
