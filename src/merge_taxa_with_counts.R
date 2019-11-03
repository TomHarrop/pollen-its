#!/usr/bin/env Rscript

library(data.table)

seqtab_file <- 'output/030_dada2/seqtab_nochim.Rds'
taxa_file <- "output/030_dada2/taxa.Rds"

seqtab <- data.table(readRDS(seqtab_file), keep.rownames = TRUE)
taxa <- data.table(readRDS(taxa_file), keep.rownames = TRUE)

# assign ASV ids in the taxa results
taxa[, asv_id := paste0("asv", 1:.N)]

# map asv ids to sequence
asv_map <- taxa[, structure(asv_id, names = rn)]

# transpose seqtab
seqtab_t <- data.table(t(data.frame(seqtab, row.names = "rn")), keep.rownames = TRUE)

# merge seqtab and taxa
taxa_with_counts <- merge(taxa, seqtab_t, by = "rn")
setnames(taxa_with_counts, "rn", "seq")

col_order <- c("asv_id",
	"Kingdom",
	"Phylum",
	"Class",
	"Order",
	"Family",
	"Genus",
	"Species",
	paste0("s", 1:80),
	"seq")
setcolorder(taxa_with_counts, col_order)

taxa_with_counts[, asv_id := factor(asv_id, levels = gtools::mixedsort(asv_id))]
setkey(taxa_with_counts, asv_id)

fwrite(taxa_with_counts, "test/taxa_with_counts.csv")

