#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)

# download taxdump (complete NCBI taxonomy)
taxdmp_url <- "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"
temp2 <- tempfile(fileext = ".zip")
temp3 <- tempdir()
download.file(taxdmp_url, temp2)

# read nodes.dmp from inside taxdmp
nodes_dmp_file <- unzip(temp2, files = "nodes.dmp", exdir = temp3)
nodes_dmp_raw <- fread(nodes_dmp_file, sep = "|", quote = "\t", header = FALSE)

# tidy up nodes.dmp
nodes_dmp_raw[, V14 := NULL]
nodes_names <- c("tax_id", "parent_tax_id", "rank", "embl_code", "division_id",
                 "inherited_div_flag", "genetic_code_id", "inherited_GC",
                 "mitochondrial_genetic_code_id", "inherited_MGC_flag",
                 "GenBank_hidden_flag", "hidden_subtree_root_flag", "comments")
setnames(nodes_dmp_raw, names(nodes_dmp_raw), nodes_names)
nodes_dmp_raw[, tax_id := as.numeric(gsub("\t", "", tax_id, fixed = TRUE))]
setkey(nodes_dmp_raw, tax_id)

# read names.dmp from inside taxdmp
names_dmp_file <- unzip(temp2, files = "names.dmp", exdir = temp3)
names_dmp_raw <- fread(names_dmp_file, sep = "|", quote = "\t", header = FALSE)

# tidy up names.dmp
names_dmp_raw[, V5 := NULL]
names_names <- c("tax_id", "name_txt", "unique_name", "name_class")
setnames(names_dmp_raw, names(names_dmp_raw), names_names)
names_dmp_raw[, tax_id := as.numeric(gsub("\t", "", tax_id, fixed = TRUE))]
setkey(names_dmp_raw, tax_id)

# save output
saveRDS(nodes_dmp_raw, snakemake@output[['nodes_dmp']])
saveRDS(names_dmp_raw, snakemake@output[['names_dmp']])

# log
sessionInfo()
