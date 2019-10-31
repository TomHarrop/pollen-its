#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)

# parse taxonomy files 
GenerateQiimeTaxonomyLine <- function(
    x,
    allowed_taxa = allowed_taxa,
    ranks = ranks) {
    # remove blank entries
    all_names <- as.character(x[x != ''])
    
    # subset allowed_taxa by taxon names in x
    matched_taxa <- unique(allowed_taxa[x],
                           by = c("name_txt", "rank"))
    
    # make sure every rank in ranks is represented
    all_ranks <- merge(data.table(rank = ranks),
                       matched_taxa,
                       by = "rank",
                       all.x = TRUE)
    all_ranks[, rank := factor(rank, levels = ranks)]
    setorder(all_ranks, rank)
    
    # replace NA with blanks
    all_ranks[is.na(name_txt), name_txt := ""]
    
    # use qiime format for taxa names (prefix with node class letter)
    all_ranks[, name_qiime := paste(substr(rank, 1, 1), name_txt, sep = "_")]
    return(all_ranks[, paste0(name_qiime, collapse = ";")])
}

# load data
tax_names <- readRDS(snakemake@input[["names_dmp"]])
tax_nodes <- readRDS(snakemake@input[["nodes_dmp"]])
acc_tax <- fread(snakemake@input[["acc_to_taxa"]])
acc_to_taxline <- snakemake@output[["acc_to_taxline"]]

# dev
# tax_names <- readRDS("data/ncbi/names.dmp.Rds")
# tax_nodes <- readRDS("data/ncbi/nodes.dmp.Rds")
# # taxonomy from genbank rbcL search results
# acc_tax <- fread("output/taxonomy/genbank_results.txt")
# acc_to_taxline <- "output/taxonomy/taxonomy.txt"

# generate list of acceptable names
ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
tax_ids_to_keep <- tax_nodes[rank %in% ranks, unique(tax_id)]
filtered_tax_names <- tax_names[tax_id %in% tax_ids_to_keep &
                                    name_class == "scientific name"]
allowed_taxa <- merge(filtered_tax_names[, .(tax_id, name_txt)],
                      tax_nodes[, .(tax_id, rank)],
                      all.x = TRUE,
                      all.y = FALSE)
setkey(allowed_taxa, name_txt)
unique_taxa <- unique(allowed_taxa, by = c("name_txt", "rank"))

# generate a taxonomy line for each accession
taxo <- acc_tax[, .(
    taxonomy = GenerateQiimeTaxonomyLine(
        taxon, unique_taxa, ranks)),
    by = .(accession, species)]

# add species from acc_tax
taxo[, species := gsub("[^[:alnum:]]+", ".", species)]
taxo[, taxonomy := paste(taxonomy, species, sep = ";s_")]
taxo[, species := NULL]

# write output
fwrite(taxo,
       file = acc_to_taxline,
       quote = FALSE,
       sep = "\t",
       col.names = FALSE)

# log
sessionInfo()

