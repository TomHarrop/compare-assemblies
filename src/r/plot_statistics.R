#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# read data
raw_data <- fread("output/assembly_statistics/stats.txt")
assembly_filenames <- fread("data/assembly_filenames.csv")
assembly_filenames[, filename := gsub(".fn?a$", "", filename)]

# tidy
setnames(raw_data, "filename", "filepath")
raw_data[, filename := gsub(".fn?a$", "", basename(filepath))]

# manually edit weird entries
raw_data[grep("mh-nopcr/output/meraculous/bbduk/run_41mer", filepath),
         filename := "mh_diploid1_41mer"]
raw_data[grep("asw-thruplex/output/soap_denovo2/decon/run_51mer", filepath),
         filename := "asw_thruplex_soap_51mer"]
raw_data[grep("asw-nopcr/output/meraculous_diploid2/bbduk/run_51mer", filepath),
         filename := "asw_diploid2_51mer"]
raw_data[grep("asw-nopcr/output/meraculous/bbduk/run_51mer", filepath),
         filename := "asw_diploid1_51mer"]

# add species and family info to results
assembly_stats_named <- merge(raw_data,
                              assembly_filenames,
                              by = "filename")
assembly_stats_named[, filepath := NULL]

# reformat some variables
assembly_stats_named[, `Scaffolds (thousands)` := n_scaffolds / 1000]
assembly_stats_named[, `Contigs (thousands)` := n_contigs / 1000]
assembly_stats_named[, `Contig length (MB)` := contig_bp / 1000000]
assembly_stats_named[, `Scaffold length (MB)` := scaf_bp / 1000000]
assembly_stats_named[, `Scaffold L50 (KB)` := scaf_L50 / 1000]
assembly_stats_named[, `Contig L50 (KB)` := ctg_L50 / 1000]
assembly_stats_named[, `GC (%)` := gc_avg * 100]
setnames(assembly_stats_named, "gap_pct", "Gaps (%)")

# go long
var_order <- c(
    "Contigs (thousands)",
    "Scaffolds (thousands)",
    "Contig length (MB)",
    "Scaffold length (MB)",
    "Contig L50 (KB)",
    "Scaffold L50 (KB)",
    "Gaps (%)",
    "GC (%)")
plot_data <- melt(assembly_stats_named,
                  id.vars = c("filename", "species_name", "family"),
                  measure.vars = var_order)
plot_data[, variable := factor(as.character(variable), levels = var_order)]

# add category
plot_data[grep("^mh_", filename), category := "Our data"]
plot_data[grep("^asw_", filename), category := "Our data"]
plot_data[is.na(category), category := "Comparison"]

# make labels
filename_order <- c("asw_thruplex_soap_51mer",
                    "asw_diploid2_51mer",
                    "asw_diploid1_51mer", 
                    "mh_diploid1_41mer",
                    "GCA_000956155.1_ASM95615v1_genomic",
                    "GCF_000355655.1_DendPond_male_1.0_genomic", 
                    "GCF_001412515.1_Dall1.0_genomic",
                    "GCF_000806365.1_ASM80636v1_genomic", 
                    "GCA_001012855.1_ASM101285v1_genomic")
fn_lab <- plyr::revalue(filename_order, c(
    "asw_thruplex_soap_51mer" = bquote(italic("Listronotus bonariensis")*" ThruPlex"),
    "asw_diploid1_51mer" = bquote(italic("Listronotus bonariensis")*" phased"),
    "asw_diploid2_51mer" = bquote(italic("Listronotus bonariensis")*" unphased"),
    "mh_diploid1_41mer" =  bquote(italic("Microctonus hyperodae")*" phased"),
    "GCA_000956155.1_ASM95615v1_genomic" = bquote(italic("Cotesia vestalis")),
    "GCF_000355655.1_DendPond_male_1.0_genomic" = bquote(italic("Dendroctonus ponderosae")), 
    "GCF_001412515.1_Dall1.0_genomic" = bquote(italic("Diachasma alloeum")),
    "GCF_000806365.1_ASM80636v1_genomic" = bquote(italic("Fopius arisanus")), 
    "GCA_001012855.1_ASM101285v1_genomic" = bquote(italic("Hypothenemus hampei"))
))
plot_data[, filename := factor(filename, levels = filename_order)]

# save plot data
saveRDS(plot_data, "output/plots/stats_pd.Rds")

# vis

