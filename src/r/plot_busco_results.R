#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# busco reading function
ReadBuscoResults <- function(filepath){
    my_busco_data <- fread(filepath,
                           header = TRUE,
                           skip = 4,
                           fill = TRUE,
                           na.strings = c(""))
    filename <- gsub("full_table_(.+).tsv", "\\1", basename(filepath))
    my_busco_data[, filename := filename]
    my_busco_data
}

# read filename mapping
assembly_filenames <- fread("data/assembly_filenames.csv")
assembly_filenames[, filename := gsub(".fn?a$", "", filename)]

# find busco results
tld <- list.dirs(path = "output/busco", recursive = FALSE)
tldp1 <- list.dirs(path = tld, recursive = FALSE)
search_dirs <- grep("run_", tldp1, value = TRUE)
busco_result_files <- list.files(path = search_dirs,
                                 pattern = "full_table_",
                                 full.names = TRUE,
                                 recursive = FALSE)

# read data
busco_result_list <- lapply(busco_result_files, ReadBuscoResults)
busco_results <- rbindlist(busco_result_list)

# add species and family info to busco results
busco_results_named <- merge(busco_results,
                             assembly_filenames,
                             by = "filename")

# count the number in each status
plot_data <- busco_results_named[, .(
    n_buscos = length(unique(`# Busco id`))),
    by = .(Status, species_name, family, filename)]
plot_data[, total_buscos := sum(n_buscos),
          by = .(species_name, family, filename)]
plot_data[, status_percent := n_buscos * 100 / total_buscos]

# separate my assemblies from others
plot_data[grep("^mh_", filename), category := "Our data"]
plot_data[grep("^asw_", filename), category := "Our data"]
plot_data[is.na(category), category := "Comparison"]
plot_data[, category := factor(category,
                               levels = c("Our data", "Comparison"))]

# order categories
status_order <- c("Complete", "Duplicated", "Fragmented", "Missing")
plot_data[, Status := factor(Status, levels = status_order)]

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


# visualise
Set1 <- RColorBrewer::brewer.pal(9, "Set1")
ggplot(plot_data[family == "Braconidae"],
       aes(x = filename, y = status_percent, fill = Status)) +
    theme(strip.background = element_blank(),
          strip.text = element_blank()) +
    facet_grid( ~ category, scales = "free_x", space = "free_x") +
    xlab(NULL) + ylab(NULL) +
    scale_x_discrete(breaks = filename_order, labels = fn_lab) +
    scale_fill_manual(values = Set1,
                      guide = guide_legend(title = NULL)) +
    geom_col(width = 0.75,
             position = position_dodge(width = 0.75))

