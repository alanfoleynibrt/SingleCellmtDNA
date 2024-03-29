---
title: "Single-cell mtDNA sequencing methods paper"
author: "Alan Foley alan_foley@live.co.uk"
date: "30/09/2023"
output:
  html_document:
    number_sections: true
    toc: true
    toc_depth_depth: 2
    toc_float: true
    theme: readable
  html_notebook: default
  pdf_document: default
---

# Overview

This R notebook enables reproduction of the R downstream analysis in: (LINK TO PAPER)

Mitochondrial DNA (mtDNA) from 4 single Chinese hamster ovary (CHO) cells and 1 bulk CHO sample are sequenced and preprocessed in linux as per (https://github.com/alanfoleynibrt/SingleCellmtDNA). The KX576660.1 mtDNA reference genome is used for mapping.

# Raw Data Availability

All raw files and intermediate mutation files available in github repository https://github.com/alanfoleynibrt/SingleCellmtDNA

# Prerequisites

```{r setup, include=FALSE}
knitr::opts_chunk$set(warning = FALSE, message = FALSE)
```

```{r setup1, warning=FALSE}
# library
library(forcats)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library("tidyverse")
library(scales)
library(xlsx)
library(colorRamps)
library(ggplot2)
library(stringr)

# theme

theme_ms <- function(base_size = 10, base_family = "Helvetica") {
  library(grid)
  (theme_bw(base_size = base_size, base_family = base_family) +
      theme(
        text = element_text(color = "black"),
        axis.title = element_text(face = "bold", size = rel(1)),
        plot.title = element_text(face = "bold", size = rel(1)),
        axis.text = element_text(size = rel(1), color = "black"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.8, "lines"),
        panel.border = element_rect(color = "black", size = 1),
        panel.grid = element_blank()
      ))
}

# create a directory for saving figures
figsdir <- "./figures/"
if (!dir.exists(figsdir)) {
  dir.create(figsdir)
}

tabsdir <- "./tables/"
if (!dir.exists(tabsdir)) {
  dir.create(tabsdir)
}

# function to save figures
saveManuscriptPlot <- function(p, width, height) {
  figfile <- file.path(figsdir, sprintf(
    "%s.png",
    gsub("\\.", "_", deparse(substitute(p)))
  ))
  ggsave(figfile, plot = p, width = width, height = height, units = "in", device = "png")
  
  fill.gradient <- scale_fill_gradientn(
    colors = matlab.like(200),
    trans = "log10"
  )
}
```

# Read mapping to KX576660.1 CHO mtDNA reference genome

```{r setup2, warning=FALSE}
# read mapping statistics

mapping_data_unshifted <- read.table(file = "plotting_data/mapping_statistics/mapping_rates.txt")

mapping_data_unshifted <- mapping_data_unshifted %>%
  mutate(
    Filename = V1,
    Unmapped = V2 - V3,
    Duplicates = V4,
    Unique = V3 - V4
  ) %>%
  mutate(sample_name = case_when(
    Filename == "SC-H-K_S4" ~ "Single Cell 1",
    Filename == "SC-H-M_S7" ~ "Single Cell 2",
    Filename == "SC-H-O_S8" ~ "Single Cell 3",
    Filename == "SC-H-P_S9" ~ "Single Cell 4",
    Filename == "MP-H-R_S5" ~ "Bulk"
  )) %>%
  select(sample_name, Unmapped, Duplicates, Unique)

write.xlsx(x = mapping_data_unshifted, 
           sheetName = "Read Mapping", 
           file = "tables/Supplementary Data.xlsx",
           append = F )

mapping_data_unshifted <- mapping_data_unshifted %>%
  pivot_longer(
    cols = c("Unmapped", "Duplicates", "Unique"),
    names_to = "Variable",
    values_to = "Reads"
  )

mapping_data_unshifted$Variable <- fct_relevel(factor(mapping_data_unshifted$Variable), c("Unmapped", "Duplicates", "Unique"))

scientific <- function(x) {
  ifelse(x == 0, "0", parse(text = gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

Fig1_read_mapping <- mapping_data_unshifted %>%
  ggplot(aes(x = sample_name, y = Reads, fill = Variable)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_ms() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 18)) +
  labs(fill = "", y = "# Reads", x = "", title = "Read mapping to the Chinese hamster mtDNA sequence") +
  theme(legend.position = "right", text = element_text(size = 18)) +
  scale_y_continuous(
    expand = c(0, 10000),
    labels = scientific
  ) +
  scale_fill_brewer(palette = "Set1")


saveManuscriptPlot(Fig1_read_mapping, 10, 5)
Fig1_read_mapping
```

# Average Sample Coverage

```{r setup3, warning=FALSE}

average_cell_line_coverage <- read.table("plotting_data/mapping_statistics/coverage/sample_depth", header = F)

average_cell_line_coverage <- average_cell_line_coverage %>%
  mutate(
    Filename = as.factor(V1),
    average_coverage = V2
  ) %>%
  mutate(sample_name = case_when(
    Filename == "SC-H-K_S4" ~ "Single Cell 1",
    Filename == "SC-H-M_S7" ~ "Single Cell 2",
    Filename == "SC-H-O_S8" ~ "Single Cell 3",
    Filename == "SC-H-P_S9" ~ "Single Cell 4",
    Filename == "MP-H-R_S5" ~ "Bulk"
  )) %>%
  select(sample_name, average_coverage)

write.xlsx(x = average_cell_line_coverage, 
           sheetName = "Average Sample Coverage", 
           file = "tables/Supplementary Data.xlsx",
           append = T )

Fig2_average_sample_coverage <- average_cell_line_coverage %>%
  ggplot(aes(x = sample_name, y = average_coverage)) +
  geom_bar(stat = "identity", fill = "#377EB8") +
  labs(x = "", y = "Depth of Coverage", title = "Average sequencing depth for each sample") +
  scale_y_continuous(expand = c(0, 100)) +
  theme_ms() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 18), text = element_text(size = 18))

saveManuscriptPlot(Fig2_average_sample_coverage, 10, 5)
Fig2_average_sample_coverage
```

# Perbase coverage

```{r setup4, warning=FALSE}
# perbase coverage graphs corrected around 0

perbase_cell_files_concat <- list.files("depthplusshifted2", full.names = T)

perbase_concat <- do.call("rbind", lapply(perbase_cell_files_concat, function(fn) {
  data.frame(Filename = fn, read.table(fn))
})) %>%
  select(Filename, V2, V3) %>%
  mutate(
    perbase = V3,
    position = V2
  ) %>%
  mutate(Filename = str_remove(Filename, "depthplusshifted2/")) %>%
  mutate(sample_name = case_when(
    Filename == "SC-H-K_S4" ~ "Single Cell 1",
    Filename == "SC-H-M_S7" ~ "Single Cell 2",
    Filename == "SC-H-O_S8" ~ "Single Cell 3",
    Filename == "SC-H-P_S9" ~ "Single Cell 4",
    Filename == "MP-H-R_S5" ~ "Bulk"
  )) %>%
  select(sample_name, position, perbase)

# average perbase of all samples

FigDepthConcatAverage <- perbase_concat %>%
  group_by(position) %>%
  summarise(mean_cov = mean(perbase)) %>%
  ggplot(aes(x = position, y = mean_cov)) +
  geom_line(color = "#377EB8") +
  labs(y = "Depth of Coverage", x = "mtDNA position (bp)", title = "Average perbase coverage for shifted mtDNA sequence") +
  theme_ms()
FigDepthConcatAverage
saveManuscriptPlot(FigDepthConcatAverage, 10,5)
```

```{r setup16, warning=FALSE}
# perbase of all samples individually

FigDepthConcatAll <- perbase_concat %>%
  ggplot(aes(x = position, y = perbase)) +
  geom_line(color = "#377EB8") +
  facet_wrap(~sample_name, ncol = 5) +
  labs(y = "Depth of Coverage", x = "mtDNA position (bp)", title = "Perbase coverage against KX576660 CHO mtDNA reference genome") +
  theme_ms() +
  theme(text = element_text(size = 8))
saveManuscriptPlot(FigDepthConcatAll, 20, 5)
FigDepthConcatAll
```

```{r setup99, warning=FALSE}
# violin plot of depths scDNAseq

Violin_Of_Depths <- perbase_concat %>%
group_by(sample_name) %>%
summarise(AverageSampleDepth = mean(perbase)) %>%
filter(!sample_name %in% c("W0", "W2", "NW3", "IW3", "NW3", "NW4"))

Violin_of_Depths <- ggplot(Violin_Of_Depths, aes(x="", y=AverageSampleDepth)) + 
    geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_hline(yintercept=2000, linetype="dashed", color = "green")
  ylab("Average Sample Depth") +
  theme(text = element_text(size = 25))
saveManuscriptPlot(Violin_of_Depths, 20,5)
Violin_of_Depths
```

# Variant Analysis
## Generating data frame
```{r setup5, warning=FALSE}
# generate variant data frame

variant_tables <- list.files("reorFINAL6/", full.names = T)

variants <- do.call("rbind", lapply(variant_tables, function(fn) {
  data.frame(Filename = fn, read.table(fn, sep = "\t", header = T))
})) %>%
  filter(POS != "POS") %>%
  mutate(
    variant_type = sapply(str_split(as.character(ANN), "\\|"), function(x) x[2]),
    impact = sapply(str_split(as.character(ANN), "\\|"), function(x) x[3])
  ) %>%
  select(-ANN) %>%
  mutate(POS = as.numeric(POS)) %>%
  mutate(gene = case_when(
    POS >= 1 & POS <= 69 ~ "tRNA-Phe",
    POS >= 71 & POS <= 1023 ~ "s-rRNA",
    POS >= 1024 & POS <= 1095 ~ "tRNA-Val",
    POS >= 1096 & POS <= 2657 ~ "l-rRNA",
    POS >= 2658 & POS <= 2732 ~ "tRNA-Leu_1",
    POS >= 2733 & POS <= 3687 ~ "ND1",
    POS >= 3688 & POS <= 3756 ~ "tRNA-Ile",
    POS >= 3754 & POS <= 3824 ~ "tRNA-Gln",
    POS >= 3829 & POS <= 3897 ~ "tRNA-Met",
    POS >= 3898 & POS <= 4930 ~ "ND2",
    POS >= 4931 & POS <= 4997 ~ "tRNA-Trp",
    POS >= 5000 & POS <= 5069 ~ "tRNA-Ala",
    POS >= 5072 & POS <= 5142 ~ "tRNA-Asn",
    POS >= 5174 & POS <= 5241 ~ "tRNA-Cys",
    POS >= 5241 & POS <= 5308 ~ "tRNA-Tyr",
    POS >= 5310 & POS <= 6854 ~ "COX1",
    POS >= 6852 & POS <= 6920 ~ "tRNA-Ser_1",
    POS >= 6924 & POS <= 6992 ~ "tRNA-Asp",
    POS >= 6994 & POS <= 7677 ~ "COX2",
    POS >= 7681 & POS <= 7744 ~ "tRNA-Lys",
    POS >= 7746 & POS <= 7949 ~ "ATP8",
    POS >= 7907 & POS <= 8587 ~ "ATP6",
    POS >= 8587 & POS <= 9370 ~ "COX3",
    POS >= 9371 & POS <= 9438 ~ "tRNA-Gly",
    POS >= 9439 & POS <= 9786 ~ "ND3",
    POS >= 9788 & POS <= 9855 ~ "tRNA-Arg",
    POS >= 9857 & POS <= 10153 ~ "ND4L",
    POS >= 10147 & POS <= 11524 ~ "ND4",
    POS >= 11525 & POS <= 11592 ~ "tRNA-His",
    POS >= 11593 & POS <= 11651 ~ "tRNA-Ser_2",
    POS >= 11651 & POS <= 11720 ~ "tRNA-Leu_2",
    POS >= 11721 & POS <= 13541 ~ "ND5",
    POS >= 13525 & POS <= 14049 ~ "ND6",
    POS >= 14050 & POS <= 14118 ~ "tRNA-Glu",
    POS >= 14123 & POS <= 15265 ~ "CYTB",
    POS >= 15267 & POS <= 15333 ~ "tRNA-Thr",
    POS >= 15337 & POS <= 15403 ~ "tRNA-Pro",
    POS >= 15404 & POS <= 16283 ~ "D-loop"
  )) %>%
  mutate(base_change = paste(REF, ALT, sep = ">")) %>%
  mutate(var_type = case_when(
    str_length(base_change) == 3 ~ "SNP",
    str_length(base_change) > 3 ~ "INDEL"
  )) %>%
  mutate(Filename = str_remove(Filename, "reorFINAL6/")) %>%
  mutate(Filename = str_remove(Filename, ".vcf")) %>%
  mutate(sample_name = case_when(
    Filename == "SC-H-K_S4" ~ "Single Cell 1",
    Filename == "SC-H-M_S7" ~ "Single Cell 2",
    Filename == "SC-H-O_S8" ~ "Single Cell 3",
    Filename == "SC-H-P_S9" ~ "Single Cell 4",
    Filename == "MP-H-R_S5" ~ "Bulk"
  )) %>%
  unite(mutation, c(POS, base_change), sep = " ", remove = FALSE) %>%
  filter(!sample_name %in% c("I1", "I8", "N1", "N32")) %>%
  unite(scFULL, c(sample_name, POS, base_change), sep = " ", remove = FALSE) %>%
  group_by(scFULL) %>%
  slice(which.max(DP))  

variants$sample_name <- factor(variants$sample_name, levels = sort(unique(variants$sample_name)))

variants2 <- variants

variants2[variants2=="intragenic_variant" | variants2=="intergenic_region"] <- "non_protein_coding"

variants2 <- filter(variants2, POS < 16276 & POS != 7968 & POS != 2948 & POS != 7961)

variants1 <- variants2

# filter only 4% heteroplasmic

variants2 <- filter(variants2, AF > 0.04, AF < 0.96)

# order of variants
variants2$sample_name <- factor(variants2$sample_name, levels=c("Bulk", "Single Cell 1", "Single Cell 2", "Single Cell 3", "Single Cell 4"))

```

## Variants per Sample

```{r setup6, warning=FALSE}

Fig7_variants_per_sample <- variants2 %>%
  group_by(sample_name) %>%
  count(var_type) %>%
  ggplot(aes(x = sample_name, y = n, fill = var_type)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "# Mutations", x = "", fill = "", title = "Mutations per sample") +
  theme(legend.position = "botton") +
  theme_ms() +
  scale_fill_brewer(palette = "Set1", ) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), text = element_text(size = 27))

saveManuscriptPlot(Fig7_variants_per_sample,20,10)
Fig7_variants_per_sample
```

## Variants per gene

```{r setup7, warning=FALSE}

Fig8_variants_per_sample <- variants2 %>%
  group_by(gene) %>%
  count(var_type) %>%
  ggplot(aes(x = gene, y = n, fill = var_type)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "# Mutations", x = "", fill = "",title = "Mutations per gene" ) +
  theme(legend.position = "botton") +
  theme_ms() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set1", ) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), text = element_text(size = 20))

saveManuscriptPlot(Fig8_variants_per_sample, 20, 10)
Fig8_variants_per_sample

```

# Allele frequency distribution

```{r setup8, warning=FALSE}
Fig9_allele_frequency_distribution <- variants2 %>%
  ggplot(aes(x = as.numeric(AF))) +
  geom_histogram(bins = 100, fill = "#377EB8", color = "black") +
  labs(x = "Allele Frequency", y = "# mtDNA mutations", 
       title = "Distribtion of all allele\nfrequencies across all samples",
       subtitle= "minimum allele frequency called > 4%") +
  theme_ms() +
  theme(text = element_text(size = 27))

saveManuscriptPlot(Fig9_allele_frequency_distribution, 10, 7)
Fig9_allele_frequency_distribution

```

# snpEff predicted impact of mutations
## per sample
```{r setup9, warning=FALSE}
# per sample
Fig10_predicted_impact_per_sample <- variants2 %>%
  select(sample_name, variant_type) %>%
  mutate(variant_type = case_when(
    variant_type == "frameshift_variant" ~ "frameshift",
    variant_type == "non_protein_coding" ~ "non protein coding",
    variant_type == "synonymous_variant" ~ "synonymous",
    variant_type == "missense_variant" ~ "missense",
    variant_type == "stop_gained" ~ "stop gained"
  )) %>%
  count(sample_name, variant_type) %>%
  ggplot(aes(x = sample_name, y = n, fill = variant_type)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "# Mutations", x = "", fill = "", title="snpEff predicted impact of mtDNA mutations per sample") +
  theme_ms() +
  scale_fill_brewer(palette = "Set1", na.value = "#F0E442") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,), text = element_text(size = 13)) 
saveManuscriptPlot(Fig10_predicted_impact_per_sample, 15, 5)
Fig10_predicted_impact_per_sample
```

## per gene

```{r setup17, warning=FALSE}
# per gene
Fig11_predicted_impact_per_gene <- variants2 %>%
  select(gene, variant_type) %>%
  count(gene, variant_type) %>%
  ggplot(aes(x = gene, y = n, fill = variant_type)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "# Variants", x = "", fill = "", title="snpEff predicted impact of mtDNA variants per gene") +
  theme_ms() +
  scale_y_continuous(expand = c(0, 1)) +
  scale_fill_brewer(palette = "Set1", na.value = "#F0E442") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "bottom")

saveManuscriptPlot(Fig11_predicted_impact_per_gene, 20, 10)
Fig11_predicted_impact_per_gene
```

# Variant heatmaps
## Base change
```{r setup10, warning=FALSE}
# base change heatmap
nb.cols <- 27
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)
Fig12_base_change_grid <- variants2 %>%
  select(sample_name, POS, base_change, AF) %>%
  ggplot(aes(x = as.factor(POS), y = as.factor(sample_name), fill = as.factor(base_change))) +
  geom_raster(hjust = 0, vjust = 0) +
  scale_fill_manual(values = mycolors) +
  theme(text = element_text(size = 20)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = -0.5, hjust = 1.25),
    axis.text.y = element_text(vjust = 2)
  ) +
  labs(x = "mtDNA position (nt)", y = "", fill = "Base \nChange",title="Mutation position base changes") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))
saveManuscriptPlot(Fig12_base_change_grid, 20,5)
Fig12_base_change_grid
```

## Allele frequency
```{r setup11, warning=FALSE}
# Allele frequency heatmap

fill.gradient <- scale_fill_gradientn(
    colors = matlab.like(200),
    trans = "log10"
  )

Fig13_allele_frequency_grid <- variants2 %>%
  select(sample_name, POS, base_change, AF) %>%
  ggplot(aes(x = as.factor(POS), y = as.factor(sample_name), fill = as.numeric(AF))) +
  scale_fill_gradient(low = "grey", high = "black") +
  geom_raster(hjust = 0, vjust = 0) +
  theme(text = element_text(size = 20)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = -0.5, hjust = 1.25),
    axis.text.y = element_text(vjust = 2)
  ) +
  labs(x = "mtDNA position (nt)", y = "", fill = "Allele \nFrequency (%)") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))
saveManuscriptPlot(Fig13_allele_frequency_grid, 20,5)
Fig13_allele_frequency_grid
```

## snpEff predicted impact
```{r setup12, warning=FALSE}
# snpEff predicted impact heatmap

Fig14_variant_impact_grid <- variants2 %>%
  select(sample_name, POS, base_change, AF, variant_type) %>%
  mutate(variant_type = case_when(
    variant_type == "frameshift_variant" ~ "frameshift",
    variant_type == "non_protein_coding" ~ "non protein coding",
    variant_type == "synonymous_variant" ~ "synonymous",
    variant_type == "missense_variant" ~ "missense",
    variant_type == "stop_gained" ~ "stop gained"
  )) %>%
  ggplot(aes(x = as.factor(POS), y = as.factor(sample_name), fill = variant_type)) +
  geom_raster(hjust = 0, vjust = 0) +
  theme(text = element_text(size = 20)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = -0.5, hjust = 1.25),
    axis.text.y = element_text(vjust = 2),
    legend.position = "bottom"
  ) +
  labs(x = "mtDNA position (nt)", y = "", fill = "") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))
saveManuscriptPlot(Fig14_variant_impact_grid, 20,5)
Fig14_variant_impact_grid

```

# Violin plots of "most variable" mutations
###### Mutation must be heteroplasmic in at least 2 single cells

```{r setup13, warning=FALSE}
# violin plots "most variable"


violinvariants <- variants2 %>%
filter(!sample_name == "bulk")

violinvariants$AF <- as.numeric(as.character(violinvariants$AF))

# only keeping mutations with at least 2 cells

by_TEST <- violinvariants %>%
   group_by(mutation) %>%
   filter(n() > 1) %>%
   ungroup

mutvariancelistabove2 <- by_TEST$mutation

violinvariants2 <- filter(violinvariants, mutation %in% mutvariancelistabove2)

# reorder variants

violinvariants2$mutation <- factor(violinvariants2$mutation, levels=c("5244 TA>T", "5462 T>C", "6262 TA>T", "7001 A>G", "7853 A>C", "14136 GA>G"))

Fig_ONLY3 <- ggplot(violinvariants2, aes(x=mutation, y=AF)) + 
    geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_hline(yintercept=0.04, linetype="dashed", color = "#CC79A7", size = 1.5) +
  geom_hline(yintercept=0.96, linetype="dashed", color = "#CC79A7", size = 1.5) +
  theme_minimal() +
  scale_fill_grey() +
  ylab("Allele Frequency") +
  scale_colour_grey() +
  xlab("Most Variable Mutations") +
  theme(text = element_text(size = 16))
saveManuscriptPlot(Fig_ONLY3, 20,5)
Fig_ONLY3
```

# Barchart of average single cell allele frequency compared to bulk allele frequency
###### Mutation must be heteroplasmic in bulk sample

```{r setup1433, warning=FALSE}
# barchart of averaged single cell AF compared to bulk AF

variants2$AF <- as.numeric(as.character(variants2$AF))

# filter only for bulk samples

variantsJUSTBULK <- variants2 %>%
filter(sample_name == "Bulk") %>%
select(mutation, AF)

variantsJUSTBULK$Type <- "Bulk"

# create list for bulk mutations

BULKvariantlist <- unique(variantsJUSTBULK$mutation)

variantsJUSTSC <- variants2 %>%
ungroup %>%
add_row(mutation = "3733 G>A", AF = 0, sample_name = "Single Cell 1")

# generate average for single cell AF

variantsJUSTSC <- variantsJUSTSC %>%
filter(!sample_name == "Bulk") %>%
filter(mutation %in% BULKvariantlist) %>%
group_by(mutation) %>%
summarise(AF = mean(AF, na.rm = TRUE))

variantsJUSTSC$Type <- "Single Cell Average"


# combine tables

variantsJOINT <- rbind(variantsJUSTBULK, variantsJUSTSC) %>%
ungroup %>%
select(mutation, AF, Type)

variantsJOINT$mutation <- factor(variantsJOINT$mutation, levels=c("3733 G>A", "5244 TA>T", "5462 T>C", "6262 TA>T", "7001 A>G", "7853 A>C", "14136 GA>G"))


# ggplot

Fig_comparebulktosc <- ggplot(variantsJOINT, aes(fill=Type, y=AF, x=mutation)) + 
    geom_bar(position="dodge", stat="identity") +
    ylab("Allele Frequency") +
    xlab("Mutation") +
    theme_minimal() +
    scale_fill_grey() +
    ylab("Allele Frequency") +
    xlab("Mutations Found in Bulk") +
    scale_colour_grey() +
    theme(text = element_text(size = 13), legend.position="bottom") +
  geom_hline(yintercept=0.04, linetype="dashed", color = "#CC79A7", size = 1.5) +
  geom_hline(yintercept=0.96, linetype="dashed", color = "#CC79A7", size = 1.5) +
  ylim(0,1)
Fig_comparebulktosc
saveManuscriptPlot(Fig_comparebulktosc, 20,5)
```

# Session info
```{r setup15, warning=FALSE}
sessionInfo()
```