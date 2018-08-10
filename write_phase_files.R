#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

parser <- ArgumentParser()

parser$add_argument("-g", "--hlagenos", type = "character", 
		    help = "HLA alleles for each individual")

parser$add_argument("-v", "--vcf", type = "character", default = NULL, 
		    help = "VCF file with known phase variants (Optional)")

parser$add_argument("-c", "--genecoord", type = "character", 
		    help = "HLA genes start position")

parser$add_argument("-o", "--outdir", type = "character", 
		    help = "Output directory")

args <- parser$parse_args()
hla <- args$hlagenos
vcf <- args$vcf 
hcoord <- args$genecoord
outdir <- args$outdir

if (!dir.exists(outdir)) dir.create(outdir)

hla_genos <- read_tsv(hla) 

hla_genos_recoded <- hla_genos %>%
    select(locus, allele) %>%
    arrange(locus, allele) %>%
    distinct() %>%
    group_by(locus) %>%
    mutate(code = seq_len(n())) %>%
    left_join(hla_genos, ., by = c("locus", "allele"))

write_tsv(hla_genos_recoded, file.path(outdir, "hla_allele_codes.txt"))

hla_coords <- read_tsv(hcoord)
    
genos_df <- hla_genos_recoded %>%
    group_by(subject, locus) %>%
    summarise(geno = paste(code, collapse = "|")) %>%
    ungroup() %>%
    left_join(hla_coords, by = "locus") %>%
    mutate(marker_type = "M") %>%
    select(subject, locus, pos, marker_type, geno)

if (!is.null(vcf)) {
    snp_genos <- read_tsv(vcf, comment = "##") %>%
	select(grep("^id$", names(.), ignore.case = TRUE),
	       grep("^pos$", names(.), ignore.case = TRUE),
	       sort(unique(hla_genos$subject))) %>%
	gather(subject, geno, -(1:2)) %>%
	mutate(marker_type = "S") %>%
	select(subject, locus = 1, pos = 2, marker_type, geno)

    genos_df <- genos_df %>% bind_rows(snp_genos)
} 

genos_df <- genos_df %>%
    arrange(subject, pos) %>%
    mutate(ix = "1|2") %>%
    separate_rows(ix, geno, sep = "\\|") %>%
    rename(allele = geno)

var_df <- genos_df %>%
    distinct(pos, locus, marker_type) %>%
    arrange(pos) %>%
    mutate(known = ifelse(marker_type == "S", 0, "*"))

write_tsv(var_df, file.path(outdir, "input_loci.tsv"))

n_ind <- n_distinct(genos_df$subject)
n_loci <- nrow(var_df)
P <- c("P", var_df$pos) %>% paste(collapse = " ")
m_types <- var_df$marker_type %>% paste(collapse = "")

haps <- genos_df %>%
    select(subject, pos, ix, allele) %>%
    arrange(subject, pos, ix) %>%
    spread(pos, allele) %>%
    select(-ix)

phaseinp <- file.path(outdir, "phase.inp")
phaseknown <- file.path(outdir, "phase.known")

file.create(phaseinp)
file.create(phaseknown)
write(n_ind, phaseinp, append = TRUE)
write(n_loci, phaseinp, append = TRUE)
write(P, phaseinp, append = TRUE)
write(m_types, phaseinp, append = TRUE)

haps %>%
    split(.$subject) %>%
    walk(~.x %>% (function(x) {
	write(paste0("#", unique(x$subject)), phaseinp, append = TRUE)
	write(unlist(x[1, -1]), phaseinp, ncolumns = n_loci, append = TRUE)
	write(unlist(x[2, -1]), phaseinp, ncolumns = n_loci, append = TRUE)
	write(paste(var_df$known, collapse = ""), phaseknown, append = TRUE)}))
