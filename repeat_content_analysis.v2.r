# Analysis in R of repetitive sequence content in Line 8 assembly/annotation

### LOAD PACKAGES ###
library(data.table)
library(GENESPACE)
library(Biostrings)

### SET INPUTS ###
# fasta file of Line 8 assembly
faFile <- 'Line8.fasta'
dnaSS <- Biostrings::readDNAStringSet(faFile)
names(dnaSS) <- sapply(names(dnaSS), function(x) strsplit(x, " ")[[1]][1])
dnaSS <- dnaSS[Biostrings::width(dnaSS) > 1e6]
seqInfo <- pull_seqInfo(dnaSS)

# index file of assembly fasta; can be made using `samtools faidx Line8.fasta`
fai_file <- 'Line8.fasta.fai'
l8_fai <- fread(fai_file)

# Line 8 gene annotation gff3 file
# for this script, need to add the following line to the gff3 header, just before first line of data:
# `seqid	source	type	start	end	score	strand	phase	attributes`
geneGffFile <- 'Line8_gene.gff3'
genes <- as.data.frame(rtracklayer::readGFF(
  geneGffFile,
  columns = c("seqid", "start", "end", "type"),
  tags = "Parent",
  filter = list(type = c("mRNA", "CDS"))))
genes <- subset(genes, seqid %in% names(dnaSS))

# Line 8 repeat annotation gff3 file
repGFFFile <- 'Line8_repeats.gff3'
repGFF <- fread(repGffFile, skip = "Arahy.01")
repeats <- repGFF[, c(1,4,5,9)]
colnames(repeats) <- c("chr", "start", "end", "type")
repeats <- subset(repeats, chr %in% names(dnaSS))

###############

bedList <- list(
  Exons = with(subset(genes, type == "CDS"), data.frame(
    chr = seqid, start = start, end = end)),
  LINE.SINE = with(subset(repeats, grepl("LINE|SINE", type)), data.frame(
    chr = chr, start = start, end = end)),
  Ty3LTR = with(subset(repeats, grepl("Gypsy", type)), data.frame(
    chr = chr, start = start, end = end)),
  OtherLTR = with(subset(repeats,
    (grepl("LTR", type) & !grepl("Gypsy", type))), data.frame(
    chr = chr, start = start, end = end)),
  Satellite = with(subset(repeats, grepl("Simple_repeat", type)), data.frame(
    chr = chr, start = start, end = end)),
  OtherRepeats = with(subset(repeats,
      !grepl("LINE|SINE|Simple_repeat|LTR", type)), data.frame(
    chr = chr, start = start, end = end)),
  Introns = with(subset(genes, type == "mRNA"), data.frame(
    chr = seqid, start = start, end = end)))

suppressWarnings(genomeClasses <- classify_genome(
  dnaSS = dnaSS, listOfBeds = bedList, verbose = T))

seq_amount <- lapply(genomeClasses, function(x)
  sum(width(x))/1e6)

seq_dt <- data.table(
  ELEMENT = names(seq_amount),
  TOT_LENGTH = unlist(seq_amount),
  PER_CONTENT = unlist(seq_amount)/sum(unlist(seq_amount)))
seq_dt[ELEMENT == 'unknown', ELEMENT := 'Unannotated']

# seq_dt can be used to analysis genome-wide patterns

#####
# Split by chromosome
gC_chr <- lapply(genomeClasses, function(x)
  split(x, seqnames(x)))
gC_chr_width <- lapply(gC_chr, function(x)
  lapply(x, function(y) sum(width(y))/1e6))

chr_cont_df <- data.frame(data =
  matrix(unlist(gC_chr_width), byrow = F, ncol = length(gC_chr_width)))
colnames(chr_cont_df) <- names(gC_chr_width)
chr_cont_dt <- data.table(chr_cont_df)
chr_cont_dt[, CHROM := names(gC_chr_width[[1]])]

# chr_cont_dt can be used to compare patterns between chromosomes

