# Script to generate filtered PAF to use for SyRi 
#   Filter PAF based on percent identity, alignment length, and primary
#     alignment
#############
# REQUIREMENTS
# pafr package
#
# USAGE
# Rscript mm2_make_filt_paf_for_syri.r PAF_FILE PER_IDY_CUT LEN_CUT PRIM_ONLY 
#
# INPUTS
# args[1] = PAF_FILE = paf_file
#   The paf alignment file from minimap2
# args[2] = PER_IDY_CUT = per_idy_cut
#   The percent identity cutoff; ex: 0.85
# args[3] = LEN_CUT = len_cut
#   The alignment length cutoff; ex: 1000 or 1e4
# args[4] = PRIM_ONLY = prim_only
#   Whether to exclude secondary alignments; ex: T or F
#
# OUTPUT
# filtered PAF file, including cigar strings, that can be used with Syri  
##########

args <- commandArgs(trailingOnly = T)

### LOAD PACKAGES ###
library(pafr)

### SET INPUTS ###
paf_file <- args[1]

### SET VARIABLES ###
# percent identity cutoff
per_idy_cut <- as.numeric(args[2])

# length cutoff
len_cut <- as.numeric(args[3])

# use only primary alignments
prim_only <- as.logical(args[4])

### SET OUTPUTS ###
per_idy_char <- as.character(per_idy_cut)
len_cut_char <- as.character(len_cut)

if(prim_only){
  prim_only_char <- 'primaryOnly'
} else{
  prim_only_char <- 'withSecondary'
}

filt_suf <- paste0('.filt_',
  paste(per_idy_char, len_cut_char, prim_only_char, sep = '_'),
  '.paf')

out_paf <- sub('.paf', filt_suf, paf_file, fixed = T)

##########

ali_0 <- read_paf(paf_file)
ali_0$PER_IDY <- ali_0$nmatch / ali_0$alen

if(prim_only){
  tp_keep <- which(ali_0$tp != 'S' & ali_0$tp != 'i')
} else{
  tp_keep <- seq(nrow(ali_0))
}

per_idy_inds <- which(ali_0$PER_IDY >= per_idy_cut)
len_inds <- which(ali_0$alen >= len_cut)

filt_inds <- sort(intersect(intersect(per_idy_inds, len_inds), tp_keep)) 

ali_keep <- ali_0[filt_inds, ]

ali_df <- as.data.frame(ali_keep[ ,c(1:12)])

tmp_cg <- paste0('cg:Z:', ali_keep$cg)

ali_df$cg <- tmp_cg

write.table(ali_df, file = out_paf, quote = F, sep = '\t',
  row.names = F, col.names = F)

quit(save = 'no')
