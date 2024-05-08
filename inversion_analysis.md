# Scripts to during analysis of inversions between peanut assemblies

## Requirements
* minimap2
* syri
* seqtk
* R
* R packages:
  * pafr
  * data.table
* Included R files
  * `mm2_make_filt_paf_for_syri.r`
    * used for filtering alignment results
  * `syri_analysis_functions_L8.r`
    * functions used for analysis  

## Input files
* fasta files for peanut assemblies
  * file names will depend on where you get the files from
  * Names used below:
    * `Line8.fa`
    * `Tifrunner.fa`
    * `Shitouqi.fa`
    * `BaileyII.fa`
  * Preparation 
    * Works best if make adjusted fasta files for each assembly 
      * Line8: remove Arahy.09_alt
      * Tifrunner, Shitouqi, and Bailey II: remove non-chromosome scaffolds
* GFF3 gene file for Line8 from annotation
  * file name may be different
  * `AhypogaeaLine8_886_v1.3.gene.gff3`

## Steps
1. Whole genome alignments with minimap2
2. Filter alignments with custom R script
3. Use SyRI to characterize differences between assemblies
4. Find locations of gaps between contigs
5. Use R to
  * Calculate distances between contig gaps and inversion boundaries
  * Identify number of genes in inversions

## Alignments using minimap2
```
minimap2 -cx asm5 --eqx -t 20 Line8.fa Tifrunner.fa > Tif_v_L8.paf
minimap2 -cx asm5 --eqx -t 20 Line8.fa Shitouqi.fa > Shi_v_L8.paf
minimap2 -cx asm5 --eqx -t 20 Line8.fa BaileyII.fa > B2_v_L8.paf
```

## Filter alignments using R script
* script requires "pafr" R package
* Filtering criteria
  * 0.85 percent identity
  * 2kb or longer alignments
  * only primary alignments
* Output files have same prefix as .paf file but with
  ".filt_0.85_2000_primaryOnly.paf" suffix
```
Rscript mm2_make_filt_paf_for_syri.r Tif_v_L8.paf 0.85 2000 T
Rscript mm2_make_filt_paf_for_syri.r Shi_v_L8.paf 0.85 2000 T
Rscript mm2_make_filt_paf_for_syri.r B2_v_L8.paf 0.85 2000 T
```

## Use SyRI to characterize differences between genomes
```
syri -c Tif_v_L8.filt_0.85_2000_primaryOnly.paf -r Line8.fa \
-q Tifrunner.fa -k -F P --prefix Tif_v_L8.
syri -c Shi_v_L8.filt_0.85_2000_primaryOnly.paf -r Line8.fa \
-q Shitouqi.fa -k -F P --prefix Shi_v_L8.
syri -c B2_v_L8.filt_0.85_2000_primaryOnly.paf -r Line8.fa \
-q BaileyII.fa -k -F P --prefix B2_v_L8.
```

## Generate .bed files for contig gaps in analyzed assebmlies
* Bailey II has some gaps that are shorter than 50bp so use shorter gap length
```
seqtk gap -l 50 Line8.fa > Line8_gaps.bed
seqtk gap -l 50 Tifrunner.fa > Tif_gaps.bed
seqtk gap -l 50 Shitouqi.fa > Shi_gaps.bed
seqtk gap -l 10 BaileyII.fa > B2_gaps.bed
```

## Analysis of contig breaks vs inversion boundaries
* in R
  * requires `data.table` package
  * requires functions in `syri_analysis_functions_L8.r`
```
### LOAD PACKAGES ###
library(data.table)
source('syri_analysis_functions_L8.r')

### SET INPUTS ###
# load contig bread .bed files
n_bedfile_list <- list()
n_bedfile_list[['line8']] <- 'Line8_gaps.bed'
n_bedfile_list[['tifrunner']] <- 'Tif_gaps.bed'
n_bedfile_list[['shitouqi']] <- 'Shi_gaps.bed'
n_bedfile_list[['bailey']] <- 'B2_gaps.bed'

n_bed_list <- lapply(n_bedfile_list, fread)

# load SyRI results
syri_resfile_list <- list()
syri_resfile_list[['tifrunner_v_line8']] <- 'Tif_v_L8.syri.out'
syri_resfile_list[['shitouqi_v_line8']] <- 'Shi_v_L8.syri.out'
syri_resfile_list[['bailey_v_line8']] <- 'B2_v_L8.syri.out'

syri_res_list <- lapply(syri_resfile_list, fread)

# load Line8 gene gff3
line8_gff_in <- 'AhypogaeaLine8_886_v1.3.gene.gff3'
l8_gff <- fread(line8_gff_in, skip = 'Arahy.01')

#########

###
# process syri results
###
syri_proc_list <- lapply(syri_res_list, process_syri_res_keep_all)

###
# adjust the contig gap bed file chromosome names
###
n_bed_list[['line8']][, SYRI_CHR := V1]
n_bed_list[['tifrunner']][, SYRI_CHR := V1]
n_bed_list[['shitouqi']][,
  SYRI_CHR := gsub('arahy.Shitouqi.gnm1.chr', 'Chr', V1)]

bailey_chrs <- sort(unique(n_bed_list[['bailey']]$V1))
names(bailey_chrs) <- paste0('Chr',
  formatC(seq(length(bailey_chrs)), width = 2, format = 'd', flag = '0'))
n_bed_list[['bailey']][, SYRI_CHR := V1]
for(i in seq(20)){
  n_bed_list[['bailey']][SYRI_CHR == bailey_chrs[i],
    SYRI_CHR := names(bailey_chrs)[i]]
}

###
# find distances between contig breaks and inversions
###
tif_l8_inv_contig_info <- get_el_gap_dists_sameChr(
  syri_dt = syri_res_list[['tifrunner_v_line8']],
  ref_bed = n_bed_list[['line8']],
  quer_bed = n_bed_list[['tifrunner']]
)

shi_l8_inv_contig_info <- get_el_gap_dists_sameChr(
  syri_dt = syri_res_list[['shitouqi_v_line8']],
  ref_bed = n_bed_list[['line8']],
  quer_bed = n_bed_list[['shitouqi']]
)

b2_l8_inv_contig_info  <- get_el_gap_dists_sameChr(
  syri_dt = syri_res_list[['bailey_v_line8']],
  ref_bed = n_bed_list[['line8']],
  quer_bed = n_bed_list[['bailey']]
)

###
# data tables can be used to look at distances between contig breakpoints
#   and inversion boundaries

###
# Filter inversions to include only those with boundaries more than 
#    2kb from contig breakpoints
### 

tif_line8_good_invs <- tif_l8_inv_contig_info[REF_MIN_GAP_DIST > 2000 &
  QUER_MIN_GAP_DIST > 2000, ]

shi_line8_good_invs <- shi_l8_inv_contig_info[REF_MIN_GAP_DIST > gap_dist_cut &
  QUER_MIN_GAP_DIST > gap_dist_cut, ]

bailey_line8_good_invs <- bailey_line8_inv_info[
  REF_MIN_GAP_DIST > gap_dist_cut & QUER_MIN_GAP_DIST > gap_dist_cut, ]

###
# find Line8 genes in filtered inversions
###

tif_line8_good_invs[, N_L8_GENES := as.numeric(NA)]
shi_line8_good_invs[, N_L8_GENES := as.numeric(NA)]
bailey_line8_good_invs[, N_L8_GENES := as.numeric(NA)]

test_gff <- l8_gff
test_inv_dt <- tif_line8_good_invs

for(i in seq(nrow(tif_line8_good_invs))){
  tmp_gene_tab <- test_gff[V3 == 'gene' & V1 == test_inv_dt$V1[i] &
      V4 >= test_inv_dt$V2[i] & V5 <= test_inv_dt$V3[i], ]
  tmp_gene_tab_2 <- test_gff[V3 == 'gene' & V1 == test_inv_dt$V1[i] &
      V4 < test_inv_dt$V2[i] & V5 > test_inv_dt$V3[i], ]
  tmp_gene_tab_3 <- test_gff[V3 == 'gene' & V1 == test_inv_dt$V1[i] &
      V4 < test_inv_dt$V2[i] & V5 > test_inv_dt$V2[i] &
      V5 < test_inv_dt$V3[i], ]
  tmp_gene_tab_4 <- test_gff[V3 == 'gene' & V1 == test_inv_dt$V1[i] &
      V4 > test_inv_dt$V2[i] & V4 < test_inv_dt$V3[i] &
      V5 > test_inv_dt$V3[i], ]
  tmp_n_genes <- nrow(tmp_gene_tab) + nrow(tmp_gene_tab_2) +
    nrow(tmp_gene_tab_3) + nrow(tmp_gene_tab_4)
  tif_line8_good_invs[i, N_L8_GENES :=
    tmp_n_genes]
}

test_inv_dt <- shi_line8_good_invs

for(i in seq(nrow(shi_line8_good_invs))){
  tmp_gene_tab <- test_gff[V3 == 'gene' & V1 == test_inv_dt$V1[i] &
      V4 >= test_inv_dt$V2[i] & V5 <= test_inv_dt$V3[i], ]
  tmp_gene_tab_2 <- test_gff[V3 == 'gene' & V1 == test_inv_dt$V1[i] &
      V4 < test_inv_dt$V2[i] & V5 > test_inv_dt$V3[i], ]
  tmp_gene_tab_3 <- test_gff[V3 == 'gene' & V1 == test_inv_dt$V1[i] &
      V4 < test_inv_dt$V2[i] & V5 > test_inv_dt$V2[i] &
      V5 < test_inv_dt$V3[i], ]
  tmp_gene_tab_4 <- test_gff[V3 == 'gene' & V1 == test_inv_dt$V1[i] &
      V4 > test_inv_dt$V2[i] & V4 < test_inv_dt$V3[i] &
      V5 > test_inv_dt$V3[i], ]
  tmp_n_genes <- nrow(tmp_gene_tab) + nrow(tmp_gene_tab_2) +
    nrow(tmp_gene_tab_3) + nrow(tmp_gene_tab_4)
  shi_line8_good_invs[i, N_L8_GENES :=
    tmp_n_genes]
}

test_inv_dt <- bailey_line8_good_invs

for(i in seq(nrow(bailey_line8_good_invs))){
  tmp_gene_tab <- test_gff[V3 == 'gene' & V1 == test_inv_dt$V1[i] &
      V4 >= test_inv_dt$V2[i] & V5 <= test_inv_dt$V3[i], ]
  tmp_gene_tab_2 <- test_gff[V3 == 'gene' & V1 == test_inv_dt$V1[i] &
      V4 < test_inv_dt$V2[i] & V5 > test_inv_dt$V3[i], ]
  tmp_gene_tab_3 <- test_gff[V3 == 'gene' & V1 == test_inv_dt$V1[i] &
      V4 < test_inv_dt$V2[i] & V5 > test_inv_dt$V2[i] &
      V5 < test_inv_dt$V3[i], ]
  tmp_gene_tab_4 <- test_gff[V3 == 'gene' & V1 == test_inv_dt$V1[i] &
      V4 > test_inv_dt$V2[i] & V4 < test_inv_dt$V3[i] &
      V5 > test_inv_dt$V3[i], ]
  tmp_n_genes <- nrow(tmp_gene_tab) + nrow(tmp_gene_tab_2) +
    nrow(tmp_gene_tab_3) + nrow(tmp_gene_tab_4)
  bailey_line8_good_invs[i, N_L8_GENES :=
    tmp_n_genes]
}

###
# objects can be used examine the number and density of genes within
#   inversions

```
