# Scripts used during analysis of tetrasomic regions

## Requirements
* minimap2
* syri
* R
* R packages:
  * pafr
  * data.table
* Included R files
  * mm2_make_filt_paf_for_syri.r
    * used for filtering alignment results
  * window_functions_L8.r
    * contains functions used for sliding window analysis of alignments

## Input files
* Line8 assembly fasta
  * use a version that has Arahy.09_alt removed
  * `Line8.fa`
* Line8 subgenome A fasta
  * contains Arahy.01-Arahy.02
  * `Line8_subA.fa`
* Line8 subgenome B fasta
  * contains Arahy.11-Arahy.20
  * `Line8_subB.fa`
* "Peanut ancestor" fasta
  * A.duranensis and A.ipaensis assemblies concatenated into single fasta
  * `A_B_ancestors.fa`

## Analysis Steps
1. Whole genome alignments with minimap2
2. Filter alignments with custom R script
3. Use SyRi to characterize differences between genomes/chromosomes
4. Calculate sliding-window values of alignment quality

## Whole genome alignments
```
minimap2 -cx asm10 --eqx -t 20 Line8.fa A_B_ancestors.fa > ancestors_v_L8.paf

minimap2 -cx asm10 --eqx -t 20 Line8_subA.fa A_B_ancestors.fa > \
ancestors_v_L8subA.paf

minimap2 -cx asm10 --eqx -t 20 Line8_subB.fa A_B_ancestors.fa > \
ancestors_v_L8subB.paf

# combine alignments against each subgenome into single file
cat ancestors_v_L8subA.paf ancestors_v_L8subB.paf > \
ancestors_v_L8subgenomes.paf

minimap2 -cx asm20 --eqx -t 20 Line8_subA.fa Line8_subB.fa > L8_subB_v_subA.paf
```

## Filter alignments
### Filter Line8 subgenome B vs subgenome A alignment
* uses `mm2_make_filt_paf_for_syri.r`
* requires `pafr` R package
* uses following filtering criteria
  * 80% sequence identity of higher
  * longer than 5kb
  * only primary alignments
* generates file called: L8_subB_v_subA.filt_0.8_5000_primaryOnly.paf
```
Rscript mm2_make_filt_paf_for_syri.r L8_subB_v_subA.paf 0.8 5000 T
```
### Make alignment files that are more easily processed by R
```
cut -f 1-12 ancestors_v_L8.paf > ancestors_v_L8.paf_info.txt

cut -f 1-12 ancestors_v_L8subgenomes.paf > \
ancestors_v_L8subgenomes.paf_info.txt

cut -f 1-12 L8_subB_v_subA.filt_0.8_5000_primaryOnly.paf > \
L8_subB_v_subA.filt_0.8_5000_primaryOnly.paf_info.txt
```

## Run Syri

## Calculate sliding window alignment quality values using R
* requires functions in `window_functions_L8.r`
* requires `data.table` R package
```
### LOAD PACKAGES ###
library(data.table)
source('window_functions_L8.r')

### SET INPUTS ###
# competitive mapping of ancestors to Line8
comp_paf <- fread('ancestors_v_L8.paf_info.txt')

# results from mapping ancestors to Line8 subgenomes separately
subg_paf <- fread('ancestors_v_L8subgenomes.paf_info.txt')

### SET VARIABLES ###
# sliding window size
window_size <- 100000

# alignment cutoff settings
perID_cut <- 0.7
len_cut <- 10000
qual_cut <- 40

######

# generate sliding window alignment values (using filtered alignments)

# competitive alignment results
paf_in <- comp_paf
paf_window_list <- list()

for(chrn in sort(unique(paf_in$V6))){
  #print(chrn)
  both_sub1 <- paf_in[V6 == chrn, ]
  for(sgv in c('aradu', 'araip')){
    #print(sgv)
    sub1 <- both_sub1[grep(sgv, V1), ]
    #
    tmp_quer_chr <- names(which.max(table(sub1$V1)))
    tmp_ref_chr <- names(which.max(table(sub1$V6)))
    tmp_quer_length <- as.integer(names(which.max(table(sub1$V2))))
    tmp_ref_length <- as.integer(names(which.max(table(sub1$V7))))
    #
    tmp_start_vals <- seq(from = 1, to = max(sub1$V7), by = window_size)
    sub2 <- lapply(tmp_start_vals, function(x) get_paf_window(sub1, x,
      x+window_size-1))
    sub3 <- fill_empty_windows(sub2, tmp_start_vals, window_size)
    sub4 <- lapply(sub3, function(x)
      if(nrow(x) > 1){
        filter_paf_window(x)
      } else{x})
    sub5 <- lapply(sub4, function(x)
      trim_window(x, x[1,WINDOW_START], x[1,WINDOW_END]))
    sub6 <- lapply(sub5, function(x) window_quality_filt(x, perID_cut, len_cut,
      qual_cut))
    # this step is for if every alignment is filtered out for a chromosome
    if(sum(unlist(lapply(sub6, nrow))) == 0){
      sub6 <- list()
      sub6[[1]] <- make_empty_paf(tmp_quer_chr, tmp_quer_length, tmp_ref_chr,
                     tmp_ref_length)
    }
    sub7 <- comb_window_pafs(sub6, tmp_start_vals, window_size)
    paf_window_list[[sgv]][[chrn]] <- sub7
  }
}

# alignments of ancestors against each subgenome separately
paf_in <- subg_paf
sub_window_list <- list()

for(chrn in sort(unique(paf_in$V6))){
  #print(chrn)
  both_sub1 <- paf_in[V6 == chrn, ]
  for(sgv in c('aradu', 'araip')){
    #print(sgv)
    sub1 <- both_sub1[grep(sgv, V1), ]
    #
    tmp_quer_chr <- names(which.max(table(sub1$V1)))
    tmp_ref_chr <- names(which.max(table(sub1$V6)))
    tmp_quer_length <- as.integer(names(which.max(table(sub1$V2))))
    tmp_ref_length <- as.integer(names(which.max(table(sub1$V7))))
    #
    tmp_start_vals <- seq(from = 1, to = max(sub1$V7), by = window_size)
    sub2 <- lapply(tmp_start_vals, function(x) get_paf_window(sub1, x,
      x+window_size-1))
    sub3 <- fill_empty_windows(sub2, tmp_start_vals, window_size)
    sub4 <- lapply(sub3, function(x)
      if(nrow(x) > 1){
        filter_paf_window(x)
      } else{x})
    sub5 <- lapply(sub4, function(x)
      trim_window(x, x[1,WINDOW_START], x[1,WINDOW_END]))
    sub6 <- lapply(sub5, function(x) window_quality_filt(x, perID_cut, len_cut,
      qual_cut))
    # this step is for if every alignment is filtered out for a chromosome
    if(sum(unlist(lapply(sub6, nrow))) == 0){
      sub6 <- list()
      sub6[[1]] <- make_empty_paf(tmp_quer_chr, tmp_quer_length, tmp_ref_chr,
                     tmp_ref_length)
    }
    sub7 <- comb_window_pafs(sub6, tmp_start_vals, window_size)
    sub_window_list[[sgv]][[chrn]] <- sub7
  }
}

# I recommend saving these lists to avoid having to make them again

```
