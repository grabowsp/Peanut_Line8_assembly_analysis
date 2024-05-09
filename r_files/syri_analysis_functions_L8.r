# Function to use when analyzing syri results

process_syri_res_keep_all <- function(raw_syri_res){
  #############
  # Process the syri results to change the position columns to numeric,
  #   calculate the length of the elements in both assembies
  #   Retains ALL elements, so can't be used for some downstream functions
  # INPUTS
  # raw_syri_res = data.table; the raw syri results
  # OUTPUTS
  # data.table of raw_syri_res with position columns adjusted to numeric
  #   and columns for the length of elements in the REF and QUER genomes
  #############
  tmp_syri_res <- raw_syri_res
  tmp_syri_res[, V2 := as.integer(V2)]
  tmp_syri_res[, V3 := as.integer(V3)]
  tmp_syri_res[, V7 := as.integer(V7)]
  tmp_syri_res[, V8 := as.integer(V8)]
  #
  tmp_syri_res[, REF_LEN := V3-V2]
  tmp_syri_res[, QUER_LEN := V8-V7]
  return(tmp_syri_res)
}

#####

calc_el_gap_dist <- function(el_start_pos, el_end_pos, gap_start_pos,
  gap_end_pos){
  ###########
  # Calculate the minimum distance from beginning and end of genomic elements 
  #       to gaps in an alignment; originaly designed for inversionts
  #  All position vectors are expected to come from the same chromosome
  #
  # INPUTS
  # el_start_pos        = integer vector, the start positions of elements
  # el_end_pos  = integer vector, the end positions of elements; should line
  #                             up with el_start_pos
  # gap_start_pos       = integer vector, the start positions of gaps in the
  #                             assembly
  # gap_end_pos         = integer vector, the end positions of gaps in the
  #                             assembly
  # OUTPUT
  # el_gap_dist_list = list with 2 elements
  # el_gap_dist_list[[1]] = vector of smallest distance of element start 
  #                             points to a gap in the assembly
  # el_gap_dist_list[[2]] = vector of smallest distance of element end 
  #                             points to a gap in the assembly
  ##########
  # inversion start positions
  el_start_gap_start <- sapply(el_start_pos,
    function(x) gap_start_pos - x)
  el_start_gap_end <- sapply(el_start_pos,
    function(x) x - gap_end_pos)
  el_start_combo <- rbind(el_start_gap_start, el_start_gap_end)
  el_start_min_dist <- apply(el_start_combo, 2, function(x) min(abs(x)))
  # check if beginning is within a gap
  el_start_gap_vec <- c()
  for(j in seq(ncol(el_start_gap_start))){
    tmp_in <- length(intersect(which(el_start_gap_start[, j] > 0),
                                which(el_start_gap_end[, j] > 0)))
    el_start_gap_vec <- c(el_start_gap_vec, tmp_in)
  }
  if(sum(el_start_gap_vec) > 0){
    el_start_in_gap <- which(el_start_gap_vec > 0)
    el_start_min_dist[el_start_in_gap] <- 0
  }
  # inversion end positions
  el_end_gap_start <- sapply(el_end_pos,
    function(x) gap_start_pos - x)
  el_end_gap_end <- sapply(el_end_pos,
    function(x) x - gap_end_pos)
  el_end_combo <- rbind(el_end_gap_start, el_end_gap_end)
  el_end_min_dist <- apply(el_end_combo, 2, function(x) min(abs(x)))
  # check if end is within a gap
  el_end_gap_vec <- c()
  for(j in seq(ncol(el_end_gap_start))){
    tmp_in <- length(intersect(which(el_end_gap_start[, j] > 0),
                                which(el_end_gap_end[, j] > 0)))
    el_end_gap_vec <- c(el_end_gap_vec, tmp_in)
  }
  if(sum(el_end_gap_vec) > 0){
    el_end_in_gap <- which(el_end_gap_vec > 0)
    el_end_min_dist[el_end_in_gap] <- 0
  }
  el_gap_dist_list <- list()
  el_gap_dist_list[[1]] <- el_start_min_dist
  el_gap_dist_list[[2]] <- el_end_min_dist
  return(el_gap_dist_list)
}

######

get_el_gap_dists_sameChr <- function(syri_dt, ref_bed, quer_bed){
  ##########
  # Get the distance of element start and end positions to gaps in assemblies
  #     for elements on the same chromosome in both assemblies
  #   Designed for looking at the same type of element at a time, but might
  #     work if looking at multiple element types at a time
  #   A small set of elements sometimes have end before start positions; not
  #     sure how this function will handle those cases - be cautious; thses
  #     show up as elements with a negative length
  #
  # INPUTS
  # syri_dt     = data.table; results from syri and processed with the 
  #                     process_syri_res() function; was designed for syri_dt
  #                     to be a subset of larger syri file and only contain
  #                     one type of element, like INV or TRANS
  # ref_bed     = data.table; bed file of position of gaps in reference
  #                     assembly used for syri; MUST have 'SYRI_CHR' column
  #                     where chromosome name is exact same format as in the
  #                     syri file
  # quer        = data.table; bed file of position of gaps in query
  #                     assembly used for syri; MUST have 'SYRI_CHR' column
  #                     where chromosome name is exact same format as in the
  #                     syri file
  #
  # OUTPUT
  # el_dist_out_dt = data.table; syri_dt with 4 extra columns showing
  #     distance of start and end of element to gaps in both assemblies;
  #  Distance of 0 means position is in a gap
  ###########
  el_dist_out_list <- list()
  #
  for(cn in unique(syri_dt$V1)){
    tmp_syri_chr <- syri_dt[V1 == cn,]
    #
    ref_gap_start_pos <- ref_bed[SYRI_CHR == cn, V2]
    ref_gap_end_pos <- ref_bed[SYRI_CHR == cn, V3]
    #
    ref_el_start_pos <- tmp_syri_chr$V2
    ref_el_end_pos <- tmp_syri_chr$V3
    #
    ref_el_gap_dists <- calc_el_gap_dist(el_start_pos = ref_el_start_pos,
                           el_end_pos = ref_el_end_pos,
                           gap_start_pos = ref_gap_start_pos,
                           gap_end_pos = ref_gap_end_pos)
    #
    quer_cn <- tmp_syri_chr[1, V6]
    quer_gap_start_pos <- quer_bed[SYRI_CHR == quer_cn, V2]
    quer_gap_end_pos <- quer_bed[SYRI_CHR == quer_cn, V3]
    #
    quer_el_start_pos <- tmp_syri_chr$V7
    quer_el_end_pos <- tmp_syri_chr$V8
    #
    quer_el_gap_dists <- calc_el_gap_dist(el_start_pos = quer_el_start_pos,
                         el_end_pos = quer_el_end_pos,
                         gap_start_pos = quer_gap_start_pos,
                         gap_end_pos = quer_gap_end_pos)
    #
    tmp_syri_chr[, `:=` (REF_START_DIST_TO_GAP = ref_el_gap_dists[[1]],
                             REF_END_DIST_TO_GAP = ref_el_gap_dists[[2]],
                             QUER_START_DIST_TO_GAP = quer_el_gap_dists[[1]],
                             QUER_END_DIST_TO_GAP = quer_el_gap_dists[[2]])]
    #
    el_dist_out_list[[cn]] <- tmp_syri_chr
    }
  el_dist_out_dt <- rbindlist(el_dist_out_list)
  #
  el_dist_out_dt[, REF_MIN_GAP_DIST := 
    apply(el_dist_out_dt[, c('REF_START_DIST_TO_GAP', 'REF_END_DIST_TO_GAP')], 
          1, min)]
  el_dist_out_dt[, QUER_MIN_GAP_DIST :=
    apply(el_dist_out_dt[, c('QUER_START_DIST_TO_GAP', 'QUER_END_DIST_TO_GAP')],
          1, min)]
  #
  return(el_dist_out_dt)
}

######

