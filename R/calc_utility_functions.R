# Collection of functions to make calculation function code neater

# function which prints sequence count and ASV count
# used in filtering taxa function
seq_check <- function(x) {
  nums <- formatC(c(x[, sum(seq_abund)], x[, uniqueN(taxon_id)]), format = 'fg', width = 12, big.mark = ',')
  cat(nums[1], ' sequences\n', nums[2], ' ASVs', sep = ''); cat('\n')
}


# Calculation of 16S copies per taxon from relative abundances and 16S copy numbers
# Note, will return a feature table with taxa as columns
copy_no <- function(data) {
  ft <- as(data@otu_table, 'matrix')
  # make relative abundances with taxa as columns and samples as rows
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  ft <- ft / base::rowSums(ft, na.rm=TRUE)
  # multiply by total 16S copy number per sample
  ft <- ft * data@sam_data[[data@qsip@abund]]
  ft[is.nan(ft)] <- 0
  return(ft)
}
