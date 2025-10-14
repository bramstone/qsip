# Collection of functions to make calculation function code neater

# function which prints sequence count and ASV count
# used in filtering taxa function
seq_summary <- function(x, reads = c(), tax_id = c()) {
  if(is.null(reads) || is.null(tax_id)) stop("Must supply columns for sequence reads and taxon ID")
  nums <- formatC(c(sum(x[[reads]], na.rm = TRUE),
                    uniqueN(x[[tax_id]])),
                  format = 'fg', width = 12, big.mark = ',')
  cat(nums[1], ' sequences\n', nums[2], ' ASVs', sep = ''); cat('\n')
}


# exclude outliers following the standard definition of above or below
# the 25th/75th quantile by a distance of 1.5X the interquartile range
neg_outlier <- function(x) {
  y <- x[!is.na(x)]
  if(length(y) == 0) warning("Outliers cannot be calculated on all-NA vector")
  neg_out <- quantile(y, .25) - (1.5 * IQR(y))
  suppressWarnings(max(y[y < neg_out]))
}

pos_outlier <- function(x) {
  y <- x[!is.na(x)]
  if(length(y) == 0) warning("Outliers cannot be calculated on all-NA vector")
  pos_out <- quantile(y, .75) + (1.5 * IQR(y))
  suppressWarnings(min(y[y > pos_out]))
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
