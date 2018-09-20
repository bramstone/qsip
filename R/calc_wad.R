# Calculation of weighted average densities

calc_wad <- function(data, transform_feature_abund=TRUE) {
  if(is(data)=='phylosip') stop('Must provide phylosip object')
  if(length(data@qsip@rep_id)==0) stop('Must specify replicate IDs')
    # maybe handle that by specifying calc as a generic with different methods
  # transform otu sequencing abundances to 16S copy numbers
  ot <- copy_no(data)
  # group and calculate by replicate ID not by taxa
  if(data@otu_table@OTUasRows==TRUE) {
    ot <- split(ot, data@sam_data@.Data[==data@qsip@rep_id])
    ot <- lapply(ot, matrix, nrow=nrow(data@otu_table), byrow=FALSE)
  } else {
    ot <- split(t(ot), data@sam_data@.Data[==data@qsip@rep_id])
    ot <- lapply(ot, matrix, ncol=ncol(data@otu_table), byrow=TRUE)
  }
  # next calculate WAD using lapply
  # wad_rep_tax = sum_for_all_fractions( density_of_fraction * (16S_frac_tax / 16S_rep_tax)
  if(transform_feature_abund==TRUE) data@otu_table@.Data <- ot
  return(data)
}

# Calculation of 16S copies per taxon from relative abundances and 16S copy numbers
copy_no <- function(data) {
  ot <- as(data@otu_table, 'matrix')
  # make relative abundances
  if(data@otu_table@OTUasRows==TRUE) {
    ot <- ot / rowSums(ot)
  } else {
    ot <- ot / colSums(ot)
  }
  # multiply by total 16S copy number per sample
  ot <- ot * data@sam_data@.Data[==data@qsip@abund]
  return(ot)
}

### Given vectors of x values and y values, calculate the weighted-average of the x-values (e.g., the weighted average density (WAD))
#
#     output = WAD.func(y, x)
#
#     y: vector of y-values (e.g., number of 16S copies)
#     x: vector of x-values (e.g., density of DNA)
#     -------------------------------------------------------
#     output: weighted-average of the x-values (single value)
#     Written by Ben Koch & Natasja van Gestel

wad <- function(y, x){
  wad <- sum(x*(y/sum(y)))
  return(wad)
}
