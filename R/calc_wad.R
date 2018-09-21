# Calculation of weighted average densities

calc_wad <- function(data) {
  if(is(data)=='phylosip') stop('Must provide phylosip object')
  if(length(data@qsip@rep_id)==0) stop('Must specify replicate IDs')
  # transform feature sequencing abundances to 16S copy numbers
  # x = each row of otu_table
  # y = vector of 16S copy number values
  data <- phyloseq::transform_sample_counts(data,
                                            function(x, y) (x / sum(x)) * y,
                                            y=data@sam_data[[data@qsip@abund]])
  ft <- as(data@otu_table, 'matrix')
  # split feature table by replicate ID
  if(data@otu_table@taxa_are_rows==TRUE) ft <- t(ft)
  ft <- split(ft, data@sam_data[[data@qsip@rep_id]])
  ft <- lapply(ft, matrix,
               byrow=FALSE,
               ncol=phyloseq::ntaxa(data))
  # next calculate WAD using lapply and combine
  density <- split(data@sam_data[[data@qsip@density]], data@sam_data[[data@qsip@rep_id]])
  ft <- lapply(ft, function(y) apply(y, 1, wad, density))
  ft <- do.call(cbind, ft)
  # rename feature names and replicate ID names
  rownames(ft) <- phyloseq::taxa_names(data)
  colnames(ft) <- data@sam_data[[data@qsip@rep_id]]
  # convert to S4 Matrix which is more memory efficient
  ft <- Matrix::Matrix(ft)
  # add wad values to data slot of qSIP portion of object
  if(!is.null(data@qsip@.Data$wad)) warning('Overwriting existing weighted average density values')
  data@qsip@.Data$wad <- ft
  return(data)
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
  wad <- sum(x * (y / sum(y, na.rm=T)))
  return(wad)
}
