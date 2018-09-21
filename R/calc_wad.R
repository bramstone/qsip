#' Calculation of weighted average densities
#' @export

calc_wad <- function(data, ...) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  if(length(data@qsip@rep_id)==0) stop('Must specify replicate IDs')
  rep_id <- data@sam_data[[data@qsip@rep_id]] # replicate IDs
  dv <- data@sam_data[[data@qsip@density]] # density values
  # transform feature sequencing abundances to 16S copy numbers
  # x = each row of otu_table
  # y = vector of 16S copy number values
  data <- phyloseq::transform_sample_counts(data,
                                            function(x, y) (x / sum(x)) * y,
                                            y=data@sam_data[[data@qsip@abund]])
  ft <- as(data@otu_table, 'matrix')
  # split feature table by replicate ID
  if(data@otu_table@taxa_are_rows==TRUE) ft <- t(ft)
  dv <- split(dv, rep_id)
  ft <- split(ft, rep_id)
  ft <- lapply(ft, matrix,
               byrow=FALSE,
               ncol=phyloseq::ntaxa(data))
  # next calculate WAD using lapply and combine
  ft <- Map(function(y, x) apply(y, 2, wad, x, ...), ft, dv)
  # combine format based on whether taxa were rows or not
  if(data@otu_table@taxa_are_rows==TRUE) {
    ft <- do.call(cbind, ft)
  } else ft <- do.call(rbind.ft)
  # add feature names back in (replicate names automatically utilized from split)
  if(is.null(rownames(ft))) {
    rownames(ft) <- phyloseq::taxa_names(data)
  } else {
    colnames(ft) <- phyloseq::taxa_names(data)
  }
  # convert to S4 Matrix which is more memory efficient
  ft <- Matrix::Matrix(ft)
  # add wad values to data slot of qSIP portion of object
  if(!is.null(data@qsip@.Data$wad)) warning('Overwriting existing weighted average density values')
  data@qsip@.Data$wad <- ft
  return(data)
}

# Given vectors of x values and y values, calculate the weighted-average of the x-values (e.g., the weighted average density (WAD))
# With na.rm=F, any samples with NA density values will be returned with all taxa's WAD values as NA. (the default action)
# With na.rm=T, all fractions with NA densities will be ignored from a taxa's weighted average density (leading to lower numbers)
#
#
#     output = WAD.func(y, x)
#
#     y: vector of y-values (e.g., number of 16S copies)
#     x: vector of x-values (e.g., density of DNA)
#     -------------------------------------------------------
#     output: weighted-average of the x-values (single value)
#     Written by Ben Koch & Natasja van Gestel

wad <- function(y, x, ...){
  wad <- sum(x * (y / sum(y, na.rm=T)), ...)
  return(wad)
}
