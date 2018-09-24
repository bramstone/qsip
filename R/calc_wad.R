#' Calculation of weighted average densities
#'
#' Calculates weighted average densities for each microbial taxa in each sample replicate
#'
#' @param data Data as a \code{phyloseq} object
# @param ... Arguments to specify \code{NA} handling in the weighted average density calculation.
#   By default, \code{na.rm=FALSE}.
#'
#' @details Specifying \code{na.rm=TRUE} will allow \code{calc_wad} to calculate weighted average density values from samples
#'   that have one or more fractions without a valid density value. The default setting, \code{na.rm=FALSE}, returns values
#'   of \code{NA} for every taxa in a sample with missing density data.
#'
#' @return \code{calc_wad} adds an S4 Matrix class (which more efficiently stores sparse matrix data) to the \code{.Data} slot within
#'   the \code{qsip} slot of the object. of weighted average density values for each taxon at each sample. The row and column
#'   specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}}, meaning if taxa are listed on the table rows,
#'   they will in the resulting S4 Matrix class.
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate weighted average densities
#'
#' @export

calc_wad <- function(data) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  if(length(data@qsip@rep_id)==0) stop('Must specify replicate IDs')
  rep_id <- data@sam_data[[data@qsip@rep_id]] # replicate IDs
  dv <- data@sam_data[[data@qsip@density]] # density values
  # transform sequencing abundances to 16S copy numbers
  # returns matrix with taxa as columns, samples as rows
  ft <- copy_no(data)
  # split feature table by replicate ID
  # the split matrix will be roughly 2x the size of the phylosip object
  dv <- split(dv, rep_id)
  ft <- split(ft, rep_id)
  ft <- base::lapply(ft, matrix,
               byrow=FALSE,
               ncol=phyloseq::ntaxa(data))
  # next calculate WAD using Map (mapply) and apply by columns
  ft <- base::Map(function(y, x) apply(y, 2, wad, x, na.rm=TRUE), ft, dv)
  # combine format based on whether taxa were rows or not
  if(phyloseq::taxa_are_rows(data)) {
    ft <- do.call(cbind, ft)
  } else ft <- do.call(rbind.ft)
  # add feature names back in (replicate names automatically utilized from split)
  if(is.null(rownames(ft))) {
    rownames(ft) <- phyloseq::taxa_names(data)
  } else {
    colnames(ft) <- phyloseq::taxa_names(data)
  }
  # convert to S4 Matrix which is more memory efficient
  # ft <- Matrix::Matrix(ft)
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

wad <- function(y, x, na.rm=TRUE){
  wad <- sum(x * (y / sum(y, na.rm=TRUE)), na.rm=na.rm)
  return(wad)
}

# Calculation of 16S copies per taxon from relative abundances and 16S copy numbers
# Note, will return a feature table with taxa as columns
copy_no <- function(data) {
  ft <- as(data@otu_table, 'matrix')
  # make relative abundances with taxa as columns and samples as rows
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  ft <- ft / base::colSums(ft, na.rm=TRUE)
  # multiply by total 16S copy number per sample
  ft <- ft * data@sam_data[[data@qsip@abund]]
  return(ft)
}
