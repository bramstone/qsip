#' Calculation of weighted average densities
#'
#' Calculates weighted average densities for each microbial taxa in each sample replicate
#'
#' @param data Data as a \code{phylosip} object
#' @param filter Logical vector specifying whether or not to filter taxa from the weighted average density calculation.
#'   This will require \code{data} to have a filter applied with \code{\link{filter_qsip}}.
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
#' @seealso \code{\link{filter_qsip}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate weighted average densities
#'
#' @export

calc_wad <- function(data, filter=FALSE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  if(length(data@qsip@rep_id)==0) stop('Must specify replicate IDs with rep_id')
  # transform sequencing abundances to 16S copy numbers
  # returns matrix with taxa as columns, samples as rows
  ft <- copy_no(data)
  if(filter && length(data@qsip@filter)!=0) {
    ft <- ft[,colnames(ft) %in% data@qsip@filter]
  }
  # manipulate data matrix and calculate
  ft <- split_data(data, ft, data@qsip@rep_id) # split by replicate IDs
  dv <- split(data@sam_data[[data@qsip@density]],
              data@sam_data[[data@qsip@rep_id]]) # split densities by replicate IDs
  ft <- base::Map(function(y, x) apply(y, 2, wad, x, na.rm=TRUE), ft, dv)
  # organize and add new data as S4 matrix
  data <- collate_results(data, ft, 'wad')
  return(data)
}

# Utility functions to be used to do actual calculations

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

wad <- function(y, x, na.rm=TRUE) {
  wad <- sum(x * (y / sum(y, na.rm=TRUE)), na.rm=na.rm)
  return(wad)
}
