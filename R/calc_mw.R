#' Calculation of molecular weights
#'
#' Calculates differences in molecular weights of microbial taxa due to isotope incorporation
#'
#' @param data Data as a \code{phyloseq} object
#'
#' @details Some details about proper isotope control-treatment factoring. If weighted average densities or the change in weighted average densities
#'   have not been calculated beforehand, \code{calc_mw} will compute those first.
#'
#' @return \code{calc_mw} adds an S4 Matrix class (which more efficiently stores sparse matrix data) to the \code{data@@qsip@@.Data} slot
#'   of differences in molecular weights for each taxon at each group of replicates. The row and column
#'   specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}}, meaning if taxa are listed on the table rows,
#'   they will in the resulting S4 Matrix class.
#'
#' @seealso \code{\link{calc_wad}}, \code{\link{calc_d_wad}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate weighted average density differences
#'
#' @export

calc_mw <- function(data) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  # if delta-WAD values don't exist, calculate those first
  # this will also handle rep_id validity (through calc_wad) and rep_group/iso_trt validity (through calc_d_wad)
  if(is.null(data@qsip[['d_wad']])) data <- calc_d_wad(data)
  # extract d_WAD / WAD-light values and convert to S3 matrices
  ft <- data@qsip[['d_wad']]
  ft <- as(ft, 'matrix')
  wl <- data@qsip[['wad_light']]
  wl <- as(wl, 'matrix)')
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft); wl <- t(wl)
  # manipulate data matrix and calculate

  # organize and add new data as S4 matrix
  data <- collate_results(data, d_ft, 'd_wad')
  return(data)
}
