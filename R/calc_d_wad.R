#' Calculation of weighted average density differences
#'
#' Calculates differences in weighted average densities across replicates due to isotope incorporation
#'
#' @param data Data as a \code{phyloseq} object
#'
#' @details Some details about proper isotope control-treatment factoring. If weighted average densities have not been calculated
#'   beforehand, \code{calc_d_wad} will compute those first.
#'
#' @return \code{calc_d_wad} adds an S4 Matrix class (which more efficiently stores sparse matrix data) to the \code{.Data} slot within
#'   the \code{qsip} slot of the object of differences in weighted average density values for each taxon at each group of replicates. The row and column
#'   specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}}, meaning if taxa are listed on the table rows,
#'   they will in the resulting S4 Matrix class.
#'
#' @seealso \code{\link{calc_wad}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate weighted average density differences
#'
#' @export

calc_d_wad <- function(data) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  # if WAD values don't exist, calculate those first, will also handle rep_id validity
  if(is.null(data@qsip[['wad']])) data <- calc_wad(data)
  if(length(data@qsip@rep_group)==0) stop('Must specify replicate groupings with rep_group')
  if(length(data@qsip@iso)==0) stop('Must specify treatment and controls with iso')
  # extract WAD values and convert to S3 matrix
  ft <- data@qsip@.Data[['wad']]
  ft <- as(ft, 'matrix')
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  # manipulate data matrix and calculate
  ft <- split_data(data, ft, data@qsip@rep_id) # split by replicate IDs
  # combine isotope control/treatment factor with replicate IDs so we know what is light and heavy
  # ...
  ft <- base::Map(function(__) apply(__, 2, __, __), ft, iso_compare)
  # organize and add new data as S4 matrix
  data <- collate_results(data, ft, 'wad')
  return(data)
}
