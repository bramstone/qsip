#' Calculation of molecular weights
#'
#' Calculates molecular weights of microbial taxa due to isotope incorporation
#'
#' @param data Data as a \code{phyloseq} object
#' @param separate_wad_light Logical value indicating whether or not WAD-light scores should be averaged across all replicate groups or not.
#'   If \code{FALSE}, WAD scores across all replicate groups will be averaged, creating a single WAD score per taxon representing it's weighted
#'   average density in the absence of isotope addition.
#' @param filter Logical vector specifying whether or not to filter taxa from the weighted average density calculation.
#'   This will require \code{data} to have a filter applied with \code{\link{filter_qsip}}.
#'
#' @details Some details about proper isotope control-treatment factoring. If weighted average densities or the change in weighted average densities
#'   have not been calculated beforehand, \code{calc_mw} will compute those first.
#'
#' @return \code{calc_mw} adds three S4 Matrix class objects (which more efficiently stores sparse matrix data) to the \code{data@@qsip@@.Data} slot
#'   of molecular weights for each taxon at each group of replicates in the labeled and unlabeled groups. The row and column
#'   specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}}, meaning if taxa are listed on the table rows,
#'   they will in the resulting S4 Matrix class.
#'
#' @seealso \code{\link{calc_wad}}, \code{\link{calc_d_wad}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate molecular weights
#'
#' @export

calc_mw <- function(data, separate_wad_light=TRUE, filter=FALSE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  # if delta-WAD values don't exist, calculate those first
  # this will also handle rep_id validity (through calc_wad) and rep_group/iso_trt validity (through calc_d_wad)
  if(is.null(data@qsip[['d_wad']])) data <- calc_d_wad(data, filter=filter)
  # extract d_WAD / WAD-light values and convert to S3 matrices
  ft <- data@qsip[['d_wad']]
  ft <- as(ft, 'matrix')
  wl <- data@qsip[['wad_light']]
  wl <- as(wl, 'matrix')
  if(phyloseq::taxa_are_rows(data)) {ft <- t(ft); wl <- t(wl)}
  tax_names <- colnames(ft)
  # calculate GC content of each taxa (averaged across all groups of samples or not)
  if(!separate_wad_light) {
    wl <- colMeans(wl, na.rm=T)
  }
  wl[is.nan(wl)] <- NA
  gc <- (1 / 0.083506) * (wl - 1.646057)
  # calculate mol. weight of taxa without isotope
  mw_l <- (0.496 * gc) + 307.691
  # calculate mol. weight of taxa in labeled treatments
  mw_lab <- ((ft/wl) + 1) * mw_l
  # organize and add new data as S4 matrices
  data <- collate_results(data, mw_lab, tax_names=tax_names, 'mw_label', sparse=TRUE)
  data <- collate_results(data, mw_l, tax_names=tax_names, 'mw_light', sparse=TRUE)
  return(data)
}
