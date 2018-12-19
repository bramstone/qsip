#' Calculation of molecular weights
#'
#' Calculates molecular weights of microbial taxa due to isotope incorporation
#'
#' @param data Data as a \code{phyloseq} object
#' @param filter Logical vector specifying whether or not to filter taxa from the weighted average density calculation.
#'   This will require \code{data} to have a filter applied with \code{\link{filter_qsip}}.
#' @param correction Logical value indicating whether or not to apply tube-level correction to labeled WAD values.
#' @param offset_taxa Value from 0 to 1 indicating the percentage of the taxa to utilize for calculating offset correction values.
#'   Taxa are ordered by lowest difference in WAD values.
#'   Default is \code{0.1} indicating 10 percent of taxa with the lowest difference in WAD values.
#' @param separate_light Logical value indicating whether or not WAD-light scores should be averaged across all replicate groups or not.
#'   If \code{FALSE}, unlabeled WAD scores across all replicate groups will be averaged, creating a single molecular weight score per taxon
#'   representing it's genetic molecular weight in the absence of isotope addition.
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

calc_mw <- function(data, filter=FALSE, correction=FALSE, offset_taxa=0.1, separate_light=FALSE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  # if delta-WAD values don't exist, calculate those first
  # this will also handle rep_id validity (through calc_wad) and rep_group/iso_trt validity (through calc_d_wad)
  if(is.null(data@qsip[['wad_light']]) || is.null(data@qsip[['wad_label']])) data <- calc_d_wad(data, filter=filter, separate_light=separate_light,
                                                                                                correction=correction, offset_taxa=offset_taxa)
  # extract WAD-heavy / WAD-light values and convert to S3 matrices
  wh <- data@qsip[['wad_label']]
  wh <- as(wh, 'matrix')
  wl <- data@qsip[['wad_light']]
  wl <- as(wl, 'matrix')
  if(phyloseq::taxa_are_rows(data)) {wh <- t(wh); wl <- t(wl)}
  tax_names <- colnames(wh)
  # calculate GC content of each taxa (averaged across all groups of samples or not)
  gc <- (1 / 0.083506) * (wl - 1.646057)
  # calculate mol. weight of taxa without isotope
  mw_l <- (0.496 * gc) + 307.691
  # calculate mol. weight of taxa in labeled treatments
  if(all(dim(wh)==dim(wl))) {
    mw_h <- (((wh - wl)/wl) + 1) * mw_l
  } else {
    mw_h <- sweep(wh, 2, wl, function(x, y) (((x - y)/y) + 1))
    mw_h <- sweep(mw_h, 2, mw_l, '*')
  }
  if(is.null(dim(data@qsip[['wad_light']]))) mw_l <- c(mw_l)
  # organize and add new data as S4 matrices
  data <- collate_results(data, mw_h, tax_names=tax_names, 'mw_label', sparse=TRUE)
  data <- collate_results(data, mw_l, tax_names=tax_names, 'mw_light', sparse=TRUE)
  return(data)
}
