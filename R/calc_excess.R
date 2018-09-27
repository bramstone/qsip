#' Calculation of atom excess
#'
#' Calculates isotope incorporation in excess of natural abundances
#'
#' @param data Data as a \code{phyloseq} object
#' @param percent Logical value indicating whether or not to calculate atom percent excess (\code{percent=TRUE}) or atom excess fraction (the default)
#'
#' @details Some details about proper isotope control-treatment factoring. If weighted average densities or the change in weighted average densities
#'   have not been calculated beforehand, \code{calc_mw} will compute those first.
#'
#' @return \code{calc_excess} adds an S4 Matrix class objects (which more efficiently stores sparse matrix data) to the \code{data@@qsip@@.Data} slot
#'   of molecular weights for each taxon at each group of replicates in the labeled and unlabeled groups. The row and column
#'   specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}}, meaning if taxa are listed on the table rows,
#'   they will in the resulting S4 Matrix class.
#'
#' @seealso \code{\link{calc_wad}}, \code{\link{calc_d_wad}}, \code{\link{calc_mw}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate atom fraction excess
#'
#' @export

calc_excess <- function(data, percent=FALSE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  # if delta-WAD values don't exist, calculate those first
  # this will also handle rep_id validity (through calc_wad) and rep_group/iso_trt validity (through calc_d_wad)
  if(is.null(data@qsip[['mw_label']] || is.null(data@qsip[['mw_light']])) data <- calc_mw(data)
  # extract MW-labeled / MW-light and convert to S3 matrices
  mw_lab <- data@qsip[['mw_label']]
  mw_lab <- as(mw_lab, 'matrix')
  mw_l <- data@qsip[['mw_light']]
  mw_l <- as(mw_l, 'matrix)')
  if(phyloseq::taxa_are_rows(data)) mw_lab <- t(mw_lab); mw_l <- t(mw_l)
  # calculate mol. weight heavy max (i.e., what is maximum possible labeling)
  # 13C labeling must account for GC content because a GC bp has 1 less C atom than an AT bp
  if(iso=='18O') adjust <- 12.07747 + mw_l; nat_abund <- 0.002000429
  else if(iso=='13C') adjust <- (-0.4987282 * gc) + 9.974564; nat_abund <- 0.01111233
  mw_max <- adjust + mw_l
  # calculate atom excess
  excess <- ((mw_lab - mw_l)/(mw_max - mw_l)) * (1- nat_abund)
  # organize and add new data as S4 matrices
  if(percent) #...
  data <- collate_results(data, excess, 'atom_excess')
  return(data)
}
