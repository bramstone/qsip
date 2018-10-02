#' Calculation of population dynamics
#'
#' Calculates population growth and death rates
#'
#' @param data Data as a \code{phyloseq} object
#'
#' @details Some details about proper isotope control-treatment factoring and timepoint specification. If weighted average densities or the
#'   change in weighted average densities have not been calculated beforehand, \code{calc_pop} will compute those first.
#'
#' @return \code{calc_pop} adds two S4 Matrix class objects (which more efficiently stores sparse matrix data) to the \code{data@@qsip@@.Data} slot
#'   of population birth rates for each taxon at each group of replicates. The row and column specifications will mirror those of the \code{phylosip}'s
#'   \code{\link{otu_table}}, meaning if taxa are listed on the table rows, they will in the resulting S4 Matrix class.
#'
#' @seealso \code{\link{calc_mw}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate population fluxes
#'
#' @export

calc_pop <- function(data) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  if(length(data@qsip@timepoint)==0) stop('Must specify different sample times with timepoint')
  if(data@qsip@iso!='18O') stop('Must use 18O-labeled treatment to calculate population flux')
  # if MW values don't exist, calculate those first
  # this will also handle rep_id validity (through calc_wad) and rep_group/iso_trt validity (through calc_d_wad)
  if(is.null(data@qsip[['mw_label']]) || is.null(data@qsip[['mw_light']])) data <- calc_mw(data)
  # transform sequencing abundances to 16S copy numbers
  # returns matrix with taxa as columns, samples as rows
  ft <- copy_no(data)
  # extract MW-labeled and convert to S3 matrix with taxa as ROWS (opposite all other calcs)
  mw_lab <- data@qsip[['mw_label']]
  mw_lab <- as(mw_lab, 'matrix')
  mw_l <- data@qsip[['mw_light']]
  if(!phyloseq::taxa_are_rows(data)) mw_lab <- t(mw_lab)
  # calculate mol. weight heavy max (i.e., what is maximum possible labeling)
  mw_max <- (12.07747 * 0.6) + mw_l
  # calculate copy number in light fraction at time t
  n_l <- ((mw_max - mw_lab)/(mw_max - mw_l)) * n_tot
  # calculate birth rate
  b <- log(ft_t / n_l) / t
  # calculate death rate
  d <- log(n_l / ft_0) / t
  # organize and add new data as S4 matrices
  data <- collate_results(data, b, 'pop_birth', sparse=TRUE)
  data <- collate_results(data, d, 'pop_death', sparse=TRUE)
  data <- collate_results(data, b - d, 'pop_flux', sparse=TRUE)
  return(data)
}
