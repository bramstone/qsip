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
  # calculate total 16S copy abundance for each sample, cbind so that taxa are rows
  ft <- split_data(data, ft, data@qsip@rep_id)
  ft <- lapply(ft, colSums, na.rm=T)
  ft <- do.call(rbind, ft)
  # separate samples based on timepoint
  time_group <- time_grouping(data,
                              as.numeric(as.character(data@qsip@timepoint)),
                              data@qsip@rep_id,
                              data@qsip@rep_group)
  # Drop any rows (probably NA) that don't appear in ft colnames, also drop any columns with NA for timepoint
  keep_cols <- (time_group$replicate %in% colnames(ft) & !is.na(time_group$time))
  if(sum(!keep_cols) > 0) {
    warning('Dropping group(s) ',
            paste(as.character(time_group$replicate[!keep_cols]), collapse=', '),
            ' from calculation', call.=FALSE)
  }
  time_group <- time_group[time_group$replicate %in% colnames(ft),]
  ft <- ft[!is.na(time_group$time),]
  time_group <- time_group[!is.na(time_group$time),]
  # split and sum 16S copy no.s across replicate groups, cbind so that taxa are rows
  ft <- split_data(data, ft, time_group$interaction, grouping_w_phylosip=F)
  ft <- lapply(ft, colSums, na.rm=T)
  ft <- do.call(cbind, ft)
  # get 16S copy numbers for different timepoints
  time_group2 <- unique(time_group[,!names(time_group) %in% 'replicate']) # only get unique elements to match levels in ft
  ft <- ft[,match(time_group2$interaction, colnames(ft))] # re-order columns to match time_group2$interaction
  timepoint_names <- paste0('n_t_',unique(time_group2$time))
  if(!any(time_group2$time==0)) {
    warning('No timepoints designated as time 0; using lowest timepoint as time 0', call.=FALSE)
  }
  for(i in 1:length(unique(time_group2$time))) {
    assign(timepoint_names[i], ft[,time_group2$time==(i-1)])
  }
  # extract MW-labeled and convert to S3 matrix with taxa as ROWS (opposite all other calcs)
  mw_lab <- data@qsip[['mw_label']]
  mw_lab <- as(mw_lab, 'matrix')
  mw_l <- data@qsip[['mw_light']]
  if(!phyloseq::taxa_are_rows(data)) mw_lab <- t(mw_lab)
  # calculate mol. weight heavy max (i.e., what is maximum possible labeling)
  mw_max <- (12.07747 * 0.6) + mw_l
  # calculate proportion in light fraction at time t
  # NEED FOR LOOP HERE TO GET N_light at each timepoint
  n_l_t <- ((mw_max - mw_lab)/(mw_max - mw_l)) * n_t_t
  # calculate birth rate
  b <- log(ft_t / n_l) / t
  # calculate death rate
  d <- log(n_l / ft_0) / t
  # organize and add new data as S4 matrices
  data <- collate_results(data, b, 'pop_birth', sparse=T)
  data <- collate_results(data, d, 'pop_death', sparse=T)
  data <- collate_results(data, b - d, 'pop_flux', sparse=T)
  return(data)
}
