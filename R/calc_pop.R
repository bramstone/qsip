#' Calculation of population dynamics
#'
#' Calculates population growth and death rates
#'
#' @param data Data as a long-format data.table where each row represents a taxonomic feature within a single fraction.
#'  Typically, this is the output from the \code{calc_wad} function.
#' @param tax_id Column name specifying unique identifier for each taxonomic feature.
#' @param sample_id Column name specifying unique identifier for each replicate.
#' @param wads Column name specifying weighted average density values.
#' @param iso_trt Column name specifying a two-level categorical column indicating whether a sample has been amended with a stable isotope (i.e., is "heavy") or if
#'  isotopic composition is at natural abundance (i.e., "light").
#'  Any terms may be applied but care should be taken for these values.
#'  If supplied as a factor, \code{calc_pop} will take the lowest level as the "light" treatment and the higher
#'  level as the "heavy" treatment.
#'  Alternatively, if supplied as a character, \code{calc_pop} will coerce the column to a factor and with the default behavior wherein the
#'  first value in alphabetical order will be assumed to be the lowest factor level (i.e. the "light" treatment).
#' @param timepoint Column name specifying the timepoint at which each sample was collected.
#'  For population rates, the lowest timepoint will be assumed to be the initial timepoints with which to base
#'  changes in abundances on.
#' @param abund Column name specifying abundance measurement for each fraction. Typically either the numbers of a target gene amplicon
#'  (e.g., 16S, ITS such as from qPCR) or DNA concentration (e.g., nanograms per microliter such as from a Qubit).
#' @param correction Whether to apply a correction to fractional enrichment values to ensure a certain proportion are positive.
#' @param rm_outlers Whether or not to remove fractional enrichment values that are 1.5X greater or lesser than the distance between the median
#'   and interquartile ranges.
#' @param non_grower_group Fractional value applied if \code{correction == TRUE} specifying the proportion of the community in each samples assumed to be
#'  non-growers and whose median enrichment values will be assumed to be zero. The adjustment necessary to place this median value at zero will be applied
#'  as a correction to all enrichment values in the sample.
#'
#' @details Details to go here.
#'
#' @return \code{calc_pop} returns a data.table where each row represents a taxonomic feature within a single replicate.
#'  The following columns are produced: growth rate (\code{growth}), mortality or turnover rate (\code{mortality}).
#'
#'
#' @seealso \code{\link{calc_wad}, \link{wad_wide}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate population fluxes
#'
#' @references
#'  Koch, Benjamin, \emph{et al.} 2018. Estimating taxon-specific population dynamics in diverse microbial communities.
#'  \emph{Ecosphere} \strong{9}.
#'
#' @export

calc_pop <- function(data, tax_id = c(), sample_id = c(), wads = 'wad',
                     iso_trt = c(), isotope = c(), timepoint = c(), abund = c(),
                     t0_grouping = c(), correction = TRUE, rm_outliers = TRUE,
                     non_grower_prop = 0.1, mu = 0.6) {
  if(is.null(timepoint) || is.null(abund)) stop("Must specify timepoint and abundances")
  # re-express the timepoint column as factor, assing lowest level as T0
  if(!is.factor(timepoint)) message('Assigned', levels(timepoint)[1], 'as the initial timepoint.')
  iso_trt <- as.factor(data$timepoint)
  # average abundances across initial timepoints
  t0 <- data[timepoint == levels(timepoint)[1]]
  if(is.null(t0_grouping)) {
    t0 <- t0[, .(abund_t0 = mean(abund)), by = tax_id]
  } else {
    t0 <- t0[, .(abund_t0 = mean(abund)), by = c(tax_id, t0_grouping)]
  }
  #
  data <- wad_wide(data[timepoint != levels(timepoint)[1]],
                   tax_id = tax_id, sample_id = sample_id, wads = wads,
                   iso_trt = iso_trt, isotope = isotope)
  #
  # correct density shifts
  if(correction) {
    label_shift <- data[, .(wad_diff = wad_label - wad_light), by = sample_id
                        ][!is.na(wad_diff)
                          ][order(wad_diff)
                            ][, .(shift = median(wad_diff[1:floor(non_grower_prop * .N)])), by = sample_id]
    data <- merge(data, label_shift, by = 'sample_id', all.x = TRUE)
    data[, wad_label := wad_label - shift][, shift := NULL]
  }
  # calculate molecular weights
  data[, gc_prop := (1 / 0.083506) * (light - 1.646057)
      ][, mw_light := (0.496 * gc_prop) + 307.691
       ][, `:=` (mw_label = (((wad_label - wad_light) / wad_light) + 1) * mw_light,
                 mw_max = mw_light + 12.07747 * mu)]
  # merge in t0 abundances
  if(is.null(t0_grouping)) {
    data <- merge(data, t0, by = tax_id, all.x = TRUE)
  } else {
    data <- merge(data, t0, by = c(tax_id, t0_grouping), all.x = TRUE)
  }
  # calculate population rates
  dat[, abund_light := abund * ((mw_max - mw_label) / (mw_max - mw_light))
      ][, `:=` (growth = log(abund / abund_light) / timepoint,
                mortality = log(abund_light / abund_t0) / timepoint)]
  # clean final data output
  if(rm_outliers) {
    r_pos_out <- pos_outlier(data$growth)
    r_neg_out <- neg_outlier(data$growth)
    m_pos_out <- pos_outlier(data$mortality)
    m_neg_out <- neg_outlier(data$mortality)
  } else {
    r_pos_out <- m_pos_out <- Inf
    r_neg_out <- m_neg_out <- -Inf
  }
  data <- data[growth > r_neg_out
               ][growth < r_pos_out
                 ][mortality > m_neg_out
                   ][mortaity < m_pos_out]
  # remove NA EAF values - this will also remove all the unlabeled samples
  data <- data[!is.na(growth) & !is.na(mortality),
               !c('wad_label', 'wad_light', 'gc_prop', 'mw_light', 'mw_label', 'mw_max', 'abund_light')]
  return(data)
}
