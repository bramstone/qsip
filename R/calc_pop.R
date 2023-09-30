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
#'  If supplied as a factor, \code{calc_excess} will take the lowest level as the "light" treatment and the higher
#'  level as the "heavy" treatment.
#'  Alternatively, if supplied as a character, \code{calc_excess} will coerce the column to a factor and with the default behavior wherein the
#'  first value in alphabetical order will be assumed to be the lowest factor level (i.e. the "light" treatment).
#' @param isotope Column name specifying the isotope applied to each replicate. For "heavy" samples, these values should be one of "13C", "15N", or "18O".
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
#' @seealso \code{\link{calc_wad}}
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
                        correction = TRUE, rm_outliers = TRUE, non_grower_prop = 0.1, mu = 0.6,
                       nat_abund_13C = 0.01111233, nat_abund_15N = 0.003663004, nat_abund_18O = 0.002011429) {
  if(is.null(timepoint) || is.null(abund)) stop("Must specify timepoint and abundances")
  t0 <- data[timepoint == min(timepoint)]
  data <- wad_wide(data, tax_id = tax_id, sample_id = sample_id, wads = wads, iso_trt = iso_trt, isotope = isotope)
  #
  # correct density shifts
  if(correction) {
    label_shift <- data[, .(wad_diff = label - light), by = sample_id
                        ][!is.na(wad_diff)
                          ][order(wad_diff)
                            ][, .(shift = median(wad_diff[1:floor(non_grower_prop * .N)])), by = sample_id]
    data <- merge(data, label_shift, by = 'sample_id', all.x = TRUE)
    data[, label := label - shift][, shift := NULL]
  }
  # calculate molecular weights
  data[, gc_prop := (1 / 0.083506) * (light - 1.646057)
      ][, mw_light := (0.496 * gc_prop) + 307.691
       ][, `:=` (mw_label = (((label - light) / light) + 1) * mw_light,
                 mw_max = mw_light + 12.07747 * mu)]
  # calculate population rates
  data <- merge(data, t0, by = tax_id, all.x = TRUE)
  #
  dat[, abund_light := abund * ((mw_max - mw_label) / (mw_max - mw_light))
      ][, `:=` (growth = log(abund / abund_light) / timepoint,
                mortality = log(abund_light / abund_t0) / timepoint)]
  # clean final data output
  if(rm_outliers) {
    pos_out <- pos_outlier(data$eaf)
    neg_out <- neg_outlier(data$eaf)
  } else {
    pos_out <- Inf
    neg_out <- -Inf
  }
  # remove NA EAF values - this will also remove all the unlabeled samples
  data <- data[!is.na(growth) & !is.na(mortality),
               !c('wad_label', 'wad_light', 'gc_prop', 'mw_light', 'mw_label', 'mw_max', 'abund_16s_light')]
  return(data)
}
