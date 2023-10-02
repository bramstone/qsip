#' Calculation of excess atom fraction
#'
#' Calculates fractional isotope incorporation in excess of natural abundances
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
#' @details \code{calc_excess} automatically averages the isotopically unamended WAD values for each taxonomic feature on the assumption that density values
#'   will be identical (or nearly identical) for those samples.
#'
#'   The equations for calculating the molecular weights of taxon \emph{i}, designated \eqn{M_{Lab,i}} for labeled and \eqn{M_{Light,i}} for
#'   unlabeled, are:
#'
#'   \deqn{M_{Light,i} = 0.496 \cdot G_{i} + 307.691}
#'   \deqn{M_{Lab,i} = \left( \frac{\Delta W}{W_{Light,i}} + 1 \right) \cdot M_{Light,i}}
#'
#'   Where
#'
#'   \deqn{G_{i} = \frac{1}{0.083506} \cdot (W_{Light,i} - 1.646057)}
#'   Which indicates the GC content of taxon \emph{i} based on the density of its DNA when unlabeled
#'
#'   The calculation for the fractional of enrichment of taxon \emph{i}, \eqn{A_{i}} is:
#'
#'   \deqn{A_{i} = \frac{M_{Lab,i} - M_{Light,i}}{M_{Heavymax,i} - M_{Light,i}} \cdot (1 - N_{x})}
#'
#'   Where
#'
#'   \eqn{N_{x}}: The natural abundance of heavy isotope. Default estimates are: \eqn{N_{18O} = 0.002000429}, \eqn{N_{13C} = 0.01111233},
#'   and \eqn{N_{15N} = 0.003663004}
#'
#'   \eqn{N_{Heavymax,i}}: The highest theoretical molecular weight of taxon \emph{i} assuming maximum labeling by the heavy isotope
#'
#'   \eqn{N_{O,Heavymax,i} = (12.07747 + M_{Light,i}) \cdot L}
#'
#'   \eqn{N_{C,Heavymax,i} = (-0.4987282 \cdot G_{i} + 9.974564 + M_{Light,i}) \cdot L}
#'
#'   \eqn{N_{N,Heavymax,i} =( 0.5024851 \cdot G_{i} + 3.517396 + M_{Light,i}) \cdot L}
#'
#'   \emph{L}: The maximum label possible based off the percent of heavy isotope making up the atoms of that element in the labeled treatment
#'
#' @return \code{calc_excess} returns a data.table where each row represents a taxonomic feature within a single replicate.
#'  The following columns are produced: excess atom fraction (\code{eaf}). \code{NA} values (usually when an organism is present in
#'  only the unlabeled or labeled samples) are removed.
#'
#'  Because fraction-level data are being condensed to replicate-level, a list of columns to keep is not necessary.
#'
#'
#' @seealso \code{\link{calc_wad}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate excess atom fraction
#'
#'  # compare
#'
#' @references
#'  Hungate, Bruce, \emph{et al.} 2015. Quantitative microbial ecology through stable isotope probing.
#'  \emph{Applied and Environmental Microbiology} \strong{81}: 7570 - 7581.
#'
#'  Morrissey, Ember, \emph{et al.} 2018. Taxonomic patterns in nitrogen assimilation of soil prokaryotes.
#'  \emph{Environmental Microbiology} \strong{20}: 1112 - 1119.
#'
#' @export

calc_excess <- function(data, tax_id = c(), sample_id = c(), wads = 'wad',
                        iso_trt = c(), isotope = c(),
                        correction = TRUE, rm_outliers = TRUE, non_grower_prop = 0.1,
                        nat_abund_13C = 0.01111233, nat_abund_15N = 0.003663004, nat_abund_18O = 0.002011429) {
  data <- wad_wide(data, tax_id = tax_id, sample_id = sample_id, wads = wads, iso_trt = iso_trt, isotope = isotope)
  # calculate molecular weights
  data[, gc_prop := (1 / 0.083506) * (wad_light - 1.646057)
       ][, mw_light := (0.496 * gc_prop) + 307.691
         ][, mw_label := (((wad_label - wad_light) / wad_light) + 1) * mw_light]
  # calculate enrichment
  data[isotope == '18O', `:=` (mw_max = mw_light + 12.07747, nat_abund = nat_abund_18O)
       ][isotope == '13C', `:=` (mw_max = mw_light + 9.974564 + (-0.4987282 * gc_prop), nat_abund = nat_abund_13C)
         ][isotope == '15N', `:=` (mw_max = mw_light + 3.517396 + (0.5024851 * gc_prop), nat_abund = nat_abund_15N)
           ][, eaf := ((mw_label - mw_light) / (mw_max - mw_light)) * (1 - nat_abund)]
  # correct enrichment values
  if(correction) {
    if(rm_outliers) {
      pos_out <- pos_outlier(data$eaf)
      neg_out <- neg_outlier(data$eaf)
    } else {
      pos_out <- Inf
      neg_out <- -Inf
    }
    shift <- data[!is.na(eaf)
                  ][order(eaf)
                    ][eaf > neg_out
                      ][eaf < pos_out
                        ][, .(shift = median(eaf[1:floor(non_grower_prop * .N)])),
                          by = sample_id]
    data <- merge(data, shift, by = 'sample_id', all.x = TRUE)
    data[, eaf_corrected := eaf - shift][, shift := NULL]
  }
  # clean final data output
  # remove NA EAF values - this will also remove all the unlabeled samples
  eaf_dat <- data[!is.na(eaf), !c('wad_label', 'wad_light', 'wvd_light', 'gc_prop',
                                  'mw_light', 'mw_label', 'mw_max', 'nat_abund', 'tube_shift')]
  setnames(eaf_dat, 'wvd_label', 'wvd')
  return(data)
}
