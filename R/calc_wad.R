#' Calculation of weighted average densities
#'
#' Calculates weighted average densities for each microbial taxa in each sample replicate
#'
#' @param data Data as a long-format data.table where each row represents a taxonomic feature within a single fraction.
#' @param tax_id Unique identifier for each taxonomic feature. Required.
#' @param sample_id Unique identifier for each replicate. Required
#' @param frac_id Fraction identifier. Does not have to be unique to each replciate because \code{calc_wad} will 
#'  combine the unique sample ID with the fraction ID to generate a unique sample-fraction code. Required
#' @param frac_dens Buoyant density value for each fraction. Typically expressed as grams per milliliter from a cesium chloride
#'  density buffer. Required
#' @param frac_abund Abundance measurement for each fraction. Typically either the numbers of a target gene amplicon
#'  (e.g., 16S, ITS such as from qPCR) or DNA concentration (e.g., nanograms per microliter such as from a Qubit). Required
#' @param rel_abund Relativized abundance of each taxonomic feature, typically calculated after the removal of non-target
#'  lineages but before frequency filtering has been applied. Required
# @param grouping_cols Additional columns that should be included as important treatment groups in the output.
#'  Not strictly necessary for the calculation, but these will be utilized next to calculate fractional isotopic enrichment.
#'  Taxonomic information may be included here as well.
#'  
#'
#' @details The weighted average buoyant density (WAD) for taxon's DNA \emph{i} in replicate \emph{j}, designated as \eqn{W_{ij}}, is:
#'
#'   \deqn{W_{ij} = \sum_{k=1}^{K} x_{jk} \cdot \left( \frac{y_{ijk}}{y_{ij}} \right)}
#'
#'   Where
#'   \eqn{y_{ij} = \sum_{k=1}^{K} y_{ijk}}: Total abundance of taxon \emph{i} in replicate \emph{j}, by summing across all fractions (\emph{K})
#'
#'   \eqn{y_{ijk} = P_{ijk} \cdot f_{jk}}: The total abundance per \eqn{\mu}L for taxon \emph{i} in density fraction
#'   \emph{k} of replicate \emph{i} as calculated by it's relative abundance, \eqn{P_{ijk}} in that fraction multiplied by the
#'   total abundance of DNA or specific amplicons in that fraction, \eqn{f_{jk}}
#'
#'   \eqn{x_{jk}}: Density of fraction \emph{k} of replicate \emph{j} in (g cm\eqn{^-3})
#'
#'
#' @return \code{calc_wad} returns a reduced data.table where each row represents a taxonomic feature within a single replicate.
#'  The following columns are produced: weighted average densities (\code{wad}), weighted variance of densities (\code{wvd}), and
#'  replicate-level abundances (\code{abund}). \code{wad} values < 0 are removed.
#'
#'  The \code{wvd} term produced is the weighted variance in density, essentially a measure of how spread-out a feature's WAD value is. 
#'  It is not strictly necessary to calculate fractional enrichment. 
#'  It can sometimes be useful as a diagnostic of the quality of a feature’s density estimate.
#'
#'  Abundances are expressed in the same unit as fraction-level abundance measures of the community. For example, if fraction-level
#'  abundances were made using qPCR of a target gene (e.g., 16S or ITS), abundances represent the proportion of that gene attributed
#'  to a given taxon. If fraction-level abundances were made using DNA concentrations (e.g., such as from a Qubit measure), then the
#'  \code{abund} column is an expression of a taxon's proportional contribution to DNA concentration.
#'
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate weighted average densities
#'
#' @references
#'  Hungate, Bruce, \emph{et al.} 2015. Quantitative microbial ecology through stable isotope probing.
#'  \emph{Applied and Environmental Microbiology} \strong{81}: 7570 - 7581.
#'
#' @export

calc_wad <- function(data, tax_id = c(), sample_id = c(), frac_id = c(), 
                     frac_dens = c(), frac_abund = c(), rel_abund = c(), 
                     grouping_cols = c()) {
  vars <- list(tax_id, sample_id, frac_id, frac_dens, frac_abund, rel_abund)
  if(any(sapply(vars, is.null))) {
    null_vars <- which(sapply(vars, is.null))
    null_vars <- paste(c('taxon IDs', 'sample IDs', 'fraction IDs', 'densities', 
                         'fraction abundances', 'relative taxon abundances')[null_vars],
                       sep = ',')
    stop("Must supply the following columns:", null_vars)
  }
  # use sample ID and fraction ID to create unique sample-fraction ID
  data[, sample_fraction := paste(sample_id, fraction_id, sep = '_')]
  # calculate WADs
  data <- data[, frac_abund_tax := rel_abund * frac_abund, by = .(taxon_id, sample_fraction)
               ][, tot_abund := sum(frac_abund_tax), by = .(taxon_id, sample_id) # group by replicate here, not sample-fraction
                 ][, weight := frac_abund_tax / tot_abund, by = .(taxon_id, sample_fraction)
                   ][, .(wad = sum(weight * density, na.rm = TRUE),
                         wvd = sum(weight * (density - sum(weight * density, na.rm = TRUE))^2, na.rm = TRUE),   # weighted variance of density
                         abund = sum(frac_abund_tax, na.rm = TRUE)),
                     by = c(taxon_id, sample_id, grouping_cols)
                     ][wad > 0]
    return(data)
  }


