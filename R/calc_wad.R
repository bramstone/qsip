#' Calculation of weighted average densities
#'
#' Calculates weighted average densities for each microbial taxa in each sample replicate
#'
#' @param data Data as a \code{phylosip} object
#'
#' @details Specifying \code{na.rm=TRUE} will allow \code{calc_wad} to calculate weighted average density values from samples
#'   that have one or more fractions without a valid density value. The default setting, \code{na.rm=FALSE}, returns values
#'   of \code{NA} for every taxa in a sample with missing density data.
#'
#'   The weighted average buoyant density for taxon's DNA \emph{i} in replicate \emph{j}, designated as \eqn{W_{ij}}, is:
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
#'  replicate-level abundances (\code{abund}). 
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

calc_wad <- function(data, sample_id = c(), fraction_id = c(), density = c(), abund = c(), grouping_cols = c()) {
  
  # calculate WADs
  dat <- dat[, abund_16s := rel_abund * avg_16S_g_soil, by = .(taxon_id, sample_fraction)
             ][, tot_abund := sum(abund_16s), by = .(taxon_id, sample_id) # group by replicate here, not sample-fraction
               ][, weight := abund_16s / tot_abund, by = .(taxon_id, sample_fraction)
                 ][, .(wad = sum(weight * Density.g.ml, na.rm = T),
                       wvd = sum(weight * (Density.g.ml - sum(weight * Density.g.ml, na.rm = T))^2, na.rm = T),   # weighted variance of density
                       abund_16s = sum(abund_16s, na.rm = T)),
                   by = c('taxon_id', 'sample_id', 'iso_trt', vars_to_keep, grouping_vars)
                   ][wad > 0]
    return(dat)
  }


