#' Calculation of atom excess
#'
#' Calculates isotope incorporation in excess of natural abundances
#'
#' @param data Data as a \code{phyloseq} object
#' @param percent Logical value indicating whether or not to calculate atom percent excess (\code{percent=TRUE}) or atom excess fraction (the default)
#' @param ci_method Character value indicating how to calculate confidence intervals of stable isotope atom excess.
#'   Options are \code{bootstrap} or \code{bayesian} (see \code{details} below for discussion on their differences).
#'   The default is blank indicating that no confidence intervals will be calculated.
#' @param ci Numeric value from 0 to 1 indicating the width of the confidence interval for bootsrapped atom excess values.
#' @param iters Number of (subsampling) iterations to conduct to calculate confidence intervals. Default is \code{999}.
#' @param filter Logical vector specifying whether or not to filter taxa from the weighted average density calculation.
#'   This will require \code{data} to have a filter applied with \code{\link{filter_qsip}}.
#' @param correction Logical value indicating whether or not to apply tube-level correction to labeled WAD values.
#' @param offset_taxa Value from 0 to 1 indicating the percentage of the taxa to utilize for calculating offset correction values.
#'   Taxa are ordered by lowest difference in WAD values.
#'   Default is \code{0.1} indicating 10 percent of taxa with the lowest difference in WAD values.
#' @param max_label Numeric value indicating the maximum possible isotope labeling in an experiment.
#'   Keeping the value at \code{1} will ensure that the maximum possible atom excess value of 1 corresponds to complete updake of the isotope.
#'   Recommended for experiments with lower atom percent enrichment treatments (see Details).
#' @param separate_light Logical value indicating whether or not WAD-light scores should be averaged across all replicate groups or not.
#'   If \code{FALSE}, unlabeled WAD scores across all replicate groups will be averaged, creating a single molecular weight score per taxon
#'   representing it's genetic molecular weight in the absence of isotope addition.
#' @param separate_label Logical value indicating whether or not WAD-label scores should be averaged across all replicate groups or not.
#'   If \code{FALSE}, labeled WAD scores across all replicate groups will be averaged, creating a single molecular weight score per taxon
#'   representing it's genetic molecular weight as a result of isotope addition. The default is \code{TRUE}.
#' @param recalc Logical value indicating whether or not to recalculate WAD and molecular weight values or use existing values. Default is \code{TRUE}.
#'   Using bootstrapped calculations will automatically recalculate all values.
#'
#' @details Some details about proper isotope control-treatment factoring. If weighted average densities or the change in weighted average densities
#'   have not been calculated beforehand, \code{calc_mw} will compute those first.
#'
#'   Atom excess values calculated when \code{max_label < 1} are \strong{not} atom percent excess values but rather \emph{percent maximum enrichment}.
#'   Setting \code{max_label < 1} will return an atom excess value higher than would otherwise be returned. A \code{max_label} value of 1 indicates
#'   that atom excess fraction ranges from 0 to 1 where 0 indicates no isotope incorporation and 1 indicates complete, or 100\%, isotope incorporation.
#'   For various reasons, complete isotope incorporation will be impossible. However, a \code{max_label} value less than 1 will indicate atom excess
#'   fraction where 0 indicates no isotope incorporation and 1 indicates the highest possible incorporation, as constrained by atom percent enrichment
#'   provided in the experiment. For example, an experiment enriching soil with 13-C at 50\% atom enrichment will want to specify \code{max_label=0.5}
#'   and an atom excess fraction value of 1 in this case corresponds to an organism that has succeeded in incorporating 13-C into it's nucleic acids
#'   at 50\%.
#'
#'   The calculation for the proportion of enrichment of taxon \emph{i}, \eqn{A_{i}} is:
#'
#'   \deqn{A_{i} = \frac{M_{Lab,i} - M_{Light,i}}{M_{Heavymax,i} - M_{Light,i}} \cdot (1 - N_{x})}
#'
#'   Where
#'
#'   \eqn{N_{x}}: The natural abundance of heavy isotope. \eqn{N_{18O} = 0.002000429}, \eqn{N_{13C} = 0.01111233},
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
#' @return \code{calc_excess} adds an S4 Matrix class objects (which more efficiently stores sparse matrix data) to the \code{data@@qsip@@.Data} slot
#'   of molecular weights for each taxon at each group of replicates in the labeled and unlabeled groups. The row and column
#'   specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}}, meaning if taxa are listed on the table rows,
#'   they will in the resulting S4 Matrix class.
#'
#'   Note that the bootstrap method produces a \emph{single} bootstrapped median (and matching confidence intervals) for groups of replicates,
#'   either grouped by isotope treatment alone, or also by some other grouping factor (if \code{data@@qsip@@rep_group} is specified).
#'   Using no bootstrap value allows separate enrichment values to be attained for each replicate, if \code{separate_label=TRUE}.
#'
#' @seealso \code{\link{calc_wad}}, \code{\link{calc_d_wad}}, \code{\link{calc_mw}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate atom fraction excess
#'
#'  # compare
#'  all.equal(aef*100, ape)
#'
#' @references
#'  Hungate, Bruce, \emph{et al.} 2015. Quantitative microbial ecology through stable isotope probing.
#'  \emph{Applied and Environmental Microbiology} \strong{81}: 7570 - 7581.
#'
#'  Morrissey, Ember, \emph{et al.} 2018. Taxonomic patterns in nitrogen assimilation of soil prokaryotes.
#'  \emph{Environmental Microbiology} \strong{20}: 1112 - 1119.
#'
#' @export

calc_excess <- function(data, tax_id = c(), sample_id = c(), iso_trt = c(),
                        wad_light = c(), wad_label = c(), abund = c(), grouping_cols = c(),
                       correction = TRUE, rm_outliers = TRUE, non_grower_prop = 0.1,
                       nat_abund_13C = 0.01111233, nat_abund_15N = 0.003663004, nat_abund_18O = 0.002011429) {
  vars <- list(tax_id, sample_id, iso_trt, wad_light, wad_label))
  if(any(sapply(vars, is.null)) {
    null_vars <- which(sapply(vars, is.null))
    null_vars <- paste(c('taxon IDs', 'sample IDs', 'isotope addition (for each row "13C" or "15N" or "18O")',
                         'unlabeled WADs', 'labeled WADs')[null_vars],
                       sep = ',')
    stop("Must supply the following columns:", null_vars)
  }
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
