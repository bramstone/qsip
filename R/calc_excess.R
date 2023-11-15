#' Calculation of excess atom fraction
#'
#' Calculates fractional isotope incorporation in excess of natural abundances
#'
#' @param data Data as a long-format data.table where each row represents a taxonomic feature within a single fraction.
#'  Typically, this is the output from the \code{calc_wad} function.
#' @param tax_id Column name specifying unique identifier for each taxonomic feature.
#' @param sample_id Column name specifying unique identifier for each replicate.
#' @param wads Column name specifying weighted average density values.
#'  It is possible to change this if preferred but not recommended.
#' @param iso_trt Column name specifying a two-level categorical column indicating whether a sample has been amended with a stable isotope (i.e., is "heavy") or if
#'  isotopic composition is at natural abundance (i.e., "light").
#'  Any terms may be applied but care should be taken for these values.
#'  If supplied as a factor, \code{calc_excess} will take the lowest level as the "light" treatment and the higher
#'  level as the "heavy" treatment.
#'  Alternatively, if supplied as a character, \code{calc_excess} will coerce the column to a factor and with the default behavior wherein the
#'  first value in alphabetical order will be assumed to be the lowest factor level (i.e. the "light" treatment).
#' @param isotope Column name specifying the isotope applied to each replicate. For "heavy" samples, these values should be one of "13C", "15N", or "18O".
#' @param bootstrap Whether to generate bootstrapped enrichment values for each taxonomic feature across groups of samples.
#'  If \code{TRUE}, replicates within specified treatment groupings will be randomly resampled with replacement and resulting
#'  WAD values used to generate a distribution of enrichment values for each taxon.
#' @param iters Integer specifying the number of bootstrap iterations to perform.
#' @param grouping_cols Column(s) used to indicate bootstrap resampling groups.
#' Within each group, replicates will be resampled with replacement.
#' Resampling will not occur across groups.
#' @param min_freq Minimum number of replicates a taxonomic feature must occur in to be kept.
#'  If treatment grouping columns are specified, frequencies will be assessed at this level.
#'  For unlabeled "light" samples, treatment groupings will be ignored when assessing adequate frequency of occurrence.
#' @param correction Whether to apply a correction to fractional enrichment values to ensure a certain proportion are positive.
#' @param rm_outlers Whether or not to remove low fractional enrichment values that are 1.5X lesser than the distance between the 25th quantile
#'   and interquartile range.
#'   If \code{bootstrap = TRUE}, outlier WAD values will be removed prior to resampling and enrichment calculation.
#' @param non_grower_prop Fractional value applied if \code{correction == TRUE} specifying the proportion of the community in each samples assumed to be
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
#'  Because fraction-level data are being condensed to replicate-level, a list of grouping columns to keep is not necessary
#'  \emph{unless} bootstrapping is specified, wherein those groupings will be used to organize the resampling effort.
#'
#'  Bootstrap functionality produces the columns: \code{2.5\%}, \code{50\%}, and \code{97.5\%}
#'  reflecting the median bootstrapped enrichment value and the 95\% confidence interval
#'  of enrichment values.
#'  In addition, probability that a taxon's enrichment was less than zero is represented
#'  in the \code{p_val} column which is based on the proportion sub-zero bootstrapped observations.
#'
#'
#' @seealso \code{\link{calc_wad}, \link{wad_wide}}
#'
#' @examples
#'  # Load in example data
#'  data(example_qsip)
#'
#'  # relativize sequence abundances (should be done after taxonomic filtering)
#'  example_qsip[, rel_abund := seq_abund / sum(seq_abund), by = sampleID]
#'
#'  # calculate weighted average densities
#'  wads <- calc_wad(example_qsip,
#'                   tax_id = 'asv_id', sample_id = 'sampleID', frac_id = 'fraction',
#'                   frac_dens = 'Density.g.ml', frac_abund = 'avg_16S_g_soil',
#'                   rel_abund = 'rel_abund',
#'                   grouping_cols = c('treatment', 'isotope', 'iso_trt', 'Phylum'))
#'
#'  # calculate fractional enrichment in excess of background
#'  eaf <- calc_excess(wads, tax_id = 'asv_id', sample_id = 'sampleID',
#'                     iso_trt = 'iso_trt', isotope = 'isotope')
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
                        bootstrap = FALSE, iters = 999L, grouping_cols = c(), min_freq = 3,
                        correction = TRUE, rm_outliers = TRUE, non_grower_prop = 0.1,
                        nat_abund_13C = 0.01111233, nat_abund_15N = 0.003663004, nat_abund_18O = 0.002011429) {
  vars <- list(tax_id, sample_id, iso_trt, isotope)
  if(any(sapply(vars, is.null))) {
    null_vars <- which(sapply(vars, is.null))
    null_vars <- paste(c('taxon IDs', 'sample IDs', 'isotope addition (for each row "13C" or "15N" or "18O")',
                         'isotope treatment (amended or unamended)')[null_vars],
                       sep = ',')
    stop("Must supply the following columns: ", null_vars)
  }
  vars <- c(vars, 'wad')
  if(any(sapply(vars, function(x) !exists(x, data)))) {
    missing_vars <- which(sapply(vars, function(x) !exists(x, data)))
    missing_vars <- paste(vars[missing_vars], sep = ',')
    stop("Missing the following column(s) in supplied data: ", missing_vars)
  }
  if(bootstrap == FALSE) {
    eafd <- wad_wide(data, tax_id = tax_id, sample_id = sample_id, wads = wads, iso_trt = iso_trt, isotope = isotope)
    # calculate molecular weights
    eafd[, gc_prop := (1 / 0.083506) * (light - 1.646057)
         ][, mw_light := (0.496 * gc_prop) + 307.691
           ][, mw_label := (((label - light) / light) + 1) * mw_light]
    # calculate enrichment
    eafd[isotope == '18O', `:=` (mw_max = mw_light + 12.07747, nat_abund = nat_abund_18O)
         ][isotope == '13C', `:=` (mw_max = mw_light + 9.974564 + (-0.4987282 * gc_prop), nat_abund = nat_abund_13C)
           ][isotope == '15N', `:=` (mw_max = mw_light + 3.517396 + (0.5024851 * gc_prop), nat_abund = nat_abund_15N)
             ][, eaf := ((mw_label - mw_light) / (mw_max - mw_light)) * (1 - nat_abund)]
    # correct enrichment values
    if(rm_outliers == TRUE) {
      neg_out <- neg_outlier(eafd$eaf)
      pos_out <- pos_outlier(eafd$eaf)
    } else {
      neg_out <- -Inf
      pos_out <- Inf
    }
    eafd <- eafd[eaf > neg_out]
    if(correction) {
      shift <- eafd[!is.na(eaf)
                    ][order(eaf)
                      ][, .(shift = median(eaf[1:floor(non_grower_prop * .N)])),
                        by = sample_id]
      eafd <- merge(eafd, shift, by = sample_id, all.x = TRUE)
      eafd[, eaf := eaf - shift][, shift := NULL]
      if(rm_outliers) {
      } else {
        neg_out <- -Inf
      }
      # clean final data output
      # remove NA EAF values - this will also remove all the unlabeled samples
      eafd <- eafd[!is.na(eaf), !c('label', 'light', 'gc_prop',
                                   'mw_light', 'mw_label', 'mw_max', 'nat_abund')]
      setnames(eafd, 'wvd_label', 'wvd', skip_absent = TRUE)
    }
    ####################
    # BOOTSTRAPPING CODE
    ####################
  } else if(bootstrap == TRUE) {
    bd <- copy(data)
    setnames(bd, old = c(sample_id, iso_trt, wads), new = c('sid', 'iso_trt', 'wad'))
    # make sure iso_trt column is a factor
    if(is.factor(bd[[iso_trt]]) == FALSE || nlevels(bd[[iso_trt]]) > 2) {
      test_trt <- factor(bd[[iso_trt]])
      light_trt <- levels(test_trt)[1]
      message('Assigned ', light_trt, ' as the unamended or "light" treatment and ',
              levels(test_trt)[!levels(test_trt) %in% light_trt],
              ' as the "heavy" treatment(s).')
      bd[[iso_trt]] <- factor(bd[[iso_trt]])
    }
    # rename factor levels
    bd[, iso_trt := factor(iso_trt, levels = levels(iso_trt), labels = c('light', 'label'))]
    # assess minimum frequency
    bd <- bd[iso_trt == 'label', rep_freq := uniqueN(sid), by = c(tax_id, grouping_cols)
             ][iso_trt == 'light', rep_freq := uniqueN(sid), by = tax_id
               ][rep_freq >= min_freq
                 ][, rep_freq := NULL]
    # remove outlier WAD values
      if(rm_outliers) {
        pos_out <- pos_outlier(bd$wad)
        neg_out <- neg_outlier(bd$wad)
      } else {
        pos_out <- Inf
        neg_out <- -Inf
      }
    bd <- bd[wad > neg_out]
    # set keys and sort for faster merging in for-loop
    bd <- setkeyv(bd, c(tax_id, grouping_cols))
    # store permutation output in a data.table
    eaf_output <- subset(bd, iso_trt == 'label', select = c(tax_id, grouping_cols))
    eaf_output <- unique(eaf_output)
    # quickly generate a list of replicate resampling permutations within your grouping variables
    reps <- subset(bd, select = c('sid', 'iso_trt', grouping_cols))
    reps <- unique(reps)
    reps[, rep_count := uniqueN(sid), by = c('iso_trt', grouping_cols)]
    resamps <- reps[, as.list(sample.int(rep_count, iters, replace = TRUE)), by = c('sid', grouping_cols)]
    #----------------
    # BEGINNING OF FOR-LOOP
    for(i in 1:iters) {
      cat('bootstrap iteration', i, 'of', iters, '\r')
      cols_for_subsample <- c('sid', grouping_cols, paste0('V', i))
      # subsample replicates
      dat_boot <- merge(bd, resamps[, ..cols_for_subsample],
                        by = c('sid', grouping_cols),
                        all.x = TRUE)
      setnames(dat_boot, paste0('V', i), 'resample_rep')
      dat_boot <- dat_boot[, wad := wad[resample_rep], by = c(tax_id, grouping_cols, 'iso_trt')]
      # convert to wide format
      dat_boot <- wad_wide(dat_boot, tax_id = tax_id, sample_id = 'sid', wads = wads, iso_trt = iso_trt, isotope = isotope)
      dat_boot <- dat_boot[isotope %in% c('18O', '13C', '15N')]
      # calculate molecular weights
      dat_boot[, gc_prop := (1 / 0.083506) * (light - 1.646057)
               ][, mw_light := (0.496 * gc_prop) + 307.691
                 ][, mw_label := (((label - light) / light) + 1) * mw_light]
      # calculate enrichment
      dat_boot[isotope == '18O', `:=` (mw_max = mw_light + 12.07747, nat_abund = nat_abund_18O)
               ][isotope == '13C', `:=` (mw_max = mw_light + 9.974564 + (-0.4987282 * gc_prop), nat_abund = nat_abund_13C)
                 ][isotope == '15N', `:=` (mw_max = mw_light + 3.517396 + (0.5024851 * gc_prop), nat_abund = nat_abund_15N)
                   ][, eaf := ((mw_label - mw_light) / (mw_max - mw_light)) * (1 - nat_abund)]
      # calculate mean enrichment for
      dat_boot <- dat_boot[, .(eaf = mean(eaf)), by = c(tax_id, grouping_cols)]
      # add iteration count to EAF columns
      setnames(dat_boot, 'eaf', paste0('eaf_', i))
      # add EAF values to output data
      cols_to_add <- c(tax_id, grouping_cols, paste0('eaf_', i))
      eaf_output <- merge(eaf_output,
                          dat_boot[, ..cols_to_add],
                          all.x = TRUE,
                          by = c(tax_id, grouping_cols),
                          sort = FALSE)
    }
    # END OF FOR-LOOP
    #----------------
    rm(dat_boot, i, cols_for_subsample, cols_to_add)
    # Calculate distribution of 18O EAF values
    eaf_output <- melt(eaf_output,
                       id.vars = c(tax_id, grouping_cols),
                       variable.name = 'iteration',
                       value.name = 'eaf',
                       na.rm = TRUE)
    # perform density correction for each tube based on distribution of EAF values
    if(correction) {
      shift <- eaf_output[!is.na(eaf)
                          ][order(eaf)
                            ][, .(shift = median(eaf[1:floor(non_grower_prop * .N)])),
                              by = c(grouping_cols, 'iteration')]
      eaf_output <- merge(eaf_output, shift, by = c(grouping_cols, 'iteration'), all.x = TRUE)
      eaf_output[, eaf := eaf - shift][, shift := NULL]
    }
    # note: your EAF confidence interval columns are called `2.5%`, `50%`, and `97.5%`
    eaf_med_ci <- eaf_output[, as.list(quantile(eaf, probs = c(0.025, .5, .975), na.rm = TRUE)),
                             by = c(tax_id, grouping_cols)]
    # calculate p-value for EAF being greater than 0
    eaf_pval <- eaf_output[, 1 - (sum(eaf > 0) / .N), by = c(tax_id, grouping_cols)]
    setnames(eaf_pval, 'V1', 'p_val')
    # combine and remove NA values -- this will also remove all the unlabeled samples
    eafd <- merge(eaf_med_ci, eaf_pval)
    eafd <- eafd[!is.na(`50%`)]
    # setnames(eafd, old = '50%', new = 'eaf_boot')
  }
  return(eafd)
}
