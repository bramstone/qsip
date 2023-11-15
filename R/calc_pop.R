#' Calculation of population dynamics
#'
#' Calculates population growth and death rates
#'
#' @param data Data as a long-format data.table where each row represents a taxonomic feature within a single fraction.
#'  Typically, this is the output from the \code{calc_wad} function.
#' @param tax_id Column name specifying unique identifier for each taxonomic feature.
#' @param sample_id Column name specifying unique identifier for each replicate.
#' @param wads Column name specifying weighted average density values.
#'  The default option is "wad" which is created by \code{calc_wad}.
#' @param iso_trt Column name specifying a two-level categorical column indicating whether a sample has been amended with a stable isotope (i.e., is "heavy") or if
#'  isotopic composition is at natural abundance (i.e., "light").
#'  Any terms may be applied but care should be taken for these values.
#'  If supplied as a factor, \code{calc_pop} will take the lowest level as the "light" treatment and the higher
#'  level as the "heavy" treatment.
#'  Alternatively, if supplied as a character, \code{calc_pop} will coerce the column to a factor and with the default behavior wherein the
#'  first value in alphabetical order will be assumed to be the lowest factor level (i.e. the "light" treatment).
#' @param time Column name specifying the timepoint at which each sample was collected.
#'  For population rates, the lowest timepoint will be assumed to be the initial timepoints with which to base
#'  changes in abundances on.
#' @param abund Column name specifying abundance measurement for each taxon in each sample.
#'  The default option is "abund" which is created by \code{calc_wad}.
#' @param mu Assumption of the fraction of oxygen atoms incorporated into new DNA that come from water (instead of from carbon sources)
#' @param bootstrap Whether to generate bootstrapped enrichment values for each taxonomic feature across groups of samples.
#'  If \code{TRUE}, replicates within specified treatment groupings will be randomly resampled with replacement and resulting
#'  WAD values used to generate a distribution of enrichment values for each taxon.
#' @param iters Integer specifying the number of bootstrap iterations to perform.
#' @param grouping_cols Column(s) used to indicate bootstrap resampling groups.
#'  Within each group, replicates will be resampled with replacement.
#'  Resampling will not occur across groups.
#' @param min_freq Minimum number of replicates a taxonomic feature must occur in to be kept.
#'  If treatment grouping columns are specified, frequencies will be assessed at this level.
#'  For unlabeled "light" samples, treatment groupings will be ignored when assessing adequate frequency of occurrence.
#' @param correction Whether to apply a correction to fractional enrichment values to ensure a certain proportion are positive.
#' @param rm_outlers Whether or not to remove fractional enrichment values that are 1.5X greater or lesser than the distance between the median
#'   and interquartile ranges.
#' @param non_grower_group Fractional value applied if \code{correction = TRUE} specifying the proportion of the community in each samples assumed to be
#'  non-growers and whose median enrichment values will be assumed to be zero. The adjustment necessary to place this median value at zero will be applied
#'  as a correction to all enrichment values in the sample.
#'
#' @details \code{calc_pop} can only calculate growth rates with samples enriched in 18O.
#'  As such, the \code{isotope} parameter is not included here, but the user must ensure
#'  that carbon and nitrogen isotopic samples are not included.
#'  The default parameter for \code{mu} is 0.6 which was determined by sensitivity analysis
#'  to be the value which produced nearly all positive growth rates in a pilot study.
#'
#'  By default, \code{calc_pop} averages the abundances of taxa from the initial time point.
#'  If \code{boostrap = TRUE}, replicates will be resampled (with replacement) then averaged so that
#'  population rates are a product of bootstrapped abundances as well as densities.
#'
#' @return \code{calc_pop} returns a data.table where each row represents a taxonomic feature within a single replicate.
#'  The following columns are produced: growth rate (\code{growth}), mortality or turnover rate (\code{mortality}).
#'
#'  Bootstrap functionality produces the columns: \code{2.5\%}, \code{50\%}, and \code{97.5\%}
#'  reflecting the median bootstrapped rate value and the 95\% confidence interval
#'  of rate values.
#'  In addition, probability that a taxon's enrichment was less than zero is represented
#'  in the \code{p_val} column which is based on the proportion sub-zero bootstrapped observations.
#'  Growth and mortality rates are differentiated by the \code{rate} column.
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
#'  # calculate population fluxes
#'  rates <- calc_pop(wads, tax_id = 'asv_id', sample_id = 'sampleID',
#'                    iso_trt = 'iso_trt', time = 'timepoint')
#'
#' @references
#'  Koch, Benjamin, \emph{et al.} 2018. Estimating taxon-specific population dynamics in diverse microbial communities.
#'  \emph{Ecosphere} \strong{9}.
#'
#' @export

calc_pop <- function(data, tax_id = c(), sample_id = c(), wads = 'wad',
                     iso_trt = c(), time = c(), abund = 'abund',
                     t0_grouping = c(), mu = 0.6,
                     bootstrap = FALSE, iters = 999L, grouping_cols = c(), min_freq = 3,
                     correction = TRUE, rm_outliers = TRUE, non_grower_prop = 0.1) {
  vars <- list(tax_id, sample_id, iso_trt, time)
  if(any(sapply(vars, is.null))) {
    null_vars <- which(sapply(vars, is.null))
    null_vars <- paste(c('taxon IDs', 'sample IDs',
                         'isotope treatment (amended or unamended)',
                         'sample timepoints')[null_vars],
                       sep = ',')
    stop("Must supply the following columns: ", null_vars)
  }
  vars <- c(vars, abund, wads)
  if(any(sapply(vars, function(x) !exists(x, data)))) {
    missing_vars <- which(sapply(vars, function(x) !exists(x, data)))
    missing_vars <- paste(vars[missing_vars], sep = ',')
    stop("Missing the following column(s) in supplied data: ", missing_vars)
  }
  # make sure timepoint is numeric
  if(!is.numeric(time)) stop('Timepoints must be supplied as numeric values')
  #
  t0 <- data[time == min(time)]
  t0_avg_cols <- c('tax_id', t0_grouping)
  #
  if(bootstrap == FALSE) {
    # average abundances across initial timepoints
    t0 <- t0[, .(abund_t0 = mean(abund)), by = t0_avg_cols]
    # convert to wide format without initial timepoints
    rated <- wad_wide(data[time > min(time)],
                      tax_id = tax_id, sample_id = sample_id, wads = wads,
                      iso_trt = iso_trt, isotope = isotope)
    #
    # correct density shifts
    if(correction) {
      label_shift <- rated[, .(wad_diff = wad_label - wad_light), by = sample_id
                           ][!is.na(wad_diff)
                             ][order(wad_diff)
                               ][, .(shift = median(wad_diff[1:floor(non_grower_prop * .N)])), by = sample_id]
      rated <- merge(rated, label_shift, by = 'sample_id', all.x = TRUE)
      rated[, wad_label := wad_label - shift][, shift := NULL]
    }
    # calculate molecular weights
    rated[, gc_prop := (1 / 0.083506) * (light - 1.646057)
          ][, mw_light := (0.496 * gc_prop) + 307.691
            ][, `:=` (mw_label = (((wad_label - wad_light) / wad_light) + 1) * mw_light,
                      mw_max = mw_light + 12.07747 * mu)]
    # merge in t0 abundances
    rated <- merge(rated, t0, by = t0_avg_cols, all.x = TRUE)
    # calculate population rates
    rated[, abund_light := abund * ((mw_max - mw_label) / (mw_max - mw_light))
          ][, `:=` (growth = log(abund / abund_light) / time,
                    mortality = log(abund_light / abund_t0) / time)]
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
    rated <- rated[growth > r_neg_out
                   ][growth < r_pos_out
                     ][mortality > m_neg_out
                       ][mortaity < m_pos_out]
    # remove NA EAF values - this will also remove all the unlabeled samples
    rated <- rated[!is.na(growth) & !is.na(mortality),
                   !c('wad_label', 'wad_light', 'gc_prop', 'mw_light', 'mw_label', 'mw_max', 'abund_light')]
  ####################
  # BOOTSTRAPPING CODE
  ####################
  } else if(bootstrap == TRUE) {
    bd <- copy(data)
    bd <- bd[time > min(time)]
    #
    setnames(bd, old = c(sample_id, iso_trt, wads), new = c('sid', 'iso_trt', 'wad'))
    setnames(t0, old = sample_id, new = 'sid')
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
    bd <- setkeyv(bd, c(tax_id, time, grouping_cols))
    t0 <- setkeyv(t0, t0_grouping)
    # store permutation output in a data.table
    rate_output <- subset(bd, iso_trt == 'label', select = c(tax_id, time, grouping_cols))
    rate_output <- unique(rate_output)
    # quickly generate a list of replicate resampling permutations within your grouping variables
    reps <- subset(bd, select = c('sid', 'iso_trt', grouping_cols))
    reps <- unique(reps)
    reps[, rep_count := uniqueN(sid), by = c('iso_trt', grouping_cols)]
    resamps <- reps[, as.list(sample.int(rep_count, iters, replace = TRUE)), by = c('sid', time, grouping_cols)]
    #
    reps_0 <- subset(t0, select = c('sid', t0_grouping))
    reps_0 <- unique(reps_0)
    reps_0[, rep_count := uniqueN(sid), by = c(t0_grouping)]
    resamps_0 <- reps[, as.list(sample.ing(rep_count, iters, replace = TRUE)), by = t0_grouping]
    #----------------
    # BEGINNING OF FOR-LOOP
    for(i in 1:iters) {
      cat('bootstrap iteration', i, 'of', iters, '\r')
      cols_for_subsample <- c('sid', time, grouping_cols, paste0('V', i))
      cols_for_subsample_0 <- c('sid', t0_grouping, paste0('V', i))
      # subsample t0 abundances
      t0_boot <- merge(t0, resamps_0[, ..cols_for_subsample_0],
                       by = c('sid', t0_grouping),
                       all.x = TRUE)
      setnames(t0_boot, paste0('V', i), 'resample_rep')
      t0_boot[, abund := abund[resample_rep], by = c(tax_id, t0_grouping)]
      # subsample replicates
      dat_boot <- merge(bd, resamps[, ..cols_for_subsample],
                        by = c('sid', time, grouping_cols),
                        all.x = TRUE)
      setnames(dat_boot, paste0('V', i), 'resample_rep')
      dat_boot <- dat_boot[, `:=` (wad = wad[resample_rep],
                                   abund = abund[resample_rep]),
                           by = c(tax_id, time, grouping_cols, 'iso_trt')]
      # convert to wide format
      dat_boot <- wad_wide(dat_boot, tax_id = tax_id, sample_id = 'sid', wads = wads, iso_trt = iso_trt, isotope = isotope)
      dat_boot <- dat_boot[iso_trt == 'label']
      # calculate molecular weights
      dat_boot[, gc_prop := (1 / 0.083506) * (light - 1.646057)
               ][, mw_light := (0.496 * gc_prop) + 307.691
                 ][, `:=` (mw_label = (((wad_label - wad_light) / wad_light) + 1) * mw_light,
                           mw_max = mw_light + 12.07747 * mu)]
      # average bootstrapped abundances across initial timepoints
      t0_boot <- t0_boot[, .(abund_t0 = mean(abund_t0)), by = t0_avg_cols]
      # merge in bootstrapped t0 abundances
      dat_boot <- merge(dat_boot, t0_boot, by = t0_avg_cols, all.x = TRUE)
      # calculate population rates
      dat_boot[, abund_light := abund * ((mw_max - mw_label) / (mw_max - mw_light))
               ][, `:=` (growth = log(abund / abund_light) / time,
                         mortality = log(abund_light / abund_t0) / time)]
      # calculate mean rates for groups
      dat_boot <- dat_boot[, .(growth = mean(growth), mortality = mean(mortality)),
                           by = c(tax_id, time, grouping_cols)]
      # add iteration count to rate columns
      setnames(dat_boot, 'growth', paste0('growth_', i))
      setnames(dat_boot, 'mortality', paste0('mortality_', i))
      # add EAF values to output data
      cols_to_add <- c(tax_id, time, grouping_cols, paste0('growth_', i), paste0('mortality_', i))
      rate_output <- merge(rate_output,
                           dat_boot[, ..cols_to_add],
                           all.x = TRUE,
                           by = c(tax_id, time, grouping_cols),
                           sort = FALSE)
    }
    # END OF FOR-LOOP
    #----------------
    rm(dat_boot, i, cols_for_subsample, cols_to_add)
    # Calculate distribution of rates
    cols <- c(tax_id, grouping_cols, sub('growth', names(rate_output)))
    growth_rate <- melt(rate_output[, ..cols],
                        id.vars = c(tax_id, time, grouping_cols),
                        variable.name = 'iteration',
                        value.name = c('growth'),
                        na.rm = TRUE)
    #
    cols <- c(tax_id, grouping_cols, sub('mortality', names(rate_output)))
    mortality_rate <- melt(rate_output[, ..cols],
                           id.vars = c(tax_id, time, grouping_cols),
                           variable.name = 'iteration',
                           value.name = c('mortality'),
                           na.rm = TRUE)
    # correct rates
    if(correction) {
      g_shift <- growth_rate[!is.na(growth)
                             ][order(growth)
                               ][, .(shift = median(growth[1:floor(non_grower_prop * .N)])),
                                 by = c(gtime, rouping_cols, 'iteration')]
      growth_rate <- merge(growth_rate, g_shift, by = c(time, grouping_cols, 'iteration'), all.x = TRUE)
      growth_rate[, growth := growth - g_shift][, g_shift := NULL]
      #
      m_shift <- mortality_rate[!is.na(mortrality)
                                ][order(mortrality)
                                  ][, .(shift = median(mortrality[1:floor(non_grower_prop * .N)])),
                                    by = c(time, grouping_cols, 'iteration')]
      mortality_rate <- merge(mortality_rate, m_shift, by = c(time, grouping_cols, 'iteration'), all.x = TRUE)
      mortality_rate[, mortality := growth + m_shift][, m_shift := NULL]
    }
    # note: your rate confidence interval columns are called `2.5%`, `50%`, and `97.5%`
    growth_med_ci <- rate_output[, as.list(quantile(growth, probs = c(0.025, .5, .975), na.rm = TRUE)),
                             by = c(tax_id, time, grouping_cols)
                             ][, rate := 'growth']
    mortality_med_ci <- rate_output[, as.list(quantile(mortality, probs = c(0.025, .5, .975), na.rm = TRUE)),
                                    by = c(tax_id, time, grouping_cols)
                                    ][, rate := 'mortality']
    # calculate p-value for rates being greater than 0
    growth_pval <- growth_rate[, .(p_val = 1 - (sum(growth > 0) / .N)), by = c(tax_id, time, grouping_cols)]
    mortality_pval <- mortality_rate[, .(p_val = 1 - (sum(mortality < 0) / .N)), by = c(tax_id, time, grouping_cols)]
    # combine and remove NA values -- this will also remove all the unlabeled samples
    growth_med_ci <- merge(growth_med_ci, growth_pval)
    mortality_med_ci <- merge(mortality_med_ci, mortality_pval)
    rated <- rbind(growth_med_ci, mortality_med_ci, use.names = T)
    rated <- rated[!is.na(`50%`)]
    # setnames(rated, old = '50%', new = 'eaf_boot')
  }
  return(rated)
}
