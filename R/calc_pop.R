#' Calculation of population dynamics
#'
#' Calculates population growth and death rates
#'
#' @param data Data as a \code{phyloseq} object
#' @param ci_method Character value indicating how to calculate confidence intervals of stable isotope atom excess.
#'   Options are \code{bootstrap} or \code{bayesian} (see \code{details} below for discussion on their differences).
#'   The default is blank indicating that no confidence intervals will be calculated.
#' @param ci Numeric value from 0 to 1 indicating the width of the confidence interval for bootsrapped atom excess values.
#' @param iters Number of (subsampling) iterations to conduct to calculate confidence intervals. Default is \code{999}.
#' @param filter Logical vector specifying whether or not to filter taxa from the weighted average density calculation.
#'   This will require \code{data} to have a filter applied with \code{\link{filter_qsip}}.
#' @param growth_model Character vector specifying whether growth rates should be calculated using exponential or linear growth.
#'   Default is \code{exponential}.
#' @param mu Numeric value from 0 to 1 indicating the expected incorporation of labeled O into DNA during replication.
#'   Default is \code{0.60}.
#' @param correction Logical value indicating whether or not to apply tube-level correction to labeled WAD values.
#' @param offset_taxa Value from 0 to 1 indicating the percentage of the taxa to utilize for calculating offset correction values.
#'   Taxa are ordered by lowest difference in WAD values.
#'   Default is \code{0.1} indicating 10 percent of taxa with the lowest difference in WAD values.
#'
#' @details Some details about proper isotope control-treatment factoring and timepoint specification. If weighted average densities or the
#'   change in weighted average densities have not been calculated beforehand, \code{calc_pop} will compute those first.
#'
#' @return \code{calc_pop} adds two S4 Matrix class objects (which more efficiently stores sparse matrix data) to the \code{data@@qsip@@.Data} slot
#'   of population birth rates for each taxon at each group of replicates. The row and column specifications will mirror those of the \code{phylosip}'s
#'   \code{\link{otu_table}}, meaning if taxa are listed on the table rows, they will in the resulting S4 Matrix class.
#'
#'   Timepoint should be in units of days, so that birth will be new 16S copies d-1 and death will be loss of light 16S copies d-1.
#'   Use of different time increments will yield growth rates (e.g. per hour), but must be appropriate for the frequency of sampling.
#'
#' @seealso \code{\link{calc_mw}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate population fluxes
#'
#' @export

calc_pop <- function(data, ci_method=c('', 'bootstrap', 'bayesian'), ci=.95, iters=999, filter=FALSE, growth_model=c('exponential', 'linear'), mu=0.6, correction=FALSE, offset_taxa=0.1) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  ci_method <- match.arg(tolower(ci_method), c('', 'bootstrap', 'bayesian'))
  growth_model <- match.arg(tolower(growth_model), c('exponential', 'linear'))
  if(data@qsip@iso!='18O') stop('Must use 18O-labeled treatment to calculate population change')
  if(length(data@qsip@timepoint)==0) stop('Must specify different sample times with timepoint')
  times <- data@sam_data[[data@qsip@timepoint]]
  if(nlevels(times)==1 || length(unique(times))==1) stop('Only one timepoint present in the data - cannot calculate population change')
  rm(times)
  #
  # -------------------------------------------------------------
  # no CI and resampling
  #
  if(ci_method=='') {
    # if MW values don't exist, calculate those first
    # this will also handle rep_id validity (through calc_wad) and rep_group/iso_trt validity (through calc_d_wad)
    if(is.null(data@qsip[['mw_label']]) || is.null(data@qsip[['mw_light']])) data <- calc_mw(data, filter=filter,
                                                                                             correction=correction, offset_taxa=offset_taxa)
    # transform sequencing abundances to 16S copy numbers
    # returns feature table (as matrix) with taxa as columns, samples as rows
    ft <- copy_no(data)
    n_taxa <- ncol(ft)
    if(filter) {
      tax_names <- data@qsip@filter
      ft <- ft[,colnames(ft) %in% tax_names]
    }
    tax_names <- colnames(ft)
    # calculate per-taxon total 16S copy abundance for each sample
    ft <- split_data(data, ft, data@qsip@rep_id)
    ft <- lapply(ft, colSums, na.rm=T)
    ft <- do.call(rbind, ft)
    # separate samples based on timepoint, keeping only valid samples
    ft <- valid_samples(data, ft, 'time')
    time_group <- ft[[2]]; ft <- ft[[1]]
    sam_names <- rownames(ft)
    # calculate per-taxon average 16S copy abundance for each group:time interaction point
    # it should be average in case time 0 data weren't fractioned
    ft <- split_data(data, ft, time_group$interaction, grouping_w_phylosip=F)
    ft <- lapply(ft, colMeans, na.rm=T)
    ft <- do.call(cbind, ft)
    ft[ft==0] <- NA
    # get 16S copy numbers for different timepoints
    time_group2 <- unique(time_group[,!names(time_group) %in% 'replicate']) # only get unique elements to match levels in ft
    ft <- ft[,match(time_group2$interaction, colnames(ft))] # re-order columns to match time_group2$interaction
    n_t_names <- paste0('n_t_',levels(time_group2$time))
    if(!any(time_group2$time==0)) {
      warning('No timepoints designated as time 0; using time ',
              levels(time_group2$time)[1],
              ' as time before isotope addition', call.=FALSE)
    }
    # NOTE: This function uses t in for-loop iterations to designate cycle through different timepoints
    for(t in 1:nlevels(time_group2$time)) {
      assign(n_t_names[t],
             ft[,time_group2$time==levels(time_group2$time)[t]])
    }
    # extract MW-labeled and convert to S3 matrix with taxa as ROWS (opposite all other calcs)
    mw_lab <- data@qsip[['mw_label']]
    mw_lab <- as(mw_lab, 'matrix')
    mw_l <- data@qsip[['mw_light']]
    if(!is.null(dim(mw_l))) mw_l <- as(mw_l, 'matrix')   # if mw_l is matrix, convert to S3 matrix, then to vector
    if(!phyloseq::taxa_are_rows(data)) mw_lab <- t(mw_lab)
    if(!phyloseq::taxa_are_rows(data) && is.matrix(mw_l)) mw_l <- t(mw_l)
    # calculate proportion in light fraction (N_light) at any time after 0
    n_l_names <- paste0('n_l_', levels(time_group2$time))
    for(t in 2:nlevels(time_group2$time)) {
      if(length(data@qsip@rep_group)==0) {
        # if there is no grouping separate from timepoints, go by the column timepoints
        mw_lab_t <- mw_lab[,t-1]
        mw_l_t <- as.matrix(mw_l)[,t-1]
        # else break mw_label and mw_light into list separated by timepoint
      } else {
        time_group_t <- time_group2[as.numeric(time_group2$time) > 1,]
        time_group_t$time <- factor(time_group_t$time)
        mw_lab_t <- split_data(data, t(mw_lab), time_group_t$time, grouping_w_phylosip=FALSE)
        mw_l_t <- suppressWarnings(split_data(data, t(mw_l), time_group_t$time, grouping_w_phylosip=FALSE))
        mw_lab_t <- t(mw_lab_t[[t - 1]])
        mw_l_t <- t(mw_l_t[[t - 1]])
      }
      # calculate mol. weight heavy max (i.e., what is maximum possible labeling)
      mw_max <- (12.07747 * mu) + mw_l_t
      if(all(dim(mw_max)==dim(mw_lab_t))) {
        n <- ((mw_max - mw_lab_t)/(mw_max - mw_l_t)) * get(n_t_names[t])
      } else {
        num <- sweep(mw_lab_t, 1, mw_max) * -1
        denom <- sweep(mw_l_t, 1, mw_max) * -1
        n <- sweep(num, 1, denom, '/') * get(n_t_names[t])
      }
      colnames(n) <- colnames(get(n_t_names[t]))
      # remove abundances less than 0 (occurs when labeled MWs are heavier than heavymax)
      n[n < 0] <- NA
      assign(n_l_names[t], n)
    }; suppressWarnings(rm(n, mw_lab_t, mw_l_t, time_group_t, mw_max))
    # calculate birth and death rate for each timepoint after 0
    b_names <- paste0('b_', levels(time_group2$time))
    d_names <- paste0('d_', levels(time_group2$time))
    for(t in 2:nlevels(time_group2$time)) {
      incubate_time <- as.numeric(levels(time_group2$time)[t])
      if(growth_model=='exponential') {
        b <- get(n_t_names[t]) / get(n_l_names[t])
        d <- get(n_l_names[t]) / get(n_t_names[1])
        b <- log(b) / incubate_time
        d <- log(d) / incubate_time
      } else {
        b <- get(n_t_names[t]) - get(n_l_names[t])
        d <- get(n_l_names[t]) - get(n_t_names[1])
        b <- (b / get(n_t_names[1])) / incubate_time
        d <- (d / get(n_t_names[1])) / incubate_time
      }
      colnames(b) <- colnames(d) <- colnames(get(n_t_names[t]))
      assign(b_names[t], b)
      assign(d_names[t], d)
    }; rm(b,d)
    # if more than two timepoints (0, and t), combine resulting matrices, but skip time 0
    if(length(n_t_names) > 2) {
      b <- t(do.call(cbind, mget(b_names[2:length(b_names)])))
      d <- t(do.call(cbind, mget(d_names[2:length(d_names)])))
    } else {
      b <- get(b_names[2])
      d <- get(d_names[2])
    }
    # organize and add new data as S4 matrices
    data <- collate_results(data, b, tax_names=tax_names, 'birth_rate', sparse=T)
    data <- collate_results(data, d, tax_names=tax_names, 'death_rate', sparse=T)
    data <- collate_results(data, b + d, tax_names=tax_names, 'growth_rate', sparse=T)
    return(data)
  #
  # -------------------------------------------------------------
  # CI values obtained through bootstrap subsampling (will need to bootstrap abundances AND WADs/MWs)
  #
  } else if(ci_method=='bootstrap') {
    # calculate and create subsampling criteria for WADs......................................
    if(is.null(data@qsip[['wad']])) data <- calc_wad(data, filter=filter)
    wads <- as(data@qsip[['wad']], 'matrix')
    if(phyloseq::taxa_are_rows(data)) wads <- t(wads)
    n_taxa <- ncol(wads)
    tax_names <- colnames(wads)
    # keep only valid samples
    wads <- valid_samples(data, wads, 'iso', quiet=TRUE)
    iso_group <- wads[[2]]; wads <- wads[[1]]
    # split by replicate groups
    sam_names_wads <- rownames(wads)
    wads <- split_data(data, wads, iso_group$interaction, grouping_w_phylosip=FALSE)
    # how many samples in each group to subsample WADS with?
    subsample_n <- base::lapply(wads, nrow)
    subsample_wads <- base::lapply(subsample_n,
                              function(x) sample.int(x, size=iters*x, replace=TRUE))
    subsample_wads <- base::mapply(matrix,
                                   subsample_wads,
                                   nrow=subsample_n,
                                   byrow=F, SIMPLIFY=FALSE)
    # calculate and create subsampling criteria for 16S gene copy number......................
    # transform sequencing abundances to 16S copy numbers (taxa as columns)
    ft <- copy_no(data)
    if(filter) ft <- ft[,colnames(ft) %in% tax_names]
    # calculate per-taxon 16S copy abundance for each sample
    ft <- split_data(data, ft, data@qsip@rep_id)
    ft <- lapply(ft, colSums, na.rm=T)
    ft <- do.call(rbind, ft)
    # separate abundances based on timepoint, keeping only valid samples
    ft <- valid_samples(data, ft, 'time')
    time_group <- ft[[2]]; ft <- ft[[1]]
    sam_names <- time_group$replicate
    ft <- split_data(data, ft, time_group$interaction, grouping_w_phylosip=F)
    # how many samples in each group to subsample with?
    subsample_n <- base::lapply(ft, nrow)
    subsample <- base::lapply(subsample_n,
                              function(x) sample.int(x, size=iters*x, replace=TRUE))
    subsample <- base::mapply(matrix,
                              subsample,
                              nrow=subsample_n,
                              byrow=F, SIMPLIFY=FALSE)
    # collect output in matrices (each column is a pop matrix from that iterations' subsampling)
    if(length(data@qsip@rep_group)==0) {
      n_groups <- 1
      boot_rnames <-  expand.grid(tax_names,
                                  levels(time_group$time)[2:nlevels(time_group$time)],
                                  stringsAsFactors=FALSE)
      boot_rnames <- interaction(boot_rnames[,1], boot_rnames[,2], sep=':')
    } else {
      n_groups <- nlevels(time_group$grouping)
      boot_rnames <- expand.grid(tax_names,
                                 levels(time_group$grouping),
                                 levels(time_group$time)[2:nlevels(time_group$time)],
                                 stringsAsFactors=FALSE)
      boot_rnames <- interaction(boot_rnames[,1], boot_rnames[,2], boot_rnames[,3], sep=':')
    }
    n_timepoints <- nlevels(time_group$time) - 1
    boot_collect_b <- matrix(0, nrow=n_taxa * n_groups * n_timepoints, ncol=iters)
    rownames(boot_collect_b) <- boot_rnames
    boot_collect_d <- boot_collect_b
    rm(boot_rnames)
    #
    for(i in 1:iters) {
      # subsample WADs
      subsample_i_wads <- lapply(subsample_wads, function(x) x[,i])
      wads_i <- mapply(function(x, y) x[y,,drop=FALSE], wads, subsample_i_wads, SIMPLIFY=FALSE)
      wads_i <- recombine_in_order(wads_i, iso_group, n_taxa)
      rownames(wads_i) <- sam_names_wads
      # calc diff_WADs, MWs, and N values
      data <- suppressWarnings(collate_results(data, wads_i, tax_names=tax_names, 'wad', sparse=TRUE))
      data <- suppressWarnings(calc_d_wad(data, correction=correction, offset_taxa=offset_taxa, separate_wad_light=TRUE))
      data <- suppressWarnings(calc_mw(data))
      mw_lab <- data@qsip[['mw_label']]
      mw_lab <- as(mw_lab, 'matrix')
      mw_l <- data@qsip[['mw_light']]
      if(!is.null(dim(mw_l))) mw_l <- as(mw_l, 'matrix')
      if(!phyloseq::taxa_are_rows(data)) mw_lab <- t(mw_lab)
      if(!phyloseq::taxa_are_rows(data) && is.matrix(mw_l)) mw_l <- t(mw_l)
      # subsample abundances
      # calculate per-taxon average 16S copy abundance for each group:time interaction point
      subsample_i <- lapply(subsample, function(x) x[,i])
      ft_i <- mapply(function(x, y) x[y,,drop=FALSE], ft, subsample_i, SIMPLIFY=FALSE)
      ft_i <- lapply(ft_i, colMeans, na.rm=T)
      ft_i <- t(recombine_in_order(ft_i, time_group, n_taxa, condensed_grouping=TRUE))
      colnames(ft_i) <- time_group$interaction[!duplicated(time_group$interaction)]
      ft_i[ft_i==0] <- NA
      # get per-taxon 16S copy numbers for different timepoints
      time_group2 <- unique(time_group[,!names(time_group) %in% 'replicate']) # only get unique elements to match levels in ft
      ft_i <- ft_i[,match(time_group2$interaction, colnames(ft_i))] # re-order columns to match time_group2$interaction
      n_t_names <- paste0('n_t_',levels(time_group2$time))
      # t represents different timepoints
      for(t in 1:nlevels(time_group2$time)) {
        assign(n_t_names[t],
               ft_i[,time_group2$time==levels(time_group2$time)[t]])
      }
      # calculate proportion in light fraction (N_light) at any time after 0
      n_l_names <- paste0('n_l_', levels(time_group2$time))
      for(t in 2:nlevels(time_group2$time)) {
        if(length(data@qsip@rep_group)==0) {
          # if there is no grouping separate from timepoints, go by the column timepoints
          mw_lab_t <- mw_lab[,t-1]
          mw_l_t <- as.matrix(mw_l)[,t-1]
          # else break mw_label and mw_light into list separated by timepoint
        } else {
          time_group_t <- time_group2[as.numeric(time_group2$time) > 1,]
          time_group_t$time <- factor(time_group_t$time)
          mw_lab_t <- split_data(data, t(mw_lab), time_group_t$time, grouping_w_phylosip=FALSE)
          mw_l_t <- suppressWarnings(split_data(data, t(mw_l), time_group_t$time, grouping_w_phylosip=FALSE))
          mw_lab_t <- t(mw_lab_t[[t - 1]])
          mw_l_t <- t(mw_l_t[[t - 1]])
        }
        # calculate mol. weight heavy max (i.e., what is maximum possible labeling)
        mw_max <- (12.07747 * mu) + mw_l_t
        if(all(dim(mw_max)==dim(mw_lab_t))) {
          n <- ((mw_max - mw_lab_t)/(mw_max - mw_l_t)) * get(n_t_names[t])
        } else {
          num <- sweep(mw_lab_t, 1, mw_max) * -1
          denom <- sweep(mw_l_t, 1, mw_max) * -1
          n <- sweep(num, 1, denom, '/') * get(n_t_names[t])
        }
        colnames(n) <- colnames(get(n_t_names[t]))
        # remove abundances less than 0 (occurs when labeled MWs are heavier than heavymax)
        n[n < 0] <- NA
        assign(n_l_names[t], n)
      }; suppressWarnings(rm(n, mw_lab_t, mw_l_t, time_group_t, mw_max))
      # calculate birth and death rate for each timepoint after 0
      b_names <- paste0('b_', levels(time_group2$time))
      d_names <- paste0('d_', levels(time_group2$time))
      for(t in 2:nlevels(time_group2$time)) {
        incubate_time <- as.numeric(levels(time_group2$time)[t])
        if(growth_model=='exponential') {
          b <- get(n_t_names[t]) / get(n_l_names[t])
          d <- get(n_l_names[t]) / get(n_t_names[1])
          b <- log(b) / incubate_time
          d <- log(d) / incubate_time
        } else {
          b <- get(n_t_names[t]) - get(n_l_names[t])
          d <- get(n_l_names[t]) - get(n_t_names[1])
          b <- (b / get(n_t_names[1])) / incubate_time
          d <- (d / get(n_t_names[1])) / incubate_time
        }
        assign(b_names[t], b)
        assign(d_names[t], d)
      }; rm(b,d)
      # if more than two timepoints (0, and t), combine resulting matrices, but skip time 0
      if(length(n_t_names) > 2) {
        b <- do.call(cbind, mget(b_names[2:length(b_names)]))
        d <- do.call(cbind, mget(d_names[2:length(d_names)]))
      } else {
        b <- get(b_names[2])
        d <- get(d_names[2])
      }
      # organize and add data as single columns in bootstrap output matrices
      boot_collect_b[,i] <- c(b)
      boot_collect_d[,i] <- c(d)
    }
    # END OF BOOTSTRAP ITERATIONS
    #
    # clean workspace
    rm(ft_i, wads_i, subsample, subsample_wads,
       subsample_n, subsample_i, subsample_i_wads, b, d,
       b_names, d_names, n_l_names, n_t_names, ft, wads, mw_lab, mw_l)
    # summarize birth, death, flux across iterations (lower CI, median, upper CI)
    ci_birth <- summarize_ci(boot_collect_b, ci,
                             grouping=time_group,
                             ncols=n_taxa,
                             data=data)
    ci_death <- summarize_ci(boot_collect_d, ci,
                             grouping=time_group,
                             ncols=n_taxa,
                             data=data)
    ci_net <- summarize_ci(boot_collect_b + boot_collect_d, ci,
                            grouping=time_group,
                            ncols=n_taxa,
                            data=data)
    rm(boot_collect_b, boot_collect_d)
    # collate results
    objects <- c('ci_birth', 'ci_death', 'ci_net')
    metric <- c('birth_rate', 'death_rate', 'growth_rate')
    ci_level <- c('ci_l', '', 'ci_u')
    for(i in 1:3) {
      for(j in 1:3) {
        data <- collate_results(data,
                                get(objects[i])[[j]],
                                tax_names=tax_names,
                                metric=sub('_{1}$', '', paste(metric[i], ci_level[j], sep='_')),
                                sparse=TRUE)
      }
    }
    # recalculate WAD, diff_WAD, and MW values (they've been replaced by bootstrapped versions)
    data <- suppressWarnings(calc_wad(data, filter=filter))
    data <- suppressWarnings(calc_d_wad(data, correction=correction, offset_taxa=offset_taxa))
    data <- suppressWarnings(calc_mw(data))
    return(data)
  #
  # -------------------------------------------------------------
  # CI values obtained through bootstrap subsampling
  #
  } else if(ci_method=='bayesian') {
    print('No Bayesian method yet, returning data unaltered')
    return(data)
    # code here.....
  }
  return(data)
}
