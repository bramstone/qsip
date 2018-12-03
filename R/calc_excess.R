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
#' @param separate_light Logical value indicating whether or not WAD-light scores should be averaged across all replicate groups or not.
#'   If \code{FALSE}, unlabeled WAD scores across all replicate groups will be averaged, creating a single molecular weight score per taxon
#'   representing it's genetic molecular weight in the absence of isotope addition.
#' @param separate_label Logical value indicating whether or not WAD-label scores should be averaged across all replicate groups or not.
#'   If \code{FALSE}, labeled WAD scores across all replicate groups will be averaged, creating a single molecular weight score per taxon
#'   representing it's genetic molecular weight as a result of isotope addition. The default is \code{TRUE}.
#'
#' @details Some details about proper isotope control-treatment factoring. If weighted average densities or the change in weighted average densities
#'   have not been calculated beforehand, \code{calc_mw} will compute those first.
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
#' @export

calc_excess <- function(data, percent=FALSE, ci_method=c('', 'bootstrap', 'bayesian'), ci=.95, iters=999, filter=FALSE,
                        correction=FALSE, offset_taxa=0.1, separate_light=FALSE, separate_label=TRUE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  ci_method <- match.arg(tolower(ci_method), c('', 'bootstrap', 'bayesian'))
  #
  # -------------------------------------------------------------
  # no CI and resampling
  #
  if(ci_method=='') {
    # if MW values don't exist, calculate those first
    # this will also handle rep_id validity (through calc_wad) and rep_group/iso_trt validity (through calc_d_wad)
    if(is.null(data@qsip[['mw_label']]) || is.null(data@qsip[['mw_light']])) {
      data <- calc_mw(data, filter=filter, correction=correction, offset_taxa=offset_taxa,
                      separate_light=separate_light, separate_label=separate_label)
    }
    # extract MW-labeled and convert to S3 matrix with taxa as ROWS (opposite all other calcs)
    mw_h <- data@qsip[['mw_label']]
    mw_h <- as(mw_h, 'matrix')
    mw_l <- data@qsip[['mw_light']]
    if(!is.null(dim(mw_l))) mw_l <- as(mw_l, 'matrix')   # if mw_l is matrix, convert to S3 matrix
    if(!phyloseq::taxa_are_rows(data)) {
      mw_h <- t(mw_h)
      if(is.matrix(mw_l)) mw_l <- t(mw_l)
    }
    tax_names <- rownames(mw_h)
    # calculate mol. weight heavy max (i.e., what is maximum possible labeling)
    if(data@qsip@iso=='18O') {
      adjust <- 12.07747
      nat_abund <- 0.002000429
    } else if(data@qsip@iso=='13C') {
      wl <- data@qsip[['wad_light']]
      wl <- as(wl, 'matrix')
      gc <- (wl - 1.646057) / 0.083506
      adjust <- (-0.4987282 * gc) + 9.974564
      nat_abund <- 0.01111233
    }
    mw_max <- adjust + mw_l
    # calculate atom excess
    if(all(dim(mw_max)==dim(mw_h))) {
      excess <- ((mw_h - mw_l)/(mw_max - mw_l)) * (1 - nat_abund)
    } else {
      num <- sweep(mw_h, 1, mw_l)
      denom <- mw_max - mw_l
      excess <- sweep(num, 1, denom, '/') * (1 - nat_abund)
    }
    # organize and add new data as S4 matrix
    if(percent) excess <- excess * 100
    data <- collate_results(data, t(excess), tax_names=tax_names, 'atom_excess', sparse=TRUE)
    return(data)
    #
    # -------------------------------------------------------------
    # CI values obtained through bootstrap subsampling
    #
  } else if(ci_method=='bootstrap') {
    # Calc WADs
    if(is.null(data@qsip[['wad']])) data <- calc_wad(data, filter=filter)
    ft <- as(data@qsip[['wad']], 'matrix')
    if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
    n_taxa <- ncol(ft)
    tax_names <- colnames(ft)
    ft <- valid_samples(data, ft, 'iso')
    iso_group <- ft[[2]]; ft <- ft[[1]]
    sam_names <- iso_group$replicate
    ft <- split_data(data, ft, iso_group$interaction, grouping_w_phylosip=FALSE)
    # how many samples in each group to subsample with?
    subsample_n <- base::lapply(ft, nrow)
    subsample <- base::lapply(subsample_n,
                              function(x) sample.int(x, size=iters*x, replace=TRUE))
    subsample <- base::mapply(matrix,
                              subsample,
                              nrow=subsample_n,
                              byrow=F, SIMPLIFY=FALSE)
    # collect output in matrix (each column is an atom excess matrix from that iterations subsampling)
    # Note: this matrix will be very large if you don't filter taxa out first
    if(length(data@qsip@rep_group)==0) {
      boot_collect <- matrix(0, nrow=n_taxa, ncol=iters)
      rownames(boot_collect) <- tax_names
      } else {
      boot_collect <- matrix(0,
                             nrow=n_taxa * nlevels(iso_group$grouping),
                             ncol=iters)
      boot_rnames <- expand.grid(tax_names,
                                 levels(iso_group$grouping),
                                 stringsAsFactors=FALSE)
      rownames(boot_collect) <- interaction(boot_rnames[,1], boot_rnames[,2], sep=':')
      rm(boot_rnames)
    }
    for(i in 1:iters) {
      # subsample WAD values, calc diff_WAD and molecular weights
      subsample_i <- lapply(subsample, function(x) x[,i])
      ft_i <- mapply(function(x, y) x[y,,drop=FALSE], ft, subsample_i, SIMPLIFY=FALSE)
      ft_i <- recombine_in_order(ft_i, iso_group, n_taxa)
      rownames(ft_i) <- sam_names
      data <- suppressWarnings(collate_results(data, ft_i, tax_names=tax_names, 'wad', sparse=TRUE))
      data <- suppressWarnings(calc_d_wad(data, correction=correction,
                                          offset_taxa=offset_taxa,
                                          separate_light=FALSE,
                                          separate_label=FALSE))
      data <- suppressWarnings(calc_mw(data))
      mw_h <- data@qsip[['mw_label']]
      mw_h <- as(mw_h, 'matrix')
      mw_l <- data@qsip[['mw_light']]
      if(!is.null(dim(mw_l))) mw_l <- as(mw_l, 'matrix')   # if mw_l is matrix, convert to S3 matrix
      if(!phyloseq::taxa_are_rows(data)) {
        mw_h <- t(mw_h)
        if(is.matrix(mw_l)) mw_l <- t(mw_l)
      }
      # calculate atom excess for this subsampling iteration
      if(data@qsip@iso=='18O') {
        adjust <- 12.07747
        nat_abund <- 0.002000429
      } else if(data@qsip@iso=='13C') {
        wl <- data@qsip[['wad_light']]
        wl <- as(wl, 'matrix')
        gc <- (wl - 1.646057) / 0.083506
        adjust <- (-0.4987282 * gc) + 9.974564
        nat_abund <- 0.01111233
      }
      mw_max <- adjust + mw_l
      # atom excess
      if(all(dim(mw_max)==dim(mw_h))) {
        excess <- ((mw_h - mw_l)/(mw_max - mw_l)) * (1 - nat_abund)
      } else {
        num <- sweep(mw_h, 1, mw_l)
        denom <- mw_max - mw_l
        excess <- sweep(num, 1, denom, '/') * (1 - nat_abund)
      }
      # organize and add data as single column in bootstrap output matrix
      boot_collect[,i] <- c(excess)
    }
    if(percent) boot_collect <- boot_collect * 100
    # clean workspace
    rm(ft_i, subsample, subsample_n, subsample_i)
    # summarize across iterations (lower CI, median, upper CI)
    ci_data <- summarize_ci(boot_collect, ci, grouping=iso_group, ncols=n_taxa, data=data)
    data <- collate_results(data, ci_data$ci_l, tax_names=tax_names, 'atom_excess_ci_l', sparse=TRUE)
    data <- collate_results(data, ci_data$med, tax_names=tax_names, 'atom_excess', sparse=TRUE)
    data <- collate_results(data, ci_data$ci_u, tax_names=tax_names, 'atom_excess_ci_u', sparse=TRUE)
    # recalculate WAD, diff_WAD, and MW values (they've been replaced by bootstrapped versions)
    data <- suppressWarnings(calc_wad(data, filter=filter))
    data <- suppressWarnings(calc_d_wad(data, correction=correction, offset_taxa=offset_taxa))
    data <- suppressWarnings(calc_mw(data))
    return(data)
    #
    # -------------------------------------------------------------
    # CI values obtained through Bayesian analysis
    #
  } else if(ci_method=='bayesian') { # method for bayesian analysis
    print('No Bayesian method yet, returning data unaltered')
    return(data)
    # code here.....
  }
}
