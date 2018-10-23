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
#'
#' @details Some details about proper isotope control-treatment factoring. If weighted average densities or the change in weighted average densities
#'   have not been calculated beforehand, \code{calc_mw} will compute those first.
#'
#' @return \code{calc_excess} adds an S4 Matrix class objects (which more efficiently stores sparse matrix data) to the \code{data@@qsip@@.Data} slot
#'   of molecular weights for each taxon at each group of replicates in the labeled and unlabeled groups. The row and column
#'   specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}}, meaning if taxa are listed on the table rows,
#'   they will in the resulting S4 Matrix class.
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

calc_excess <- function(data, percent=FALSE, ci_method=c('', 'bootstrap', 'bayesian'), ci=.95, iters=999, filter=FALSE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  ci_method <- match.arg(tolower(ci_method), c('', 'bootstrap', 'bayesian'))
  #
  # -------------------------------------------------------------
  # no CI and resampling
  #
  if(ci_method=='') {
    # if MW values don't exist, calculate those first
    # this will also handle rep_id validity (through calc_wad) and rep_group/iso_trt validity (through calc_d_wad)
    if(is.null(data@qsip[['mw_label']]) || is.null(data@qsip[['mw_light']])) data <- calc_mw(data, filter=filter)
    # extract MW-labeled and convert to S3 matrix with taxa as ROWS (opposite all other calcs)
    mw_lab <- data@qsip[['mw_label']]
    mw_lab <- as(mw_lab, 'matrix')
    mw_l <- data@qsip[['mw_light']]
    if(!is.null(dim(mw_l))) mw_l <- as(mw_l, 'matrix')   # if mw_l is matrix, convert to S3 matrix
    if(!phyloseq::taxa_are_rows(data)) {
      mw_lab <- t(mw_lab)
      if(is.matrix(mw_l)) mw_l <- t(mw_l)
    }
    tax_names <- rownames(mw_lab)
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
    excess <- ((mw_lab - mw_l)/(mw_max - mw_l)) * (1 - nat_abund)
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
    iso_group <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
    ft <- ft[match(iso_group$replicate, rownames(ft)),] # match row orders to replicate IDs
    # keep only valid rows
    keep_rows <- (iso_group$replicate %in% rownames(ft) & !is.na(iso_group$iso))
    if(sum(!keep_rows) > 0 && is.null(data@qsip[['d_wad']])) {
      warning('Dropping group(s): ',
              paste(as.character(iso_group$replicate[!keep_rows]), collapse=', '),
              ' - from calculation', call.=FALSE)
    }
    iso_group <- iso_group[iso_group$replicate %in% rownames(ft),]
    ft <- ft[!is.na(iso_group$iso),]
    iso_group <- iso_group[!is.na(iso_group$iso),]
    # split by replicate groups
    sam_names <- rownames(ft)
    iso_group$interaction <- factor(iso_group$interaction) # limit to existing combinations only
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
    if(isTRUE(all.equal(iso_group$iso, iso_group$grouping))) {
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
      ft_i <- mapply(function(x, y) x[y,], ft, subsample_i, SIMPLIFY=FALSE)
      ft_i <- do.call(rbind, ft_i)
      rownames(ft_i) <- sam_names
      data <- suppressWarnings(collate_results(data, ft_i, tax_names=tax_names, 'wad', sparse=TRUE))
      data <- suppressWarnings(calc_d_wad(data))
      data <- suppressWarnings(calc_mw(data))
      mw_lab <- data@qsip[['mw_label']]
      mw_lab <- as(mw_lab, 'matrix')
      mw_l <- data@qsip[['mw_light']]
      if(!is.null(dim(mw_l))) mw_l <- as(mw_l, 'matrix')   # if mw_l is matrix, convert to S3 matrix
      if(!phyloseq::taxa_are_rows(data)) {
        mw_lab <- t(mw_lab)
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
      excess <- ((mw_lab - mw_l)/(mw_max - mw_l)) * (1 - nat_abund)
      # organize and add data as single column in bootstrap output matrix
      if(percent) excess <- excess * 100
      boot_collect[,i] <- c(excess)
    }
    # clean workspace
    rm(ft_i, subsample, subsample_n, subsample_i)
    # summarize across iterations (lower CI, median, upper CI)
    ci_data <- summarize_ci(boot_collect, ci, grouping=iso_group, ncols=n_taxa)
    data <- collate_results(data, ci_data$ci_l, tax_names=tax_names, 'atom_excess_ci_l', sparse=TRUE)
    data <- collate_results(data, ci_data$med, tax_names=tax_names, 'atom_excess', sparse=TRUE)
    data <- collate_results(data, ci_data$ci_u, tax_names=tax_names, 'atom_excess_ci_u', sparse=TRUE)
    # recalculate WAD, diff_WAD, and MW values (they've been replaced by bootstrapped versions)
    data <- suppressWarnings(calc_wad(data, filter=filter))
    data <- suppressWarnings(calc_d_wad(data))
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
