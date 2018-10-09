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
  if(missing(ci_method)) ci_method <- ''
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
    if(!phyloseq::taxa_are_rows(data)) mw_lab <- t(mw_lab)
    # calculate mol. weight heavy max (i.e., what is maximum possible labeling)
    if(data@qsip@iso=='18O') {
      adjust <- 12.07747 + mw_l
      nat_abund <- 0.002000429
    }
    else if(data@qsip@iso=='13C') {
      adjust <- (-0.4987282 * gc) + 9.974564
      nat_abund <- 0.01111233
    }
    mw_max <- adjust + mw_l
    # calculate atom excess
    excess <- ((mw_lab - mw_l)/(mw_max - mw_l)) * (1 - nat_abund)
    # organize and add new data as S4 matrix
    if(percent) excess <- excess * 100
    data <- collate_results(data, t(excess), 'atom_excess', sparse=TRUE)
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
    iso_group <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
    # keep only valid rows
    keep_rows <- (iso_group$replicate %in% rownames(ft) & !is.na(iso_group$iso))
    if(sum(!keep_rows) > 0) {
      warning('Dropping group(s): ',
              paste(as.character(iso_group$replicate[!keep_rows]), collapse=', '),
              ' - from calculation', call.=FALSE)
    }
    iso_group <- iso_group[iso_group$replicate %in% rownames(ft),]
    ft <- ft[!is.na(iso_group$iso),]
    iso_group <- iso_group[!is.na(iso_group$iso),]
    iso_group <- iso_group[match(rownames(ft), iso_group$replicate),] # match row orders to ft
    # split by replicate groups
    ft <- split_data(data, ft, iso_group$interaction, grouping_w_phylosip=FALSE)
    # how many samples in each group to subsample with?
    subsample_n <- base::lapply(ft, nrow)
    subsample <- base::lapply(subsample_n,
                              function(x) sample.int(x, size=iters*x, replace=TRUE))
    subsample <- base::mapply(matrix,
                              subsample,
                              nrow=subsample_n,
                              byrow=FALSE)
    # collect output in matrix (each column is an atom excess matrix from that iterations subsampling)
    # Note: this matrix will be very large if you don't filter taxa out first
    boot_collect <- matrix(0,
                           nrow=ncol(ft) * nrow(ft),
                           ncol=iters)
    boot_rnames <- expand.grid(rownames(excess_i),
                               colnames(excess_i),
                               stringsAsFactors=FALSE)
    rownames(boot_collect) <- interaction(boot_rnames[,1], boot_rnames[,2], sep=':')
    rm(boot_rnames)
    for(i in 1:iters) {
      # subsample WAD values, calc diff_WAD and molecular weights
      subsample_i <- lapply(subsample, function(x) x[,i])
      ft_i <- mapply(function(x, y) x[y,], ft, subsample_i)
      ft_i <- do.call(rbind, ft_i)
      data <- collate_results(data, ft_i, 'wad', sparse=TRUE)
      data <- calc_mw(data)
      mw_lab <- data@qsip[['mw_label']]
      mw_lab <- as(mw_lab, 'matrix')
      mw_l <- data@qsip[['mw_light']]
      if(!is.null(dim(mw_l))) mw_l <- as(mw_l, 'matrix')   # if mw_l is matrix, convert to S3 matrix
      if(!phyloseq::taxa_are_rows(data)) mw_lab <- t(mw_lab)
      # calculate atom excess for this subsampling iteration
      if(data@qsip@iso=='18O') {
        adjust <- 12.07747 + mw_l
        nat_abund <- 0.002000429
      }
      else if(data@qsip@iso=='13C') {
        adjust <- (-0.4987282 * gc) + 9.974564
        nat_abund <- 0.01111233
      }
      mw_max <- adjust + mw_l
      # atom excess
      excess <- ((mw_lab - mw_l)/(mw_max - mw_l)) * (1 - nat_abund)
      # organize and add data as single column in bootstrap output matrix
      if(percent) excess <- excess * 100
      colnames(excess) <- colnames(ft)
      boot_collect[,i] <- c(excess)
    }
    rm(ft_i, subsample, subsample_n)
    # summarize across iterations (lower CI, median, upper CI)
    summaries <- t(apply(boot_collect, 1,
                         quantile,
                         c((1 - ci)/2, .5, (1 - ci)/2 + ci)))
    ci_l <- summaries[,1]
    med <- summaries[,2]
    ci_u <- summaries[,3]
    rm(boot_collect, summaries)
    # reconstruct matrices with taxa as columns and output
    ci_l <- matrix(ci_l, ncol=phyloseq::ntaxa(data), byrow=TRUE)
    med <- matrix(med, ncol=phyloseq::ntaxa(data), byrow=TRUE)
    ci_u <- matrix(ci_u, ncol=phyloseq::ntaxa(data), byrow=TRUE)
    rownames(ci_l) <- rownames(med) <- rownames(ci_u) <- rownames(excess)
    ci_l_name <- paste0('ae_',ci*100,'_ci_l')
    ci_u_name <- paste0('ae_',ci*100,'_ci_u')
    data <- collate_results(data, ci_l, ci_l_name, sparse=TRUE)
    data <- collate_results(data, med, 'atom_excess', sparse=TRUE)
    data <- collate_results(data, ci_l, ci_u_name, sparse=TRUE)
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
