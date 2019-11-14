#' Calculation of weighted average densities
#'
#' Calculates weighted average densities for each microbial taxa in each sample replicate
#'
#' @param data Data as a \code{phylosip} object
#' @param filter Logical vector specifying whether or not to filter taxa from the weighted average density calculation.
#'   A hard filter is applied for \code{calc_wad} (meaning taxa are removed from all replicates/groups if they don't meet threshold).
#'   Note, however, that this will have different outcomes depending on whether or not \code{filter_qsip} has been called
#'   before calculations or not (see \code{details})
#' @param pool_unlabeled Logical vector specifying if unlabeled replicates should be pooled together across any grouping factor prior
#'   to filtering. If \code{TRUE} (the default), unlabeled replicates will be pooled, and any soft filtering threshold will be applied to
#'   \emph{all} unlabeled replicates together.
#' @param calc_wvd Logical vector specifying whether or not to calculate weighted variance of densities (WVD).
#'   The weighted variance of densities of a taxon within a replicate is the complement to its weighted average.
#'   It describes how much variance exists around every given weighted average density.
#'   If \code{TRUE}, another output matrix (\code{wvd}) will be returned with identical dimensions to the weighted average
#'   density output.
#'
#' @details Specifying \code{na.rm=TRUE} will allow \code{calc_wad} to calculate weighted average density values from samples
#'   that have one or more fractions without a valid density value. The default setting, \code{na.rm=FALSE}, returns values
#'   of \code{NA} for every taxa in a sample with missing density data.
#'
#'   Filtering for \code{calc_wad} utilized \code{filter_qsip} and employs a hard filter. If filtering is specified using
#'   \code{filter=TRUE}, then filtering is employed on fractions \emph{only} (\emph{i.e.}, \code{filter_qsip} is implemented with
#'   \code{replicate=1}). In order to set more stringent hard filters, \code{filter_qsip} must be employed before calculating WAD values.
#'   We \emph{strongly} recommend that \code{rm_low_freq} be set as \code{TRUE}, as this should reduce the occurrence of
#'   negative shifts in density.
#'
#'   The weighted average density for taxon \emph{i} in replicate \emph{j}, designated as \eqn{W_{ij}}, is:
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
#' @return \code{calc_wad} adds an S4 Matrix class (which more efficiently stores sparse matrix data) to the \code{.Data} slot within
#'   the \code{qsip} slot of the object. of weighted average density values for each taxon at each sample. The row and column
#'   specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}}, meaning if taxa are listed on the table rows,
#'   they will in the resulting S4 Matrix class.
#'
#' @seealso \code{\link{filter_qsip}}
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

calc_wad <- function(data, filter=FALSE, pool_unlabeled=TRUE, calc_wvd=FALSE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  if(length(data@qsip@rep_id)==0) stop('Must specify replicate IDs with rep_id')
  # transform sequencing abundances to 16S copy numbers
  # returns matrix with taxa as columns, samples as rows
  ft <- pa <- copy_no(data)
  pa <- ceiling(pa/max(pa, na.rm=T))
  storage.mode(pa) <- 'integer'
  n_taxa <- ncol(ft)
  # manipulate data matrix and calculate
  ft <- split_data(data, ft, data@qsip@rep_id) # split by replicate IDs
  dv <- split(data@sam_data[[data@qsip@density]],
              data@sam_data[[data@qsip@rep_id]]) # split densities by replicate IDs
  # identify NA densities and make the corresponding abundance calcs NA for those fractions
  dv_nas <- lapply(dv, function(x) which(is.na(x)))
  ft <- base::Map(function(x, y) {x[y,] <- NA; x}, ft, dv_nas)
  # calculate relative abundances and WADs
  ft <- base::lapply(ft, function(x) {x <- t(x); x <- t(x / rowSums(x, na.rm=T)); x[is.nan(x)] <- 0; x}) # create relative abundances
  ft <- base::Map(function(y, x) sweep(y, 1, x, '*'), ft, dv)
  ft <- base::lapply(ft, colSums, na.rm=T)
  # if specified, calculate weighted variance of densities
  if(calc_wvd) {
    wvd <- base::Map(function(x, y) outer(x, y, '-')^2, dv, ft)
    wvd <- base::Map(function(y, x) y * x, wvd, ft)
    wvd <- base::lapply(wvd, colSums, na.rm=T)
  }
  # apply filtering first if desired.
  # 1. Soft filter
  # Here, taxa who do not meet the threshold(s) have their group-specific WAD values converted to 0
  if(filter && any(data@qsip@filter_levels$soft)) {
    filter_levels <- data@qsip@filter_levels
    fraction <- filter_levels$frac_freq[which(filter_levels$soft)[1]]
    replicate <- filter_levels$rep_freq[which(filter_levels$soft)[1]]
    pa <- split_data(data, pa, data@qsip@rep_id)
    pa <- base::lapply(pa, colSums, na.rm=T)
    pa <- do.call(rbind, pa)
    sam_names <- rownames(pa)
    # fraction filtering
    pa <- ifelse(pa >= fraction, 1, 0)
    # Apply Alicia's filtering criteria (between-replicate fraction filtering):
    # remove taxa WAD/WVD values only in the replicates where they fail to meet fraction thresholds
    rm_low_freq <- filter_levels$rm_low_freq[which(filter_levels$soft)[1]]
    if(rm_low_freq) {
      pa <- split(pa, sam_names)
      pa <- pa[names(ft)]
      ft <- base::Map(function(x,y) {x[y==0] <- 0; x}, ft, pa)
      if(calc_wvd) {
        wvd <- base::Map(function(x,y) {x[y==0] <- 0; x}, wvd, pa)
      }
    }
    # replicate-treatment grouping filtering
    iso_group <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
    iso_group <- iso_group[match(rownames(pa), iso_group$replicate),]
    pa <- pa[!is.na(iso_group$iso),]
    iso_group <- iso_group[!is.na(iso_group$iso),]
    pa <- split_data(data, pa, iso_group$interaction, grouping_w_phylosip=F)
    pa <- base::lapply(pa, colSums, na.rm=T)
    # if pool unlabeled, combine frequencies across all unlabeled replicate groups
    # then replace frequencies in pa corresponding to those unlabeled groups
    pa <- do.call(rbind, pa)
    if(pool_unlabeled) {
      unlabeled <- rownames(pa) %in% iso_group$interaction[as.numeric(iso_group$iso)==1]
      pa_unlab <- pa[unlabeled,,drop=FALSE]
      pa_unlab <- colSums(pa_unlab, na.rm=T)
      pa[unlabeled,] <- rep(pa_unlab, each=sum(unlabeled))
    }
    pa <- ifelse(pa >= replicate, 1, 0)
    # regroup in order to match ft
    sf <- matrix(1L, nrow=length(sam_names), ncol=n_taxa) # sf = soft filter
    rownames(sf) <- sam_names
    groups <- levels(iso_group$interaction)
    for(i in 1:length(groups)) {
      relevant_samples <- as.character(iso_group$replicate[iso_group$interaction==groups[i]])
      sf[relevant_samples,] <- rep(pa[groups[i],], each=length(relevant_samples))
    }
    sf <- split_data(data, sf, rownames(sf), grouping_w_phylosip=FALSE)
    sf <- lapply(sf, c)
    # apply filter to WAD values
    ft <- base::Map('*', ft, sf) # or sf <- Map('*', ft, sf)
    if(calc_wvd) {
      wvd <- base::Map('*', wvd, sf)
    }
  }
  # 2. Hard filter
  # Unless taxa are removed beforehand, filtering here is hard filter on fractions only (i.e., replicate threshold is 1)
  if(filter && length(data@qsip@filter) > 0 && any(data@qsip@filter_levels$hard)) {
    ft <- do.call(rbind, ft)
    colnames(ft) <- phyloseq::taxa_names(data)
    ft <- ft[,colnames(ft) %in% data@qsip@filter]
    #
    if(calc_wvd) {
      wvd <- do.call(rbind, wvd)
      colnames(wvd) <- phyloseq::taxa_names(data)
      wvd <- wvd[,colnames(wvd) %in% data@qsip@filter]
    }
  } else if(filter && length(data@qsip@filter)==0 && any(data@qsip@filter_levels$hard)) {
    data <- filter_qsip(data, replicate=1)
    ft <- do.call(rbind, ft)
    colnames(ft) <- phyloseq::taxa_names(data)
    ft <- ft[,colnames(ft) %in% data@qsip@filter]
    #
    if(calc_wvd) {
      wvd <- do.call(rbind, wvd)
      colnames(wvd) <- phyloseq::taxa_names(data)
      wvd <- wvd[,colnames(wvd) %in% data@qsip@filter]
    }
  }
  # Regardless of hard filter, remove taxa that don't occur at all (all WADs=0) after applying soft filter
  if(filter && any(data@qsip@filter_levels$soft)) {
    #sf <- do.call(rbind, sf)
    ft <- do.call(rbind, ft)
    rownames(ft) <- sam_names
    labeled <- ft[as.character(iso_group$replicate[as.numeric(iso_group$iso)==2]),]
    colnames(ft) <- phyloseq::taxa_names(data)
    ft <- ft[, colSums(labeled) > 0] # or colSums(sf)
    data@qsip@filter <- colnames(ft)
    if(calc_wvd) {
      wvd <- do.call(rbind, wvd)
      rownames(wvd) <- sam_names
      colnames(wvd) <- phyloseq::taxa_names(data)
      wvd <- wvd[, colnames(ft)] # or colSums(sf)
    }
  # if filtering, use reduced tax name list
    tax_names <- data@qsip@filter
  # if not filtering, use full tax name list
  } else if(!filter) {
    tax_names <- phyloseq::taxa_names(data)
  }
  # return tax names
  # WAD values of 0 indicate no taxa present
  # if(class(ft)=='list') {
  #   ft <- base::lapply(ft, function(x) {x[x==0] <- NA; x})
  # } else ft[ft==0] <- NA
  # organize and add new data as S4 matrix
  data <- collate_results(data, ft, tax_names=tax_names, 'wad', sparse=TRUE)
  if(calc_wvd) {
    data <- collate_results(data, wvd, tax_names=tax_names, 'wvd', sparse=TRUE)
  }
  return(data)
}

# Utility functions to be used to do actual calculations

# Given vectors of x values and y values, calculate the weighted-average of the x-values (e.g., the weighted average density (WAD))
# With na.rm=F, any samples with NA density values will be returned with all taxa's WAD values as NA. (the default action)
# With na.rm=T, all fractions with NA densities will be ignored from a taxa's weighted average density (leading to lower numbers)
#
#
#     output = WAD.func(y, x)
#
#     y: vector of y-values (e.g., number of 16S copies)
#     x: vector of x-values (e.g., density of DNA)
#     -------------------------------------------------------
#     output: weighted-average of the x-values (single value)
#     Written by Ben Koch & Natasja van Gestel

wad <- function(y, x, na.rm=TRUE) {
  wad <- sum(x * (y / sum(y, na.rm=TRUE)), na.rm=na.rm)
  return(wad)
}
