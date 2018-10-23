#' Calculation of weighted average density differences
#'
#' Calculates differences in weighted average densities across replicates due to isotope incorporation
#'
#' @param data Data as a \code{phyloseq} object
#' @param filter Logical vector specifying whether or not to filter taxa from the weighted average density calculation.
#'   This will require \code{data} to have a filter applied with \code{\link{filter_qsip}}.
#'
#' @details Some details about proper isotope control-treatment factoring. If weighted average densities have not been calculated
#'   beforehand, \code{calc_d_wad} will compute those first.
#'
#' @return \code{calc_d_wad} adds two S4 Matrix objects to the \code{data@@qsip@@.Data} slot, one for differences
#'   in weighted average density, and the other for weighted average density values of light treatments only (to be used in
#'   future calculations). The row and column specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}},
#'   meaning if taxa are listed on the table rows, they will in the resulting S4 Matrix objects
#'
#' @seealso \code{\link{calc_wad}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate weighted average density differences
#'
#' @export

calc_d_wad <- function(data, filter=FALSE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  # if WAD values don't exist, calculate those first, this will also handle rep_id validity
  if(is.null(data@qsip[['wad']])) data <- calc_wad(data, filter=filter)
  if(length(data@qsip@rep_group)==0) stop('Must specify replicate groupings with rep_group')
  if(length(data@qsip@iso_trt)==0) stop('Must specify treatment and controls with iso_trt')
  # extract WAD values and convert to S3 matrix
  ft <- data@qsip[['wad']]
  ft <- as(ft, 'matrix')
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  n_taxa <- ncol(ft)
  tax_names <- colnames(ft)
  # manipulate data matrix and calculate
  # split by replicate groups, but keep track of light and heavy fractions
  ft <- valid_samples(data, ft, 'iso')
  iso_group <- ft[[2]]; ft <- ft[[1]]
  ft <- split_data(data, ft, iso_group$interaction, grouping_w_phylosip=F)
  # WAD values of 0 indicate no taxa presence in that replicate, convert to NA
  # so that mean WAD values are not pulled down
  ft <- base::lapply(ft, function(x) {x[x==0] <- NA; x})
  # calculate average WAD per taxa for each replicate group
  ft <- base::lapply(ft, colMeans, na.rm=T)
  # remove any NaNs resulting from when a taxon is missing in all replicates
  ft <- base::lapply(ft, function(x) {x[is.nan(x)] <- NA; x})
  # If there is no replicate grouping (i.e., all replicates in a treatment are grouped)...
  iso_group2 <- unique(iso_group[,!names(iso_group) %in% 'replicate']) # only get unique elements to match levels in ft
  iso_group2 <- iso_group2[match(names(ft), iso_group2$interaction),]
  if(isTRUE(all.equal(iso_group2$iso, iso_group2$grouping))) {
    d_ft <- ft[[2]] - ft[[1]]
    d_ft <- matrix(d_ft, nrow=1)
    rownames(d_ft) <- iso_group2$iso[as.numeric(iso_group2$iso)==2]
    keep_groups <- !logical(2)
  } else { # use a for-loop to subtract heavy from light fraction in each group
  # create a new list to add results of mean WAD difference into
  d_ft <- as.list(rep(0, nlevels(iso_group$grouping)))
  d_ft <- base::lapply(d_ft, matrix,
                       rep(0, n_taxa),
                       nrow=1,
                       ncol=n_taxa)
  names(d_ft) <- levels(iso_group$grouping)
  # For each repliate group: identify which elements of ft are light and which are heavy, then get difference
  keep_groups <- !logical(length(d_ft))
    for(i in 1:length(d_ft)) {
      # use numbers to reference non-labeled additions since they're element agnostic
      # any NA values result here when a taxa is completely missing from a heavy or light treatment in a replicate group
      which_light <- which(as.numeric(iso_group2$grouping)==i &
                             as.numeric(iso_group2$iso)==1)
      which_heavy <- which(as.numeric(iso_group2$grouping)==i &
                             as.numeric(iso_group2$iso)==2)
      # If there's no light OR no heavy treatment for a group of replicates, remove them
      if(length(which_light)==0 || length(which_heavy)==0) {
        warning('Unpaired isotope treatment in replicate group(s): ', names(d_ft)[1],
                '\nRemoving sample(s): ', paste(as.character(iso_group[iso_group$grouping==names(d_ft)[i], 'replicate']), collapse=', '),
                ' - from calculation', call.=FALSE)
        keep_groups[i] <- FALSE
        next
      }
      light <- ft[[which_light]]
      heavy <- ft[[which_heavy]]
      d_ft[[i]] <- heavy - light
    }
  d_ft <- d_ft[keep_groups]
  }
  # organize and add new data as S4 matrix
  data <- collate_results(data, d_ft, tax_names=tax_names, 'd_wad', sparse=TRUE)
  # return weighted average densities of light calcs only
  ft <- ft[iso_group2$grouping %in% levels(iso_group2$grouping)[keep_groups]] # remove unpaired & dropped groups
  iso_group2 <- iso_group2[iso_group2$grouping %in% levels(iso_group2$grouping)[keep_groups],] # remove unpaired & dropped groups
  wl <- ft[which(as.numeric(iso_group2$iso)==1)]
  if(isFALSE(all.equal(iso_group2$iso, iso_group2$grouping))) names(wl) <- unique(iso_group2$grouping)
  data <- collate_results(data, wl, tax_names=tax_names, 'wad_light', sparse=TRUE)
  return(data)
}
