#' Melt qSIP data to long form
#'
#' Converts data in phylosip object's qsip slot from several matrices to single, long-form data.frame
#'
#' @param data Data as a \code{phylosip} object
#' @param taxonomy Logical tag specifying whether or not to include taxonomic information of each feature. Default is \code{TRUE}.
#' @param abundance Logical tage specifying whether or not to include sequencing abundance of each feature. Default is \code{TRUE}.
#' @param relativize Logical tag specifying whether or not to present feature abundances as relative from 0 to 1 or to use those in @@otu_table slot as is.
#' @param exclude Character vector specifying any qSIP measures to exclude from the returned data frame.
#'
#' @details Details here on different specifications
#'
#'   Note that WAD calculations will always be performed on a per-replicate basis. Thus, when mapped against grouped calculations, duplicate values
#'   will exist for non-WAD measures. To avoid this issue, exclude WAD calculations from the output with \code{exlude='wad'}.
#'
#' @return \code{qsmelt} returns a single data frame with all qSIP related values for each taxon in each replicate, or group of replicates.
#'   Any replicate grouping (specified by \code{data@@qsip@@rep_group}) or distinct timepoints (specified by \code{data@@qsip@@timepoint}) will
#'   be represented in separate columns.
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate weighted average densities
#'
#'  # Melt data
#'
#' @export

qsmelt <- function(data, taxonomy=FALSE, abundance=FALSE, relativize=FALSE, exclude=c()) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  if(length(data@qsip@.Data)==0) stop('No values to combine')
  qsip <- as.list(data@qsip@.Data)
  names(qsip) <- names(data@qsip)
  # transform to long-format data frames
  for(val in names(qsip)) {
    # make sure taxa are rows
    if(!is.numeric(qsip[[val]])) {
      if(!phyloseq::taxa_are_rows(data)) qsip[[val]] <- t(qsip[[val]])
      # convert numeric vectors to dgCMatrix
    } else {
      qsip[[val]][is.na(qsip[[val]])] <- 0
      qsip[[val]] <- Matrix::Matrix(qsip[[val]], sparse=TRUE, dimnames=list(names(qsip[[val]]), NULL))
    }
    # convert all to dgTMatrix, which stores indexing information more intuitively
    qsip[[val]] <- as(qsip[[val]], 'dgTMatrix')
    # construct data frames
    if(ncol(qsip[[val]]) > 1) {
      # if tax ID and rep ID exist
      qsip[[val]] <- data.frame(tax_id=qsip[[val]]@Dimnames[[1]][qsip[[val]]@i + 1],
                              rep_id=qsip[[val]]@Dimnames[[2]][qsip[[val]]@j + 1],
                              value=qsip[[val]]@x,
                              stringsAsFactors=TRUE)
    } else {
      # if only tax ID exists
      qsip[[val]] <- data.frame(tax_id=qsip[[val]]@Dimnames[[1]][qsip[[val]]@i + 1],
                              value=qsip[[val]]@x,
                              stringsAsFactors=TRUE)
    }
    # change "value" column name
    names(qsip[[val]])[which(names(qsip[[val]])=='value')] <- val
  }
  # accommodate grouped values
  if(length(data@qsip@rep_group)!=0 | length(data@qsip@timepoint)!=0) {
    # get data frame of replicate IDs mapped to replicate groups
    rep_map <- unique(data@sam_data[,c(data@qsip@rep_id, data@qsip@rep_group)])
    rep_map <- as(rep_map, 'data.frame')
    rep_map <- rep_map[complete.cases(rep_map),]
    # identify when group IDs are present vs. sample IDs, add additional columns if rep_group and timepoint are present
    for(val in names(qsip)) {
      # Do sample IDs from a qSIP data frame match with those of the whole dataset?
      # If not, then they are likely group IDs. Rename them and try to map rep_group and/or timepoint
      # Ignore single vector values (i.e., ones that don't have a replicate ID associated with them)
      if(!any(unique(qsip[[val]]$rep_id) %in% rep_map[,data@qsip@rep_id]) && !is.null(qsip[[val]]$rep_id)) {
        names(qsip[[val]])[which(names(qsip[[val]])=='rep_id')] <- 'group_id'
        # if it's the case that they're group IDs, need to map them to the appropriate groups
        if(length(data@qsip@rep_group)!=0) {
          # use regexpr to return match start point, and endpoint, if such exists
          rep_group_levels <- paste(unique(rep_map[,data@qsip@rep_group]), collapse='|')
          rep_group_ids <- regexpr(rep_group_levels, qsip[[val]]$group_id)
          # use substr to return matching portions, NA for blank response
          qsip[[val]]$rep_group <- substr(qsip[[val]]$group_id, rep_group_ids, attributes(rep_group_ids)$match.length)
          qsip[[val]]$rep_group[nchar(qsip[[val]]$rep_group)==0] <- NA
          names(qsip[[val]])[which(names(qsip[[val]])=='rep_group')] <- data@qsip@rep_group
        }
        if(length(data@qsip@timepoint)!=0) {
          #
          timepoint_levels <- paste(unique(rep_map[,data@qsip@timepoint]), collapse='|')
          timepoint_ids <- regexpr(timepoint_levels, qsip[[val]]$group_id)
          #
          qsip[[val]]$timepoint <- substr(qsip[[val]]$group_id, timepoint_ids, attributes(timepoint_ids)$match.length)
          qsip[[val]]$timepoint[nchar(qsip[[val]]$timepoint)==0] <- NA
          names(qsip[[val]])[which(names(qsip[[val]])=='timepoint')] <- data@qsip@timepoint
        }
        # if instead you used replicate IDs the whole way through, you can still get grouping info
      } else if(any(unique(qsip[[val]]$rep_id) %in% rep_map[,data@qsip@rep_id]) && !is.null(qsip[[val]]$rep_id)) {
        names(qsip[[val]])[which(names(qsip[[val]])=='rep_id')] <- data@qsip@rep_id
        qsip[[val]] <- merge(qsip[[val]], rep_map, all.x=TRUE)
      }
    }
  }
  # combine into single data frame
  if(length(exclude) > 0) qsip <- qsip[!names(qsip) %in% exclude]
  comb_qsip <- Reduce(function(x, y) merge(x, y, all=TRUE), qsip)
  # identify if all numeric values are NA, remove those rows
  all_nas <- apply(comb_qsip[,sapply(comb_qsip, is.numeric)], 1, function(x) all(is.na(x)))
  comb_qsip <- comb_qsip[!all_nas,]
  rownames(comb_qsip) <- NULL
  # include taxonomy
  if(taxonomy) {
    tax <- as(data@tax_table, 'matrix')
    tax <- data.frame(tax_id=rownames(tax), tax)
    rownames(tax) <- NULL
    for(i in 1:ncol(tax)) attr(tax[,i], 'names') <- NULL
    comb_qsip <- merge(comb_qsip, tax, all.x=TRUE)
  }
  if(abundance) {
    # get 16S abundances
    ft <- copy_no(data)
    tax_names <- colnames(ft)
    # split by replicate IDs
    ft <- split_data(data, ft, data@qsip@rep_id)
    # limit to samples in qSIP data
    if(!is.null(comb_qsip[,data@qsip@rep_id])) ft <- ft[match(unique(comb_qsip[,data@qsip@rep_id]), names(ft))]
    # calculate totals
    if(relativize) {
      ft <- lapply(ft, function(x) colSums(x) / sum(x))
    } else {
      ft <- lapply(ft, colSums)
    }
    # recombine and remove taxa
    ft <- do.call(rbind, ft)
    ft <- ft[,match(unique(comb_qsip$tax_id), tax_names)]
    colnames(ft) <- unique(comb_qsip$tax_id)
    # convert to data frame
    ft <- Matrix::Matrix(ft, sparse=TRUE)
    ft <- as(ft, 'dgTMatrix')
    ft <- data.frame(rep_id=ft@Dimnames[[1]][ft@i + 1],
                     tax_id=ft@Dimnames[[2]][ft@j + 1],
                     abund=ft@x,
                     stringsAsFactors=TRUE)
    # rename
    names(ft)[which(names(ft)=='rep_id')] <- data@qsip@rep_id
    # if comb_qsip doesn't have a rep_id column (i.e., all groups), need to aggregate
    if(is.null(comb_qsip[,data@qsip@rep_id])) {
      ft <- merge(ft, rep_map, all.x=TRUE)
      # if you only have data at different replicate groups (no timepoints)
      if(is.null(ft[,data@qsip@timepoint])) {
        ft <- aggregate(ft, list(ft$tax_id,
                                 ft[,data@qsip@rep_group]),
                        mean)
        names(ft) <- c(data@qsip@tax_id, data@qsip@rep_id, 'abund')
        # or if you only have data at different timepoints (no replicate groups)
      } else if(is.null(ft[,data@qsip@rep_group])) {
        ft <- aggregate(ft, list(ft$tax_id,
                                 ft[,data@qsip@timepoint]),
                        mean)
        names(ft) <- c(data@qsip@tax_id, data@qsip@timepoint, 'abund')
        # or if you have data for both
      } else if(!is.null(ft[,data@qsip@rep_group]) & is.null(!ft[,data@qsip@timepoint])) {
        ft <- aggregate(ft, list(ft$tax_id,
                                 ft[,data@qsip@rep_group],
                                 ft[,data@qsip@timepoint]),
                        mean)
        names(ft) <- c(data@qsip@tax_id, data@qsip@rep_group, data@qsip@timepoint, 'abund')
      }
    }
    comb_qsip <- merge(comb_qsip, ft, all.x=TRUE)
  }
  # return final data frame
  return(comb_qsip)
}
