# Collection of functions to make calculation function code neater

# Calculation of 16S copies per taxon from relative abundances and 16S copy numbers
# Note, will return a feature table with taxa as columns
copy_no <- function(data) {
  ft <- as(data@otu_table, 'matrix')
  # make relative abundances with taxa as columns and samples as rows
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  ft <- ft / base::rowSums(ft, na.rm=TRUE)
  # multiply by total 16S copy number per sample
  ft <- ft * data@sam_data[[data@qsip@abund]]
  ft[is.nan(ft)] <- 0
  return(ft)
}

# Function to generate small reference data frame of unique isotope and isotope x grouping combinations
# to use with grouped calculations where the light and heavy treatments must be identified
iso_grouping <- function(data, iso, rep_id, grouping) {
  if(length(data@qsip@rep_group)==0) grouping <- 1 else grouping <- data@sam_data[[grouping]]
  output <- data.frame(iso=data@sam_data[[iso]],
                       replicate=data@sam_data[[rep_id]],
                       grouping=grouping)
  #if(length(data@qsip@rep_group)==0) {
 #   output$interaction <- output$iso
 #   output$full.interaction <- output$interaction
 # } else {
    output$interaction <- interaction(output$iso, output$grouping)
    output$full.interaction <- interaction(output$replicate, output$interaction)
#  }
  output <- output[!duplicated(output$full.interaction) | is.na(output$full.interaction),]
  output$full.interaction <- NULL
  rownames(output) <- NULL
  return(output)
}

# Function to generate small reference data frame of unique replicate and timepoint combinations
# to use with pop calculations where the time 0 and time t must be identified
time_grouping <- function(data, timepoint, rep_id, grouping) {
  if(length(data@qsip@rep_group)==0) grouping <- 1 else grouping <- data@sam_data[[grouping]]
  output <- data.frame(time=data@sam_data[[timepoint]],
                       replicate=data@sam_data[[rep_id]],
                       grouping=grouping)
  #if(length(data@qsip@rep_group)==0) {
  #  output$interaction <- output$time
  #  output$full.interaction <- output$interaction
  #} else {
    output$interaction <- interaction(output$time, output$grouping)
    output$full.interaction <- interaction(output$replicate, output$interaction)
  #}
  output <- output[!duplicated(output$full.interaction) | is.na(output$full.interaction),]
  output$full.interaction <- NULL
  rownames(output) <- NULL
  return(output)
}

# Function used to split qSIP data into list of sub-matrices
split_data <- function(data, new_data, grouping, grouping_w_phylosip=T, keep_names=0) {
  if(grouping_w_phylosip) grouping <- data@sam_data[[grouping]]
  if(keep_names==1) {
    sam_names <- rownames(new_data)
    sam_names <- split(sam_names, grouping)
  } else if(keep_names==2) {
    sam_names <- rownames(new_data)
    sam_names <- split(sam_names, grouping)
  }
  n_taxa <- ncol(new_data)
  new_data <- split(new_data, grouping)
  new_data <- base::lapply(new_data, matrix,
                     byrow=FALSE,
                     ncol=n_taxa)

  if(keep_names==1) {
    new_data <- base::Map(function(x,y) {rownames(x) <- y; x}, new_data, sam_names)
  } else if(keep_names==2) {
    new_data <- base::Map(function(x,y) {colnames(x) <- y; x}, new_data, sam_names)
  }
  return(new_data)
}

# Function used to handle adding new data to phylosip object .Data slot (data@qsip@.Data)
# parameter ... Indicates options to pass to Matrix, for specifying whether it should be sparse or not
collate_results <- function(data, new_data, tax_names=NULL, metric, ...) {
  if(is.null(tax_names)) {
    if(ncol(new_data)==length(phyloseq::taxa_names(data))) {
      tax_names <- phyloseq::taxa_names(data)
    } else if(length(data@qsip@filter) > 0) {
      tax_names <- data@qsip@filter
    }
  }
  # combine format based on whether taxa were rows or not
  if(class(new_data)=='list') {
    if(phyloseq::taxa_are_rows(data)) {
      new_data <- do.call(cbind, new_data)
      if(is.null(rownames(new_data))) {
        rownames(new_data) <- tax_names
      }
    } else {
      new_data <- do.call(rbind, new_data)
      if(is.null(colnames(new_data))) {
        colnames(new_data) <- tax_names
      }
    }
    new_data <- Matrix::Matrix(new_data, ...)
    # add feature names back in (replicate names automatically utilized from split)
  } else if(class(new_data)=='matrix') {
    if(phyloseq::taxa_are_rows(data)) new_data <- t(new_data)
    if(is.null(rownames(new_data))) {
      rownames(new_data) <- tax_names
    } else if(is.null(colnames(new_data))) {
      colnames(new_data) <- tax_names
    }
    # convert to S4 Matrix which is more memory efficient
    new_data <- Matrix::Matrix(new_data, ...)
  } else if(class(new_data)=='numeric' && is.null(names(new_data))) names(new_data) <- tax_names
  # add wad values to data slot of qSIP portion of object
  if(any(attributes(data@qsip)$names %in% metric)) { # if wad alreay exists, replace
    warning('Overwriting existing ', metric, ' values', call.=FALSE)
    replace_num <- which(attributes(data@qsip)$names %in% metric)
    data@qsip@.Data[[replace_num]] <- new_data
  }
  else { # else append it to the list and update the list names
    data@qsip@.Data[[length(data@qsip@.Data) + 1]] <- new_data
    attributes(data@qsip)$names[length(data@qsip@.Data)] <- metric
  }
  return(data)
}

# function that generates confidence intervals for a metric (lower CI, median, upper CI)
# also handles column names. Returns a list
summarize_ci <- function(bootstraps, ci, grouping, ncols, list_names=c('ci_l', 'med', 'ci_u'), data) {
  summaries <- t(apply(bootstraps, 1,
                       quantile,
                       c((1 - ci)/2, .5, (1 - ci)/2 + ci),
                       na.rm=TRUE))
  ci_l <- summaries[,1]
  med <- summaries[,2]
  ci_u <- summaries[,3]
  # arrange into matrices
  ci_l <- matrix(ci_l, ncol=ncols, byrow=T)
  med <- matrix(med, ncol=ncols, byrow=T)
  ci_u <- matrix(ci_u, ncol=ncols, byrow=T)
  # assign rownames (colnames assigned in collate_results)
  if(length(data@qsip@rep_group)==0) {
    rownames(ci_l) <- rownames(med) <- rownames(ci_u) <- levels(grouping[,1])[2:nlevels(grouping[,1])]
  } else {
    rnames <- expand.grid(levels(grouping$grouping),
                          levels(grouping[,1])[2:nlevels(grouping[,1])],
                          stringsAsFactors=FALSE)
    rownames(ci_l) <- rownames(med) <- rownames(ci_u) <- interaction(rnames[,1], rnames[,2], sep=':')
  }
  # put into list and return
  output <- list(ci_l, med, ci_u)
  names(output) <- list_names
  return(output)
}

# function that removes invalid (missing) rows between matrix (ft, WADS, diff_WADS, etc.) and grouping data.frame
# returns a list of 2; [[1]] is the matrix, [[2]] is the grouping data.frame with invalid rows removed
# samples should be rows in feature table
valid_samples <- function(data, feature_table, grouping=c('iso', 'time'), quiet=FALSE, match_replicate=FALSE) {
  grouping <- match.arg(grouping, c('iso', 'time'))
  # grouping by isotope or timepoint?
  if(grouping=='iso') {
    group_data <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
  } else {
    group_data <- time_grouping(data, data@qsip@timepoint, data@qsip@rep_id, data@qsip@rep_group)
  }
  # if matching replicates, identify any samples that are unmatched between time 0 and time t
  if(match_replicate) {
    rep_nums <- unique(as(data@sam_data[,c(data@qsip@rep_id, data@qsip@rep_num)], 'data.frame'))
    names(rep_nums) <- c('replicate', 'replicate_num')
    matching <- merge(group_data, rep_nums, all.x=TRUE)
    matching <- split(matching, interaction(matching$grouping, matching$replicate_num))
    matching <- lapply(matching, function(x) {if(nrow(x)==1) x$replicate_num <- NA; x})
    matching <- do.call(rbind, matching)
    group_data <- merge(group_data, matching[,c('replicate', 'replicate_num')], all=TRUE)
    group_data <- group_data[order(group_data$interaction, group_data$replicate_num),]
  }
  # identify any NAs that occur across iso/timepoint column, replicate IDs, or grouping variable
  invalid <- sapply(group_data, is.na)
  invalid <- apply(invalid, 1, any)
  # issue warning of samples lost
  if(sum(invalid) > 0 && quiet==FALSE) {
    warning('Dropping sample(s): ',
            paste(as.character(group_data$replicate[invalid]), collapse=', '),
            ' - from calculation', call.=FALSE)
  }
  # match matrix rows with samples with complete data
  group_data <- group_data[!invalid,]
  # re-level factors
  group_data <- data.frame(lapply(group_data, factor))
  feature_table <- feature_table[match(group_data$replicate, rownames(feature_table)),]
  # return data
  return(list(feature_table, group_data))
}


# function which recombines randomly subsampled matrix in the correct order so that downstream calculations will
# be performed on the correct subsets
recombine_in_order <- function(ft, grouping, n_taxa, condensed_grouping=FALSE) {
  if(condensed_grouping) {
    # pop value abundances are summarized, reducing number of samples of each group to 1
    grouping <- grouping[!duplicated(grouping$interaction),]
    group_counts <- nlevels(grouping$interaction)
  } else group_counts <- sapply(ft, nrow)
  which_group <- as.numeric(grouping$interaction)
  temp <- matrix(nrow=sum(group_counts), ncol=n_taxa)
  for(i in 1:nlevels(grouping$interaction)) {
    temp[which_group==i,] <- ft[[i]]
  }
  return(temp)
}
