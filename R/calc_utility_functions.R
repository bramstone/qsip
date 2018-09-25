# Collection of functions to make calculation function code neater

# Calculation of 16S copies per taxon from relative abundances and 16S copy numbers
# Note, will return a feature table with taxa as columns
copy_no <- function(data) {
  ft <- as(data@otu_table, 'matrix')
  # make relative abundances with taxa as columns and samples as rows
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  ft <- ft / base::colSums(ft, na.rm=TRUE)
  # multiply by total 16S copy number per sample
  ft <- ft * data@sam_data[[data@qsip@abund]]
  return(ft)
}

# Function to generate small reference data frame of unique isotope and isotope x grouping combinations
# to use with grouped calculations where the light and heavy treatments must be identified
iso_grouping <- function(data, iso, grouping) {
  output <- data.frame(iso=data$sam_dat[[iso]],
                       grouping=data$sam_data[[grouping]])
  output$interaction <- interaction(output$iso, output$grouping)
  output <- unique(output[,1:3])
  return(output)
}

# Function used to split qSIP data into list of sub-matrices
split_data <- function(data, new_data, grouping, grouping_w_phylosip=T) {
  if(grouping_w_phylosip) grouping <- data@sam_data[[grouping]]
  new_data <- split(new_data, grouping)
  new_data <- base::lapply(new_data, matrix,
                     byrow=FALSE,
                     ncol=phyloseq::ntaxa(data))
  return(new_data)
}

# Function used to handle adding new data to phylosip object .Data slot
collate_results <- function(data, new_data, metric) {
  # combine format based on whether taxa were rows or not
  if(phyloseq::taxa_are_rows(data)) {
    new_data <- do.call(cbind, new_data)
  } else new_data <- do.call(rbind, new_data)
  # add feature names back in (replicate names automatically utilized from split)
  if(is.null(rownames(new_data))) {
    rownames(new_data) <- phyloseq::taxa_names(data)
  } else {
    colnames(new_data) <- phyloseq::taxa_names(data)
  }
  # convert to S4 Matrix which is more memory efficient
  new_data <- Matrix::Matrix(new_data)
  # add wad values to data slot of qSIP portion of object
  if(any(attributes(data@qsip)$names %in% metric)) { # if wad alreay exists, replace
    warning('Overwriting existing ', metric, ' values', call.=FALSE)
    replace_num <- which(attributes(data@qsip)$names %in% metric)
    data@qsip@.Data[[replace_num]] <- new_data
  }
  else { # else append it to the list and update the list names
    data@qsip@.Data[[length( data@qsip@.Data) + 1]] <- new_data
    attributes(data@qsip)$names <- c(attributes(data@qsip)$names, metric)
  }
  return(data)
}