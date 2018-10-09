#' Specify multiple taxa filting options for qSIP calculations
#'
#' Creates a data frame of per-replicate and per-fraction minimum incidence frequencies to retain microbial taxa in the dataset
#'
#' @param replicate Numeric vector specifying the minimum frequency of occurrence of microbial taxa across replicates.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the replicate level
#' @param fraction Numeric vector specifying the minimum frequency of occurrence of microbial taxa across fractions within a sample.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the fraction level
#' @param rm_combn Optional character vector specifying any particular combinations of replicate and fraction frequency to remove.
#'   Replicate and frequency combinations should be specified by separation with \code{:} (\emph{e.g.}, \code{'3:12'})
#'
#' @details \code{create_filters} should be used to specify the desired frequencies of microbial taxa \emph{prior} to calculation of
#'   atom excess fraction / atom percent excess. Filters must also be specified so that diagnostic functions can be run as well.
#'   Because \code{create_filters} uses \code{\link[base]{expand.grid}}, all possible combinations of the first two arguments will
#'   be used to construct the resulting data frame If certain combinations are not desirable, the \code{rm_combn} argument can be
#'   supplied to remove them.
#'
#' @return \code{create_filters} produces a data frane of replicate frequencies and within-replicate-fraction frequencies to be
#'   investigated in downstream analyses.
#'
#' @examples
#' # Only filter on samples
#'  create_filters(2:3)
#'
#' # Filter on samples and fractions
#' create_filters(2:3, 10:12)
#'
#' # remove certain combinations
#' create_filters(2:3, 10:13, c('3:13', '2:10'))
#'
#' @export

create_filters <- function(replicate=0, fraction=0, rm_combn=character()) {
  filters <- expand.grid(replicate, fraction)
  names(filters) <- c('rep_freq', 'frac_freq')
  filters <- filters[order(filters$rep_freq),]
  if(!missing(rm_combn)) {
    rm_combn <- strsplit(rm_combn, ':')
    rm_combn <- do.call(rbind, rm_combn)
    storage.mode(rm_combn) <- 'integer'
    for(i in 1:nrow(rm_combn)) {
      filters <- filters[filters$rep_freq!=rm_combn[i,1] |
                           filters$frac_freq!=rm_combn[i,2],]
    }
  }
  rownames(filters) <- NULL
  return(filters)
}

#' Creates filtering criteria
#'
#' Returns binary vector specifying whether taxa from a phylosip object satisfy minimum frequency specifications
#'
#' @param data \code{Phylosip}-class object to pull feature taxa table from.
#' @param replicate Numeric vector specifying the minimum frequency of occurrence of microbial taxa across replicates.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the replicate level
#' @param fraction Numeric vector specifying the minimum frequency of occurrence of microbial taxa across fractions within a sample.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the fraction level
#' @param code Optional character vector specifying a particular combination of replicate and fraction frequency to test.
#'   Replicate and frequency combinations should be specified by separation with \code{:} (\emph{e.g.}, \code{'3:12'})
#'
#' @details \code{impose_filter} is primarily utilized within other functions.
#'
#' @return Returns a
#'
#' @seealso \code{\link{create_filters}}
#'
#' @examples
#' # Filter a tax table

impose_filter <- function(data, replicate=0, fraction=0, code=character()) {
  # if !is.null(code) parse and use code
  if(!is.null(code)) {
    if(length(code) > 1) warning('More than one filter code provided; will only use ', code[1], call.=F)
    code <- strsplit(code, ':')[[1]]
    replicate <- code[1]
    fraction <- code[2]
  }
  # extract feature table and convert to matrix with taxa as columns
  ft <- as(data@otu_table, 'matrix')
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  # make presence-absence
  ft <- ceiling(ft / max(ft))
  storage.mode(ft) <- 'integer'
  # split by replicates
  ft <- split_data(data, ft, data@qsip@rep_id)
  # calculate within replicate (i.e., fraction frequency)
  ft <- lapply(ft, colSums, na.rm=T)
  # apply filter, combine
  ft <- lapply(ft, function(x) ifelse(x >= fraction, 1, 0))
  ft <- do.call(rbind, ft)
  # split by replicate group
  iso_group <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
  # Drop any rows (probably NA) that don't appear in ft rownames
  iso_group <- iso_group[match(rownames(ft), iso_group$replicate),]
  ft <- ft[!is.na(iso_group$iso),]
  iso_group <- iso_group[!is.na(iso_group$iso),]
  ft <- split_data(data, ft, iso_group$interaction, grouping_w_phylosip=F)
  # apply filter, combine
  ft <- lapply(ft, colSums, na.rm=T)
  ft <- lapply(ft, function(x) ifelse(x >= replicate, 1, 0))
  ft <- do.call(rbind, ft)
  # add in taxa names and calculate presence or absence after filters
  colnames(ft) <- phyloseq::taxa_names(data)
  output <- colSums(ft)
  output <- ifelse(output > 0, 1L, 0L)
}

#' Explore different filtering levels
#'
#' Generates an S4 \code{filter} class object which stores filtering statistics
#'
#' @param data \code{Phylosip}-class object to pull feature taxa table from.
#' @param filters Optional data frame specifying filtering levels to explore.
#'   The default is to utilize specifications from \code{data@@qsip@@filter}.
#'   If \code{NULL}, \code{explore_filters} will generate statistics for all possible frequencies at both the replicate and fraction level.
#'
#' @details Some text here
#'
#' @return Some text here
#'
#' @seealso \code{\link{create_filters}}
#'
#' @examples
#' # Filter a tax table
#'
#' @export

explore_filters <- function(data, filters=data@qsip@filter) {
  # if no filters provided, will use all
  if(is.null(filters)) {
    # number of fractions in data?
    fractions <- table(data@sam_data[[data@qsip@rep_id]])
    # group by isotope treatment, replicate ID, replicate group
    group <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
    group <- group[match(names(fractions), iso_group$replicate),]
    replicates <- split(fractions, group$interaction) # will drop NA isotope values
    # number of replicates in data?
    relicates <- sapply(replicates, length)
    filters <- create_filters(1:max(replicates), 1:max(fractions))
  }
  # create blank matrix (taxa as rows, number of filter combos as columns)
  output <- matrix(0L, nrow=phyloseq::ntaxa(data),
                   ncol=nrow(filters),
                   dimnames=list(phyloseq::tax_names(data0),
                                 interaction(filters, sep=':')))
  # for every combination, run impose filters, stick output to column i of matrix
  for(i in 1:nrow(filters)) {
    output[,i] <- impose_filter(data,
                                replicate=filters$rep_freq[i],
                                fraction=filters$frac_freq[i])
  }
  # convert to S4 Matrix class
  output <- Matrix::Matrix(output, sparse=T)
  # add other statistics into output
  # NEED TO CREATE NEW S4 CLASS FOR FILTERING STATISTICS
  # NEED TO UPDATE SET_CLASS AND SET_INIT_METHODS SCRIPTS
  # maybe one of the slots could be an extension of the input filtering data frame
  #  with columns for the proportion (or number) of remaining taxa and another for remaining samples
}

#' Filter taxa
#'
#' Filters taxa from qSIP portion of phyloseq object
#'
#' @param data \code{Phylosip}-class object to pull feature taxa table from.
#' @param replicate Numeric vector specifying the minimum frequency of occurrence of microbial taxa across replicates.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the replicate level
#' @param fraction Numeric vector specifying the minimum frequency of occurrence of microbial taxa across fractions within a sample.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the fraction level
#' @param code Optional character vector specifying a particular combination of replicate and fraction frequency to test.
#'   Replicate and frequency combinations should be specified by separation with \code{:} (\emph{e.g.}, \code{'3:12'})
#'
#' @details \code{impose_filter} is primarily utilized within other functions.
#'
#' @return Returns a
#'
#' @seealso \code{\link{create_filters}}, \code{\link{explore_filters}}
#'
#' @examples
#' # Filter a phyloseq object
#'
#' @export

filter_qsip <- function(data, replicate=0, fraction=0, code=character()) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  tax_filter <- impose_filter(data, replicate=replicate, fraction=fraction, code=code)
  tax_filter <- names(tax_filter[tax_filter > 0])
  types <- sapply(data@qsip@.Data, function(x) class(x)[1])
  # remove taxa from each qSIP-related output object in the data@qsip@.Data list
  for(i in 1:length(dat@qsip)) {
    x <- data@qsip@.Data[[1]]
    if(types=='dgCMatrix') {
      if(phyloseq::taxa_are_rows(data)) {
        x <- x[rownames(x) %in% tax_filter,]
      }
      else {
        x <- x[,colnames(x) %in% tax_filter]
      }
    } else if(types=='numeric') {
      x <- x[names(x) %in% tax_filter]
    }
    data@qsip@.Data[[1]] <- x
  }
  return(data)
}
