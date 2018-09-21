#' Read QIIME 2 data
#'
#' Reads TSV tables of microbial features exported from QIIME 2's \code{biom} format
#'
#' @param file The name of the file which the feature data are to be read from. Each row of the table appears as one line of the file.
#'   If it does not contain an absolute path, the file name is relative to the current working directory, \code{getwd()}.
#'   Tilde-expansion is performed where supported.
#'
#'   This can be a compressed file (see \code{\link{file}}).
#' @param quiet Whether or not to produce messages notifying the user about file loading progress. Recommended for large tables
#'   (larger than 300 Mb, uncompressed).
#'
#' @details QIIME 2 uses the term features to refer to microbial taxonomic units in a way that is agnostic to the user's decision to use
#'   ASVs or OTUs. OTUs must first be clustered in QIIME 2, commonly at 97\% sequence similarity. ASVs represent unique sequences and so these tables will be
#'   larger than OTU tables for the same data set. QIIME 2 favors ASVs by default. For a discussion on the benefits of using ASVs over
#'   OTUs see: \url{https://www.nature.com/articles/ismej2017119} (referenced below).
#'
#'   The input table is created using \code{qiime tools export} using a feature table \code{qza} to create a \code{biom} file
#'   followed by \code{biom convert} to produce a tab-delimited feature table. \code{read_qiime2_table} uses \code{\link{scan}} to
#'   read in the data in two parts, first the column headers, and then the rest of the data.
#'
#' @return \code{read_qiime2_table} returns an integer matrix of microbial features (ASVs or OTUs) along
#'   rows and samples along columns. Feature ID codes are stored as row names, and sample codes as column names.
#'
#' @seealso \code{\link{read_qiime2_tax}}
#'
#' @examples
#' data(example_data)
#'
#' dat <- read_qiime2_table(data_name)
#'
#' dim(dat)
#' dat[1:10,1:10]
#'
#' #' @references Callahan, B.J., P.J. McMurdie, S.P. Holmes. 2017. Exact sequence variants should replace operational
#'   taxonomic units in marker-gene data analysis. \emph{ISME Journal} 11: 2639-2643.
#'
#' @export

read_qiime2_table <- function(file, quiet=FALSE) {
  if(is.null(file)) stop('Must provide file for input')
  sample_cols <- scan(file, what=character(), nlines=1, skip=1, sep='\t', quiet=T)
  sample_cols <- sample_cols[-1]
  data <- scan(file, what=character(), skip=2, sep='\t', quiet=quiet)
  data <- matrix(data, ncol=length(sample_cols) + 1, byrow=T)
  if(quiet==FALSE) message('Converting to integer matrix')
  rownames(data) <- data[,1]
  data <- data[,-1]
  storage.mode(data) <- 'integer'
  colnames(data) <- sample_cols
  return(data)
}

#' Read QIIME 2 data
#'
#' Reads TSV tables of microbial taxonomic information exported from QIIME 2's \code{biom} format
#'
#' @param file The name of the file which the feature data are to be read from. Each row of the table appears as one line of the file.
#'   If it does not contain an absolute path, the file name is relative to the current working directory, \code{getwd()}.
#'   Tilde-expansion is performed where supported.
#'
#'   This can be a compressed file (see \code{\link[base]{file}}).
# @param confidence Logical value indicating whether to keep any confidence values created by QIIME 2.
#' @param feature_type Single character indicating whether the table is composed of ASVs (amplicon sequence variants) or
#'   OTUs (operational taxonomic units). feature_type accepts only one argument.
#' @param fill_na Logical value indicating whether or not to return unclassified taxonomic values as blank (\code{fill_na=FALSE}) or NA.
#'
#' @details QIIME 2 uses the term features to refer to microbial taxonomic units in a way that is agnostic to the user's decision to use
#'   ASVs or OTUs. OTUs must first be clustered in QIIME 2, commonly at 97\% sequence similarity. ASVs represent unique sequences and so these tables will be
#'   larger than OTU tables for the same data set. QIIME 2 favors ASVs by default. For a discussion on the benefits of using ASVs over
#'   OTUs see: \url{https://www.nature.com/articles/ismej2017119} (referenced below).
#'
#'   The input tab-separated TSV table is created using \code{qiime tools export} using a taxonomy \code{qza} file.
#'
#' @return \code{read_qiime2_tax} returns a character matrix of taxonomic assignments from Kingdom to Species, identified to each
#'   microbial taxa by the row names of the matrix.
#'
#' @seealso \code{\link{read_qiime2_table}}, \code{\link{specify_unclassifieds}}
#'
#' @examples
#' data(example_data)
#'
#' dat <- read_qiime2_tax(data_name, feature_type='otu', fill_na=TRUE)
#' head(dat)
#'
#' dat <- read_qiime2_tax(data_name, feature_type='otu', fill_na=FALSE)
#' head(dat)
#'
#' @references Callahan, B.J., P.J. McMurdie, S.P. Holmes. 2017. Exact sequence variants should replace operational
#'   taxonomic units in marker-gene data analysis. \emph{ISME Journal} 11: 2639-2643.
#'
#' @export

read_qiime2_tax <- function(file, confidence=TRUE, feature_type=c('OTU', 'ASV'), fill_na=TRUE) {
  if(is.null(file)) stop('Must provide file for input')
  feature_type <- match.arg(toupper(feature_type), c('OTU', 'ASV'))
  featureName <- paste0(feature_type, '_ID')
  data <- read.table(file, sep='\t', header=T, stringsAsFactors=F)
  if(fill_na==FALSE) {
    tax <- read.table(text=as.character(data$Taxon), sep=';', fill=T, stringsAsFactors=F)
    tax <- sapply(tax, function(x) sub('\\s?\\w{1}__', '', x, perl=T))
  } else {
    na.strings <- c('', ' s__', ' g__', ' f__', ' o__', ' c__', ' p__', ' k__')
    tax <- read.table(text=as.character(data$Taxon), sep=';', fill=T, stringsAsFactors=F, na.strings=na.strings)
    tax <- sapply(tax, function(x) sub('\\s?\\w{1}__', '', x, perl=T))
  }
  colnames(tax) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus' ,'Species')
#  if(confidence) {
#    data <- data.frame(data[,1], tax, Confidence=data[,3], stringsAsFactors=F)
#  } else {
#    data <- data.frame(data[,1], tax, stringsAsFactors=F)
#  }
  rownames(tax) <- data[,1]
#  names(data)[1] <- featureName
#  return(data)
  return(tax)
}

#' Read QIIME 2 data
#'
#' Reads TSV or txt files used as QIIME 2 metadata
#'
#' @param file The name of the file(s) which the data are to be read from. After ignoring any redundant headers, each row of the table appears
#'   as one line from the file(s). If there is no absolute path, the file name(s) will be relative to the current working directory, \code{getwd()}.
#'   Tilde-expansion is performed where supported.
#'
#'   This can be a compressed file or files (see \code{\link[base]{file}}).
#' @param barcodes Logical value indicating whether or not to keep barcodes and linker sequences.
#' @param stringsAsFactors Logical: should the character vector be converted to a factor? This tag excludes the first column of sample ID names and any possible
#'   description columns (see details below).
#' @param run_names Character vector specifying an identifying name for each run to be read in. This will create a \code{Run} column in the resulting data frame
#'
#' @details Some details
#'
#' @return \code{read_qiime2_metadata} returns a data frame from either a single or multiple QIIME2-formatted metadata tables.
#'
#' @examples
#'
#' @export

read_qiime2_metadata <- function(file, barcodes=FALSE, stringsAsFactors=TRUE, run_names=c()) {
  if(is.null(file)) stop('Must provide file for input')
  data <- vector(mode='list', length(file))
  for(i in 1:length(file)) {
    col_names <- scan(file[i], what=character(), nlines=1, sep='\t', quiet=T)
    info <- scan(file[i], what=character(), skip=2, sep='\t', quiet=T)
    info <- matrix(info, ncol=length(col_names), byrow=T)
    colnames(info) <- col_names
    if(is.null(run_names)==FALSE) {
      if(length(data)!=length(run_names)) stop('The number of sequencing runs to read and the length of the naming vector must be equal')
      run_names <- as.character(run_names)
      info <- cbind(info, run_names[i])
      colnames(info)[ncol(info)] <- 'Run'
    }
    if(stringsAsFactors==FALSE) {
      data[[i]] <- as.data.frame(info, stringsAsFactors=F)
    } else {
      data[[i]] <- as.data.frame(info, stringsAsFactors=T)
      data[[i]][,1] <- as.character(data[[i]][,1])
      try(data[[i]][,'Description'] <- as.character(data[[i]][,'Description']), silent=T)
    }
  }
  if(length(file) > 1) {
    column_names <- lapply(data, colnames)
    names_equal <- Reduce(identical, column_names)
    if(names_equal==TRUE) {
      data <- do.call(rbind, data)
    } else {
      warning('Column names are not equal across all files. Merging will produce missing data', call.=F)
      data <- Reduce(function(x, y) merge(x, y, all=T, sort=F), data)
    }
  } else data <- data[[1]]
  if(barcodes==FALSE) {
    data <- data[,!names(data) %in% c('BarcodeSequence', 'LinkerPrimerSequence')]
  }
  rownames(data) <- data[,1]
  return(data)
}

# function that changes designation of unclassified sequences to 'Unclassified x' or 'Unidentified x' depending on user input
# can either use the lowest classifcation available, or some set level (perhaps both, like maybe ID to family, but use next best if family is unavailable)
# arguments: x, use_lowest, use_same) {
#
#}
