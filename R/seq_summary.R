#' Sequence and taxon count summaries
#'
#' Returns the number of sequences and taxonomic features in the data. Useful for filtering
#'
#' @param file The name of the file which the feature data are to be read from. Each row of the table appears as one line of the file.
#'   If it does not contain an absolute path, the file name is relative to the current working directory, \code{getwd()}.
#'   Tilde-expansion is performed where supported.
#'
#'   This can be a compressed file (see \code{\link[base]{connections}}).
#' @param list Whether or not to return the structure of the HDF5 file.
#'    If TRUE, returns only the list of the file contents.
#'    If FALSE (the default), returns the actual BIOM data.
#'
#' @details QIIME 2 uses the term features to refer to microbial taxonomic units in a way that is agnostic to the user's decision to use
#'   ASVs or OTUs. OTUs must first be clustered in QIIME 2, commonly at 97\% sequence similarity. ASVs represent unique sequences and so these tables will be
#'   larger than OTU tables for the same data set. QIIME 2 favors ASVs by default. For a discussion on the benefits of using ASVs over
#'   OTUs see: \url{https://www.nature.com/articles/ismej2017119} (referenced below).
#'
#'
#' @return \code{read_qiime2_table} returns a long-form data frame of microbial features (ASVs or OTUs) and sample IDs. 
#'
#' @seealso \code{\link{read_qiime2_tax}}
#'
#' @examples
#' data(example_data)
#'
#' # initial sequence and ASV count? 
#' seq_check(dat)
#'
#' Identify lineages to remove
#' unassign <- tax[Kingdom == 'Unassigned', taxon_id]
#' euk <- tax[Kingdom == 'Eukaryota', taxon_id]
#' arch <- tax[Kingdom == 'Archaea', taxon_id]
#' mito_chloro <- tax[, tax_string := paste(Phylum, Class, Order, Family, Genus, Species, sep = ';')
#'                    ][grepl('mitochond|chloroplast', tax_string, ignore.case = T), taxon_id]
#' 
#' cat('Unassigned\n'); seq_check(dat[taxon_id %in% unassign]) 
#' cat('Eukarya\n'); seq_check(dat[taxon_id %in% euk])  
#' cat('Archaea\n');seq_check(dat[taxon_id %in% arch]) 
#' cat('Mitochondria, Chloroplasts\n');seq_check(dat[taxon_id %in% mito_chloro])  
#'
#' Filter out lineages
#' dat <- dat[!taxon_id %in% c(arch, mito_chloro, euk, unassign)]
#' cat('\n\nFinal after filtering\n'); seq_check(dat)
#'
#' @export

seq_summary <- function(x, reads = c(), tax_id = c()) {
  if(is.null(reads) || is.null(tax_id)) stop("Must supply columns for sequence reads and taxon ID")
  nums <- formatC(c(sum(x[[reads]], na.rm = TRUE), 
                    uniqueN(x[[tax_id]])), 
                  format = 'fg', width = 12, big.mark = ',')
  cat(nums[1], ' sequences\n', nums[2], ' ASVs', sep = ''); cat('\n')
}
