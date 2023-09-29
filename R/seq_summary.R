#' Sequence and taxon count summaries
#'
#' Returns the number of sequences and taxonomic features in the data. Useful for filtering.
#'
#' @param data data.table object in long-format where each row represents a sequence feature in a given sequenced sample.
#' @param reads Column name specifying sequence read data.
#'   Should be sequence read counts and not relativized abundances.
#' @param tax_id Column name specifying a unique identifier for each sequence feature.
#'
#' @details QIIME 2 uses the term features to refer to microbial taxonomic units in a way that is agnostic to the user's decision to use
#'   ASVs or OTUs. OTUs must first be clustered in QIIME 2, commonly at 97\% sequence similarity. ASVs represent unique sequences and so these tables will be
#'   larger than OTU tables for the same data set. QIIME 2 favors ASVs by default. For a discussion on the benefits of using ASVs over
#'   OTUs see: \url{https://www.nature.com/articles/ismej2017119} (referenced below).
#'
#'
#' @return Produces formatted count of total sequence reads and unique taxonomic features across the data.
#'
#' @seealso \code{\link{freq_filter}}
#'
#' @examples
#' data(example_data)
#'
#' # initial sequence and ASV count? 
#' seq_summary(dat, 'seq_abund', 'taxon_id')
#'
#' # Remove global singletons and doubletons
#' dat <- dat[seq_abund > 2]
#' seq_summary(dat, 'seq_abund', 'taxon_id')
#'
#' # Identify lineages to remove
#' unassign <- tax[Kingdom == 'Unassigned', taxon_id]
#' euk <- tax[Kingdom == 'Eukaryota', taxon_id]
#' arch <- tax[Kingdom == 'Archaea', taxon_id]
#' mito_chloro <- tax[, tax_string := paste(Phylum, Class, Order, Family, Genus, Species, sep = ';')
#'                    ][grepl('mitochond|chloroplast', tax_string, ignore.case = T), taxon_id]
#' 
#' cat('Unassigned\n'); seq_summary(dat, 'seq_abund', 'taxon_id')
#' cat('Eukarya\n'); seq_summary(dat, 'seq_abund', 'taxon_id')
#' cat('Archaea\n'); seq_summary(dat, 'seq_abund', 'taxon_id')
#' cat('Mitochondria, Chloroplasts\n'); seq_summary(dat, 'seq_abund', 'taxon_id')
#'
#' # Filter out lineages
#' dat <- dat[!taxon_id %in% c(arch, mito_chloro, euk, unassign)]
#' cat('\n\nFinal after filtering\n'); seq_summary(dat, 'seq_abund', 'taxon_id')
#'
#' @export

seq_summary <- function(data, reads = c(), tax_id = c()) {
  if(is.null(reads) || is.null(tax_id)) stop("Must supply columns for sequence reads and taxon ID")
  if(length(reads) > 1 || length(tax_id) > 1) stop("Must supply single column name")
  nums <- formatC(c(sum(x[[reads]], na.rm = TRUE), 
                    uniqueN(x[[tax_id]])), 
                  format = 'fg', width = 12, big.mark = ',')
  cat(nums[1], ' sequences\n', nums[2], ' Features', sep = ''); cat('\n')
}
