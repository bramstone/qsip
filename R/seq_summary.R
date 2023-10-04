#' Sequence and taxon count summaries
#'
#' Returns the number of sequences and taxonomic features in the data. Useful for filtering.
#'
#' @param data data.table object in long-format where each row represents a sequence
#'  feature in a given sequenced sample.
#'  Alternatively, a data.table containing only two columns corresponding to sequence
#'   reads and taxonomic ID may be supplied (see \code{details}).
#' @param reads Column name specifying sequence read data.
#'   Should be sequence read counts and not relativized abundances.
#' @param tax_id Column name specifying a unique identifier for each sequence feature.
#'
#' @details QIIME 2 uses the term features to refer to microbial taxonomic units in
#'   a way that is agnostic to the user's decision to use ASVs or OTUs.
#'   OTUs must first be clustered in QIIME 2, commonly at 97\% sequence similarity.
#'   ASVs represent unique sequences and so these tables will be larger than OTU tables for the same data set.
#'   QIIME 2 favors ASVs by default.
#'   For a discussion on the benefits of using ASVs over OTUs see:
#'   \url{https://www.nature.com/articles/ismej2017119} (referenced below).
#'
#'   If supplying a data.table containing only two columns, \code{seq_summary} will
#'   attempt to identify which column specifies the taxonomic IDs and which specifies
#'   the sequence read count. In this case, it is not necessary to specify these
#'   columns in the other arguments.
#'
#'
#' @return Produces formatted count of total sequence reads and unique taxonomic features across the data.
#'
#' @seealso \code{\link{freq_filter}}
#'
#' @examples
#' data(example_qsip)
#'
#' # initial sequence and ASV count?
#' seq_summary(example_qsip, 'seq_abund', 'asv_id')
#' seq_summary(example_qsip[, c('seq_abund', 'asv_id')]
#' seq_summary(example_qsip[, c(1, 12)]) # should be the same as above
#'
#' # Identify lineages to remove
#' tx <- unique(example_qsip[, c('asv_id', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')])
#' unassign <- tx[Kingdom == 'Unassigned', asv_id]
#' euk <- tx[Kingdom == 'Eukaryota', asv_id]
#' arch <- tx[Kingdom == 'Archaea', asv_id]
#' mito_chloro <- tx[, tax_string := paste(Phylum, Class, Order, Family, Genus, sep = ';')
#'                   ][grepl('mitochond|chloroplast', tax_string, ignore.case = T), asv_id]
#'
#' cat('Unassigned\n'); seq_summary(example_qsip[asv_id %in% unassign, c(1, 12)])
#' cat('Eukarya\n'); seq_summary(example_qsip[asv_id %in% euk, c(1, 12)])
#' cat('Archaea\n'); seq_summary(example_qsip[asv_id %in% arch, c(1, 12)])
#' cat('Mitochondria, Chloroplasts\n'); seq_summary(example_qsip[asv_id %in% mito_chloro, c(1, 12)])
#'
#' # Filter out lineages
#' example_qsip <- example_qsip[!asv_id %in% c(arch, mito_chloro, euk, unassign)]
#' cat('\n\nFinal after filtering\n'); seq_summary(example_qsip, 'seq_abund', 'asv_id')
#' rm(tx, unassign, euk, arch, mito_chloro)
#'
#' @export

seq_summary <- function(data, reads = c(), tax_id = c()) {
  if(length(data) == 2) {
    tax_id <- which(sapply(data, is.character))
    reads <- which(sapply(data, function(x) is.numeric(x) | is.integer(x)))
    if(length(tax_id) != 1 || length(reads) != 1) stop("A single numeric and character column could not be identified")
  } else {
    if(is.null(reads) || is.null(tax_id)) stop("Must supply columns for sequence reads and taxon ID")
    if(length(reads) > 1 || length(tax_id) > 1) stop("Must supply single column name")
  }
  nums <- formatC(c(sum(data[[reads]], na.rm = TRUE),
                    uniqueN(data[[tax_id]])),
                  format = 'fg', width = 12, big.mark = ',')
  cat(nums[1], ' sequences\n', nums[2], ' tax features', sep = ''); cat('\n')
}
