% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seq_summary.R
\name{seq_summary}
\alias{seq_summary}
\title{Sequence and taxon count summaries}
\usage{
seq_summary(data, reads = c(), tax_id = c())
}
\arguments{
\item{data}{data.table object in long-format where each row represents a sequence
feature in a given sequenced sample.
Alternatively, a data.table containing only two columns corresponding to sequence
 reads and taxonomic ID may be supplied (see \code{details}).}

\item{reads}{Column name specifying sequence read data.
Should be sequence read counts and not relativized abundances.}

\item{tax_id}{Column name specifying a unique identifier for each sequence feature.}
}
\value{
Produces formatted count of total sequence reads and unique taxonomic features across the data.
}
\description{
Returns the number of sequences and taxonomic features in the data. Useful for filtering.
}
\details{
QIIME 2 uses the term features to refer to microbial taxonomic units in
  a way that is agnostic to the user's decision to use ASVs or OTUs.
  OTUs must first be clustered in QIIME 2, commonly at 97\% sequence similarity.
  ASVs represent unique sequences and so these tables will be larger than OTU tables for the same data set.
  QIIME 2 favors ASVs by default.
  For a discussion on the benefits of using ASVs over OTUs see:
  \url{https://www.nature.com/articles/ismej2017119} (referenced below).

  If supplying a data.table containing only two columns, \code{seq_summary} will
  attempt to identify which column specifies the taxonomic IDs and which specifies
  the sequence read count. In this case, it is not necessary to specify these
  columns in the other arguments.
}
\examples{
data(example_qsip)

# initial sequence and ASV count?
seq_summary(example_qsip, 'seq_abund', 'asv_id')
seq_summary(example_qsip[, c('seq_abund', 'asv_id')]
seq_summary(example_qsip[, c(1, 12)]) # should be the same as above

# Identify lineages to remove
tx <- unique(example_qsip[, c('asv_id', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')])
unassign <- tx[Kingdom == 'Unassigned', asv_id]
euk <- tx[Kingdom == 'Eukaryota', asv_id]
arch <- tx[Kingdom == 'Archaea', asv_id]
mito_chloro <- tx[, tax_string := paste(Phylum, Class, Order, Family, Genus, sep = ';')
                  ][grepl('mitochond|chloroplast', tax_string, ignore.case = T), asv_id]

cat('Unassigned\n'); seq_summary(example_qsip[asv_id \%in\% unassign, c(1, 12)])
cat('Eukarya\n'); seq_summary(example_qsip[asv_id \%in\% euk, c(1, 12)])
cat('Archaea\n'); seq_summary(example_qsip[asv_id \%in\% arch, c(1, 12)])
cat('Mitochondria, Chloroplasts\n'); seq_summary(example_qsip[asv_id \%in\% mito_chloro, c(1, 12)])

# Filter out lineages
example_qsip <- example_qsip[!asv_id \%in\% c(arch, mito_chloro, euk, unassign)]
cat('\n\nFinal after filtering\n'); seq_summary(example_qsip, 'seq_abund', 'asv_id')
rm(tx, unassign, euk, arch, mito_chloro)

}
\seealso{
\code{\link{freq_filter}}
}
