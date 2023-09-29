# Add documentation here

seq_summary <- function(x, reads = c(), tax_id = c()) {
  if(is.null(reads) || is.null(tax_id)) stop("Must supply columns for sequence reads and taxon ID")
  nums <- formatC(c(sum(x[[reads]], na.rm = TRUE), 
                    uniqueN(x[[tax_id]])), 
                  format = 'fg', width = 12, big.mark = ',')
  cat(nums[1], ' sequences\n', nums[2], ' ASVs', sep = ''); cat('\n')
}
