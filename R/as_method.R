#' @exportMethod
setAs('phyloseq', 'phylosip', function(from, to) {
  to <- new(to,
            otu_table=from@otu_table,
            tax_table=from@tax_table,
            sam_data=from@sam_data,
            phy_tree=from@phy_tree,
            refseq=from@refseq,
            qsip=new('qsip'))
})
# Conversion method to convert phyloseq to phylosip object
