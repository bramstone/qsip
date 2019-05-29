#' @rdname show-methods
setMethod("show", "qsip", function(object){
  # print qsip list (always there).
  cat(paste("qSIP data:          [ list of ",
            length(object), " qSIP measures ]", sep = ""), fill = TRUE)
  if( length(object@rep_group) == 1 ){
    cat("    replicate groups are specified", fill=TRUE)
  } else if( length(object@rep_num) == 1 ){
    cat("    replicate matching is specified", fill=TRUE)
  } else if( length(object@timepoint) == 1 ){
    cat("    different timepoints are specified", fill=TRUE)
  }
  show(as(object, "list"))
})

#' method extensions to show for phylosip objects.
#' Most code is taken directly from the phyloseq package
#'
#' See the general documentation of \code{\link[methods]{show}} method for
#' expected behavior.
#'
#' @seealso \code{\link[methods]{show}}
#'
#' @inheritParams methods::show
#' @export
#' @rdname show-methods
setMethod("show", "phylosip", function(object){
  cat("phylosip-class experiment-level object", fill=TRUE)
  # print qsip data (always there).
  cat(paste("specify_qsip() qSIP Data:         [ ", nrep(object), " true replicates and ",
            nmeasure(object), " qSIP measures ]", sep = ""), fill = TRUE)

  # print otu_table (always there).
  cat(paste("otu_table()    OTU Table:         [ ", ntaxa(otu_table(object)), " taxa and ",
            nsamples(otu_table(object)), " samples ]", sep = ""), fill = TRUE)

  # print Sample Data if there
  if(!is.null(sample_data(object, FALSE))){
    cat(paste("sample_data()  Sample Data:       [ ", dim(sample_data(object))[1], " samples by ",
              dim(sample_data(object))[2],
              " sample variables ]", sep = ""), fill = TRUE)
  }

  # print tax Tab if there
  if(!is.null(tax_table(object, FALSE))){
    cat(paste("tax_table()    Taxonomy Table:    [ ", dim(tax_table(object))[1], " taxa by ",
              dim(tax_table(object))[2],
              " taxonomic ranks ]", sep = ""), fill = TRUE)
  }

  # print tree if there
  if(!is.null(phy_tree(object, FALSE))){
    cat(paste("phy_tree()    Phylogenetic Tree: [ ", ntaxa(phy_tree(object)), " tips and ",
              phy_tree(object)$Nnode,
              " internal nodes ]", sep = ""),
        fill = TRUE
    )
  }

  # print refseq summary if there
  if(!is.null(refseq(object, FALSE))){
    cat(paste("refseq()      ", class(refseq(object))[1], ":      [ ", ntaxa(refseq(object)), " reference sequences ]", sep = ""), fill=TRUE)
  }

})
