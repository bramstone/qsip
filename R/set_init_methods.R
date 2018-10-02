#' @export

setGeneric('initialize', function(.Object, ...) {
  standardGeneric('initialize')})

#' The S4 method for initializing a qsip object
setMethod('initialize', 'qsip',
          function(.Object,
                   density=character(0),
                   abund=character(0),
                   rep_id=character(0),
                   rep_group=character(0),
                   iso=character(0),
                   timepoint=character(0),
                   filter=NULL) {
            .Object <- callNextMethod()
            # code here for any QC if-stop statements
            .Object
          })
# Maybe filter='data.frameOrNULL')

# The S4 method for initializing a phylosip object
# NOTE: doesn't seem to be needed, as you can change classes with setAs() not new()
#
# @import methods
# @importClassesFrom phyloseq phyloseq
# setMethod('initialize,', 'phylosip',
#          function(.Object='phyloseq', ...) {
#            .Object <- callNextMethod()
 #           .Object
  #        })
