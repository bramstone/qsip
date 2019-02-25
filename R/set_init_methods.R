#' @export

setGeneric('initialize', function(.Object, ...) {
  standardGeneric('initialize')})

#' The S4 method for initializing a qsip object
setMethod('initialize', 'qsip',
          function(.Object,
                   density=character(0),
                   abund=character(0),
                   rep_id=character(0),
                   rep_num=character(0),
                   rep_group=character(0),
                   iso=character(0),
                   iso_trt=character(0),
                   timepoint=character(0),
                   filter_levels=NULL,
                   filter=character(0)) {
            .Object <- callNextMethod()
            # code here for any QC if-stop statements
            .Object
          })
# Maybe filter='data.frameOrNULL')
