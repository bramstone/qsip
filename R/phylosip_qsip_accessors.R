# Accessor and subset methods for phylosip and qsip data

#' @rdname nrep-methods
#' @export
setGeneric('nrep', function(data) standardGeneric('nrep'))
#' @rdname nrep-methods
setMethod('nrep', 'phylosip', function(data) {
  length(unique(sample_data(data)[,data@qsip@rep_id]))
})

#' @rdname nmeasure-methods
#' @export
setGeneric('nmeasure', function(data) standardGeneric('nmeasure'))
#' @rdname nmeasure-methods
setMethod('nmeasure', 'phylosip', function(data) {
  length(data@qsip)
})
#' @rdname nmeasure-methods
setMethod('nmeasure', 'qsip', function(data) {
  length(data)
})
