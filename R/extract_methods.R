# Subsetting functions

setMethod("[[", "qsip", function(x, i, ...) {
  list(x)[[which(attributes(x)$names==i)]]
})

setMethod("$", "qsip", function(x, i, ...) {
  list(x)[[which(attributes(x)$names==i)]]
})
