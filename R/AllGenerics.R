## repaired a sever bug in expectedArea

setGeneric("Y", function(object) {
  standardGeneric("Y")
})

setGeneric("X", function(object) {
  standardGeneric("X")
})

setGeneric("u_Y", function(object) {
  standardGeneric("u_Y")
})

setGeneric("u_X", function(object) {
  standardGeneric("u_X")
})

setGeneric("weight", function(object) {
  standardGeneric("weight")
})

setGeneric("stdErrors", function(object, ...) {
  standardGeneric("stdErrors")
})

setGeneric("pValues", function(object, ...) {
  standardGeneric("pValues")
})

setGeneric("dates", function(object, sample) {
  standardGeneric("dates")
})

setGeneric("P", function(object, sample) {
  standardGeneric("P")
})

setGeneric("rawAreas", function(object, sample) {
  standardGeneric("rawAreas")
})

setGeneric("correctedAreas", function(object, sample) {
  standardGeneric("correctedAreas")
})

setGeneric("compensatedRawAreas", function(object, sample) {
  standardGeneric("compensatedRawAreas")
})

setGeneric("compensatedCorrectedAreas", function(object, sample) {
  standardGeneric("compensatedCorrectedAreas")
})

setGeneric("expectedAreas", function(object, sample) {
  standardGeneric("expectedAreas")
})

setGeneric("covMatrix", function(object, sample, type) {
  standardGeneric("covMatrix")
})

setGeneric("corMatrix", function(object, sample, type) {
  standardGeneric("corMatrix")
})

setGeneric("residualAreas", function(object, sample,
                                     type = c(
                                       "raw.area",
                                       "corrected.area",
                                       "compensated.corrected.area"
                                     )) {
  standardGeneric("residualAreas")
})

setGeneric("variance", function(object, sample, type) {
  standardGeneric("variance")
})

setGeneric("driftFactor", function(object, date, area) {
  standardGeneric("driftFactor")
})

setGeneric("plotGlobalFit", function(object, size = 4, ...)
  standardGeneric("plotGlobalFit"),
  signature = "object"
)

setGeneric("plotControlChart", function(object, sample, peak, size = 4, ...)
  standardGeneric("plotControlChart"),
  signature = c("object", "sample", "peak")
)

setGeneric("plotLocalFit", function(object, sample, peak, size = 4, ...)
  standardGeneric("plotLocalFit"),
  signature = c("object", "sample", "peak")
)
