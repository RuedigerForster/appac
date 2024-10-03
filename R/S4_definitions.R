library(data.table)

source("./R/plot.R")

#------------- Class definitions --------------------------


setClass("Drift",
  slots = c(
    drift.factors = "data.frame",
    bias = "list",
    samples = "list"
  ),
  prototype = c(
    drift.factors = data.frame(),
    bias = list(),
    samples = list()
  )
)

setClass("Correction",
  slots = c(
    global.fit = "list",
    local.fits = "list",
    samples = "list"
  ),
  prototype = c(
    global.fit = list(),
    local.fits = list(),
    samples = list()
  )
)

setClass("Appac",
  slots = c(
    drift = "Drift",
    correction = "Correction"
  )
)

#------------- not exported methods -----------------------

setGeneric("Y", function(object) {
  standardGeneric("Y")
})

setMethod("Y", signature(object = "Correction"), function(object) {
  return(lapply(object@local.fits, function(x) x$slope))
})

setGeneric("X", function(object) {
  standardGeneric("X")
})

setMethod("X", signature(object = "Correction"), function(object) {
  return(lapply(object@local.fits, function(x) x$area.ref))
})

setGeneric("u_Y", function(object) {
  standardGeneric("u_Y")
})

setMethod("u_Y", signature(object = "Correction"), function(object) {
  return(lapply(object@local.fits, function(x) x$se.slope))
})

setGeneric("u_X", function(object) {
  standardGeneric("u_X")
})

setMethod(
  "u_X",
  signature(object = "Correction"), function(object) {
    return(lapply(object@local.fits, function(x) x$se.area.ref))
  }
)

setGeneric("weight", function(object) {
  standardGeneric("weight")
})

setMethod(
  "weight",
  signature(object = "Correction"),
  function(object) {
    return(lapply(object@local.fits, function(x) x$n))
  }
)

#------------- Correction methods -------------------------

setMethod("coefficients", signature(object = "Correction"), function(object, ...) {
  return(object@global.fit$coefficients)
})

setGeneric("stdErrors", function(object, ...) {
  standardGeneric("stdErrors")
})

setMethod(
  "stdErrors",
  signature(object = "Correction"),
  function(object, ...) {
    return(object@global.fit$standard.errors)
  }
)

setGeneric("pValues", function(object, ...) {
  standardGeneric("pValues")
})

setMethod(
  "pValues",
  signature(object = "Correction"),
  function(object, ...) {
    return(object@global.fit$p.values)
  }
)

setGeneric("dates", function(object, sample) {
  standardGeneric("dates")
})

setMethod(
  "dates",
  signature(
    sample = "character",
    object = "Correction"
  ),
  function(object, sample) {
    return(object@samples[[sample]]$date)
  }
)

setGeneric("P", function(object, sample) {
  standardGeneric("P")
})

setMethod(
  "P",
  signature(
    object = "Correction",
    sample = "character"
  ),
  function(object, sample) {
    return(object@samples[[sample]]$pressure)
  }
)

setGeneric("rawAreas", function(object, sample) {
  standardGeneric("rawAreas")
})

setMethod(
  "rawAreas",
  signature(
    object = "Correction",
    sample = "character"
  ),
  function(object, sample) {
    return(object@samples[[sample]]$raw.area)
  }
)

setGeneric("correctedAreas", function(object, sample) {
  standardGeneric("correctedAreas")
})

setMethod(
  "correctedAreas",
  signature(
    object = "Correction",
    sample = "character"
  ),
  function(object, sample) {
    return(object@samples[[sample]]$corrected.area)
  }
)

setGeneric("compensatedRawAreas", function(object, sample) {
  standardGeneric("compensatedRawAreas")
})

setMethod(
  "compensatedRawAreas",
  signature(
    object = "Drift",
    sample = "character"
  ),
  function(object, sample) {
    return(
      object@samples[[sample]]$compensated.raw.area
    )
  }
)

setGeneric("compensatedCorrectedAreas", function(object, sample) {
  standardGeneric("compensatedCorrectedAreas")
})

setMethod(
  "compensatedCorrectedAreas",
  signature(
    object = "Correction",
    sample = "character"
  ),
  function(object, sample) {
    return(
      object@samples[[sample]]$compensated.corrected.area
    )
  }
)

setGeneric("expectedAreas", function(object, sample) {
  standardGeneric("expectedAreas")
})

setMethod(
  "expectedAreas",
  signature(
    object = "Correction",
    sample = "character"
  ),
  function(object, sample) {
    return(
      object@samples[[sample]]$expected.area
    )
  }
)

setGeneric("covMatrix", function(object, sample, type) {
  standardGeneric("covMatrix")
})

setMethod(
  "covMatrix",
  signature(
    object = "Correction",
    sample = "character",
    type = "character"
  ),
  function(object,
           sample,
           type = c(
             "raw.area",
             "compensated.corrected.area"
           )) {
    type <- match.arg(type)
    return(cov(object@samples[[sample]][[type]]))
  }
)

setGeneric("corMatrix", function(object, sample, type) {
  standardGeneric("corMatrix")
})

setMethod(
  "corMatrix",
  signature(
    object = "Correction",
    sample = "character",
    type = "character"
  ),
  function(object,
           sample,
           type =
             c(
               "raw.area",
               "compensated.corrected.area"
             )) {
    type <- match.arg(type)
    return(cor(object@samples[[sample]][[type]]))
  }
)

setGeneric("residualAreas", function(object, sample,
                                     type = c(
                                       "raw.area",
                                       "corrected.area",
                                       "compensated.corrected.area"
                                     )) {
  standardGeneric("residualAreas")
})

setMethod(
  "residualAreas",
  signature(
    object = "Correction",
    sample = "character",
    type = "character"
  ),
  function(object,
           sample,
           type = c(
             "raw.area",
             "corrected.area",
             "compensated.corrected.area"
           )) {
    type <- match.arg(type)
    return(object@samples[[sample]][[type]] - object@samples[[sample]]$expected.area)
  }
)

setGeneric("variance", function(object, sample, type) {
  standardGeneric("variance")
})

setMethod(
  "variance",
  signature(
    object = "Correction",
    sample = "character",
    type = "character"
  ),
  function(object,
           sample,
           type = c(
             "raw.area",
             "compensated.corrected.area"
           )) {
    type <- match.arg(type)
    return(diag(covMatrix(object, sample, type))) # %/% cor.mat( object, sample,type)))
  }
)


#------------- Drift methods ------------------------------

setGeneric("driftFactor", function(object, date, area) {
  standardGeneric("driftFactor")
})

setMethod(
  "driftFactor",
  signature(
    object = "Drift",
    date = "IDate",
    area = "matrix"
  ),
  function(object, date, area) {
    if (nrow(area) != length(date)) stop("Inconsistent sizes of 'date' and 'area'")
    date <- data.table::as.IDate(date)
    coefs <- lapply(date, function(d) object@drift.factors[object@drift.factors$date == d, c(3, 4, 5)])
    coefs <- do.call(rbind, coefs)
    factor <- coefs[, 1] + coefs[, 2] * area + coefs[, 3] * area^2
    return(factor)
  }
)


#------------- Appac methods ------------------------------

setGeneric("plotGlobalFit", function(object, size, colors) {
  standardGeneric("plotGlobalFit")
})

setMethod(
  "plotGlobalFit",
  signature(
    object = "Appac",
    size = "numeric",
    colors = "list"
  ),
  function(object, size, colors) {
    if (length(colors) != 4 || !all(names(colors) %in% c(
      "highlight_color",
      "lowlight_color",
      "line_color",
      "fill_color"
    ))) {
      stop("Incorrect colors.")
    }
    data <- data.frame(
      x = unname(unlist(X(object@correction))),
      y = unname(unlist(Y(object@correction))),
      u_x = unname(unlist(u_X(object@correction))),
      u_y = unname(unlist(u_Y(object@correction))),
      n = unname(unlist(weight(object@correction)))
    )
    fit.coefs <- coefficients(object@correction)
    p <- plot_global_fit(data, fit.coefs, colors)
    return(p)
  }
)

setGeneric("plotControlChart", function(object, sample, peak, size, colors) {
  standardGeneric("plotControlChart")
})

setMethod(
  "plotControlChart",
  signature(
    object = "Appac",
    sample = "character",
    peak = "character",
    size = "numeric",
    colors = "list"
  ),
  function(object, sample, peak, size, colors) {
    if (!sample %in% names(object@correction@samples)) {
      stop("Unknown sample '", sample, "'. Sample may be any of: ", paste0("'", names(object@correction@samples), sep = "' "))
    }
    if (!(peak %in% colnames(object@correction@samples[[sample]]$raw.area))) {
      # cat("I was here _o_")
      stop(
        "Sample '", sample, "' does not include peak '", peak, "'. Peak may be any of: ",
        paste0("'", colnames(object@correction@samples[[sample]]$raw.area), sep = "' ")
      )
    }
    if (length(colors) != 4 || !all(names(colors) %in% c(
      "highlight_color",
      "lowlight_color",
      "line_color",
      "fill_color"
    ))) {
      stop("Incorrect colors.")
    }
    p <- plot_control_chart(
      data = object, sample = sample,
      peak = peak, colors = colors, size = size,
      show.compensated.areas = TRUE,
      plot.residuals = TRUE, bins = 50
    )
    return(p)
  }
)

setGeneric("plotLocalFit", function(object, sample, peak, size, colors) {
  standardGeneric("plotLocalFit")
})

setMethod(
  "plotLocalFit",
  signature(
    object = "Appac",
    sample = "character",
    peak = "character",
    size = "numeric",
    colors = "list"
  ),
  function(object, sample, peak, size, colors) {
    if (!sample %in% names(object@correction@samples)) {
      stop("Unknown sample '", sample, "'. Sample may be any of: ", paste0("'", names(object@correction@samples), sep = "' "))
    }
    if (!(peak %in% colnames(object@correction@samples[[sample]]$raw.area))) {
      # cat("I was here _o_")
      stop(
        "Sample '", sample, "' does not include peak '", peak, "'. Peak may be any of: ",
        paste0("'", colnames(object@correction@samples[[sample]]$raw.area), sep = "' ")
      )
    }
    if (length(colors) != 4 || !all(names(colors) %in% c(
      "highlight_color",
      "lowlight_color",
      "line_color",
      "fill_color"
    ))) {
      stop("Incorrect colors.")
    }
    p <- plot_local_fit(
      data = object, sample = sample,
      peak = peak, coefs = coefficients(object@correction),
      colors = colors, size = size, show.compensated.areas = TRUE,
      plot.residuals = TRUE, bins = 50
    )
    return(p)
  }
)


#------------- TODO methods ------------------------------

setMethod(
  "summary",
  signature(object = "Appac"),
  function(object) {
    return(0)
  }
)





# setGeneric("aref", function(sample, peak, object) {
#   standardGeneric("aref")
# })
#
# #
# # Not yet converted to object Correction ???
# #
# setMethod("aref",
#   signature(sample = "character",
#             peak = "character",
#             object = "Appac"),
#   function(sample, peak, object) {
#   idx <- object@correction@samples[[sample]]$local.fits$peak == peak
#   # first <- object@correction@samples[[sample]]$local.fits$begin[idx]
#   # last <- object@correction@samples[[sample]]$local.fits$end[idx]
#   # area.ref <-
#   object@correction@samples[[sample]]$local.fits$area.ref[idx]
#   # dates <- seq(first[1], last[length(last)], 1)
#   # aref <-
#   # return(object@correction@samples[[sample]]$residuals[, peak])
# })
