#------------- not exported methods -----------------------

setMethod("Y", signature(object = "Correction"), function(object) {
  res <- lapply(object@local.fits, function(x) x$slope)
  nam <- lapply(object@local.fits, function(x) x$peak)
  for (i in seq_along(res)) {
    names(res[[i]]) <- nam[[i]]
  }
  return(res)
})

setMethod("X", signature(object = "Correction"), function(object) {
  res <- lapply(object@local.fits, function(x) x$area.ref)
  nam <- lapply(object@local.fits, function(x) x$peak)
  for (i in seq_along(res)) {
    names(res[[i]]) <- nam[[i]]
  }
  return(res)
})

setMethod("u_Y", signature(object = "Correction"), function(object) {
  res <- lapply(object@local.fits, function(x) x$se.slope)
  nam <- lapply(object@local.fits, function(x) x$peak)
  for (i in seq_along(res)) {
    names(res[[i]]) <- nam[[i]]
  }
  return(res)
})

setMethod("u_X", signature(object = "Correction"), function(object) {
  res <- lapply(object@local.fits, function(x) x$se.area.ref)
  nam <- lapply(object@local.fits, function(x) x$peak)
  for (i in seq_along(res)) {
    names(res[[i]]) <- nam[[i]]
  }
  return(res)
})

setMethod("weight", signature(object = "Correction"), function(object) {
  res <- lapply(object@local.fits, function(x) x$n)
  nam <- lapply(object@local.fits, function(x) x$peak)
  for (i in seq_along(res)) {
    names(res[[i]]) <- nam[[i]]
  }
  return(res)
})

#------------- Correction methods -------------------------

setMethod(
  "coefficients",
  signature(object = "Correction"),
  function(object, ...) {
    return(object@global.fit$coefficients)
  }
)

setMethod(
  "stdErrors",
  signature(object = "Correction"),
  function(object, ...) {
    return(object@global.fit$standard.errors)
  }
)

setMethod(
  "pValues",
  signature(object = "Correction"),
  function(object, ...) {
    return(object@global.fit$p.values)
  }
)

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
    res <- object@samples[[sample]][[type]] -
      object@samples[[sample]]$expected.area
    return(res)
  }
)

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
    return(diag(covMatrix(object, sample, type)))
  }
)


#------------- Drift methods ------------------------------

setMethod(
  "driftFactor",
  signature(
    object = "Drift",
    date = "IDate",
    area = "matrix"
  ),
  function(object, date, area) {
    if (nrow(area) != length(date)) {
      stop("Inconsistent sizes of 'date' and 'area'")
    }
    date <- data.table::as.IDate(date)
    coefs <- lapply(
      date,
      function(d) {
        object@drift.factors[object@drift.factors$date == d, c(3, 4, 5)]
      }
    )
    coefs <- do.call(rbind, coefs)
    factors <- coefs[, 1] + coefs[, 2] * area + coefs[, 3] * area^2
    return(factors)
  }
)


#------------- Appac methods ------------------------------

setMethod(
  "plotGlobalFit",
  signature(
    object = "Appac"
  ),
  function(object, size = 4, ...) {
    if (!hasArg(colors) || (hasArg(colors) & (length(colors) != 4 || !all(names(colors) %in% c(
      "highlight_color",
      "lowlight_color",
      "line_color",
      "fill_color"
    ))))) {
      colors <- default_palette
      message("Using default color palette.")
    }
    data <- data.frame(
      x = unname(unlist(X(object@correction))),
      y = unname(unlist(Y(object@correction))),
      u_x = unname(unlist(u_X(object@correction))),
      u_y = unname(unlist(u_Y(object@correction))),
      n = unname(unlist(weight(object@correction)))
    )
    fit_coefs <- coefficients(object@correction)
    p <- .plot_global_fit(data, fit_coefs, colors, size)
    return(p)
  }
)

setMethod(
  "plotControlChart",
  signature(
    object = "Appac",
    sample = "character",
    peak = "character"
  ),
  function(object, sample, peak, size = 4, ...) {
    if (!sample %in% names(object@correction@samples)) {
      stop(
        "Unknown sample '", sample,
        "'. Sample may be any of: ",
        paste0(
          "'",
          names(object@correction@samples), sep = "' "
        )
      )
    }
    if (!(peak %in% colnames(object@correction@samples[[sample]]$raw.area))) {
      stop(
        "Sample '", sample,
        "' does not include peak '", peak, "'. Peak may be any of: ",
        paste0(
          "'",
          colnames(object@correction@samples[[sample]]$raw.area), sep = "' "
        )
      )
    }
    if (!hasArg(colors) || (hasArg(colors) & (length(colors) != 4 || !all(names(colors) %in% c(
      "highlight_color",
      "lowlight_color",
      "line_color",
      "fill_color"
    ))))) {
      colors <- default_palette
      message("Using default color palette.")
    }
    p <- .plot_control_chart(
      data = object, sample = sample,
      peak = peak, colors = colors, size = size,
      show.compensated.areas = TRUE,
      plot.residuals = TRUE, bins = 50
    )
    return(p)
  }
)

setMethod(
  "plotLocalFit",
  signature(
    object = "Appac",
    sample = "character",
    peak = "character"
  ),
  function(object, sample, peak, size = 4, ...) {
    if (!sample %in% names(object@correction@samples)) {
      stop(
        "Unknown sample '", sample,
        "'. Sample may be any of: ",
        paste0("'", names(object@correction@samples), sep = "' ")
      )
    }
    if (!(peak %in% colnames(object@correction@samples[[sample]]$raw.area))) {
      stop(
        "Sample '", sample, "' does not include peak '", peak,
        "'. Peak may be any of: ",
        paste0("'", colnames(object@correction@samples[[sample]]$raw.area),
               sep = "' ")
      )
    }
    if (!hasArg(colors) || (hasArg(colors) & (length(colors) != 4 || !all(names(colors) %in% c(
      "highlight_color",
      "lowlight_color",
      "line_color",
      "fill_color"
    ))))) {
      colors <- default_palette
      message("Using default color palette.")
    }
    p <- .plot_local_fit(
      data = object,
      sample = sample,
      peak = peak,
      colors = colors,
      size = size,
      show.compensated.areas = TRUE,
      plot.residuals = TRUE,
      bins = 50
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
