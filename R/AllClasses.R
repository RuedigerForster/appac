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
