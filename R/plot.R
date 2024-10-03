library(ggplot2)
library(ggpubr)
library(grDevices)
library(patchwork)
library(scales)

# mm <- 4
pts <- function(mm) {
  mm * ggplot2::.pt
}

main_plot_theme <- function(size = 4) {
  theme(
    text = element_text(size = rel(size)),
    plot.title = element_text(size = rel(4)),
    plot.subtitle = element_text(size = rel(3.5)),
    axis.title.x = element_text(size = rel(3.5)),
    axis.title.y = element_text(size = rel(3.5)),
    axis.text.x = element_text(size = rel(3.7)),
    axis.text.y = element_text(size = rel(3.7)),
    legend.position = "none"
  )
}

residuals_plot_theme <- function(size = 4) {
  theme(
    text = element_text(size = rel(size)),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title.x = element_text(size = rel(3.5)),
    axis.title.y = element_text(size = rel(3.5)),
    axis.text.x = element_text(size = rel(3.7)),
    axis.text.y = element_text(size = rel(3.7)),
    legend.position = "none"
  )
}

histogram_plot_theme <- function(size = 4) {
  theme(
    text = element_text(size = rel(size)),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title.x = element_text(size = rel(3.5)),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = rel(3.7)),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0)
  )
}


plot_global_fit <- function(data = NA, fit.coefs = NA,
                            colors = list(
                              highlight_color = NA,
                              lowlight_color  = NA,
                              line_color      = NA,
                              fill_color      = NA
                            ),
                            size = 4) {
  # for devtools::check()s sake:
  # alt: changing x -> X$x or y -> X$y in ggplot function will throw a warning
  x <- NULL
  y <- NULL
  if (!all(c("x", "y", "u_x", "u_y", "n") %in% colnames(data))) {
    stop("Must provide a data frame with columns: 'x', 'y', 'u_x', 'u_y', 'n'")
  }
  if (any(is.na(colors))) {
    stop("Colors missing in arguments")
  }
  if (is.na(size)) size <- 4
  dots_color <- unname(unlist(match.arg(colors$dots_color, colors)))
  errorbar_color <- unname(unlist(match.arg(colors$errorbar_color, colors)))
  line_color <- unname(unlist(match.arg(colors$line_color, colors)))
  fill_color <- unname(unlist(match.arg(colors$fill_color, colors)))
  title <- "Atmospheric Pressure Peak Area Correction (APPAC) function"
  subtitle <- "Global Fit"
  xlab <- "peak area at reference pressure"
  ylab <- "slope of peak area vs. pressure"
  min <- function(x, ux) {
    ifelse(x >= 0, x - abs(ux), x + abs(ux))
  }
  max <- function(x, ux) {
    ifelse(x >= 0, x + abs(ux), x - abs(ux))
  }
  xmin <- min(data$x, data$u_x) # plot with k=2
  xmax <- max(data$x, data$u_x)
  ymin <- min(data$y, data$u_y)
  ymax <- max(data$y, data$u_y)
  X <- data.frame(x = data$x, y = data$y, xmin, xmax, ymin, ymax, n = data$n)
  if (is.na(fit.coefs[2])) fit.coefs[2] <- 0
  draw_function <- function(x) {
    return(fit.coefs[1] * x + fit.coefs[2] * x * x)
  }
  plot <- ggplot(data = X, aes(x = x, y = y)) +
    stat_smooth(
      method = "lm",
      data = X,
      formula = y ~ poly(x, 1),
      level = 0.99, color = "NA", fill = colors$fill_color, alpha = 0.5
    ) +
    geom_function(fun = draw_function, linewidth = 1, color = colors$line_color) +
    geom_point(size = 3, shape = 20, color = colors$highlight_color) +
    geom_errorbar(X,
      mapping = aes(x = x, ymin = ymin, ymax = ymax),
      linewidth = 0.7, color = colors$lowlight_color
    ) +
    geom_errorbarh(X,
      mapping = aes(y = y, xmin = xmin, xmax = xmax),
      linewidth = 0.7, color = colors$lowlight_color
    ) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    main_plot_theme(size = size)

  return(plot)
}


plot_local_fit <- function(data, sample, peak, coefs, colors, size = 4, show.compensated.areas = TRUE, plot.residuals = TRUE, bins = 50) {
  if (any(is.na(colors))) {
    stop("Colors missing in arguments")
  }
  if (is.na(size)) size <- 4
  lowlight_color <- unname(unlist(match.arg(colors$lowlight_color, colors)))
  highlight_color <- unname(unlist(match.arg(colors$highlight_color, colors)))
  line_color <- unlist(match.arg(colors$line_color, colors))
  fill_color <- unlist(match.arg(colors$fill_color, colors))
  xlab <- "atmospheric pressure (hPa)"
  ylab <- "peak area"
  title <- "Local Fit"
  subtitle <- paste(sample, "-", peak)
  # plot data
  x <- P(object = data@correction, sample = sample)
  y1 <- rawAreas(object = data@correction, sample = sample)[, peak]
  y2 <- compensatedRawAreas(object = data@drift, sample = sample)[, peak]
  # plot margins
  area.ref <- data@correction@local.fits[[sample]]$area.ref[data@correction@local.fits[[sample]]$peak == peak]
  span <- 0.025 * area.ref
  l.y.limit <- area.ref - span
  h.y.limit <- area.ref + span
  # outliers
  outliers <- NULL
  if (!is.null(outliers)) {
    subset <- ifelse(outliers, 2, 1)
  } else {
    subset <- rep(1, length(dates))
  }
  X <- data.frame(x = x, y1 = y1, y2 = y2, s = subset)
  aref <- data@correction@local.fits[[sample]][data@correction@local.fits[[sample]]$peak == peak, "area.ref"]
  slop <- data@correction@local.fits[[sample]][data@correction@local.fits[[sample]]$peak == peak, "slope"]
  P_ref <- data@correction@global.fit$P.ref
  draw_function <- function(x) {
    return(aref + slop * (x - P_ref))
  }
  p1 <- ggplot(X) +
    # raw areas
    geom_point(
      data = na.omit(X), mapping = aes(x = x, y = y1),
      color = lowlight_color,
      size = 1,
      alpha = 1,
      shape = 20
    )
  # drift + bias compensated areas
  if (show.compensated.areas) {
    p1 <- p1 + geom_point(
      data = X, mapping = aes(x = x, y = y2),
      color = highlight_color, # factor(s)),
      size = 1,
      shape = 20,
      alpha = 0.5
    )
  }
  p1 <- p1 +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    ylim(l.y.limit, h.y.limit) +
    geom_function(fun = draw_function, linewidth = 1, color = line_color) +
    main_plot_theme(size = size)
  # residuals
  ylab <- "residuals"
  # plot data
  x <- P(object = data@correction, sample = sample)
  y1 <- residualAreas(object = data@correction, sample = sample, type = "raw")[, peak]
  y2 <- compensatedRawAreas(object = data@drift, sample = sample)[, peak] - expectedAreas(object = data@correction, sample = sample)[, peak]
  l.y.limit <- -span / 2.5
  h.y.limit <- span / 2.5
  # outliers
  outliers <- NULL
  if (!is.null(outliers)) {
    subset <- ifelse(outliers, 2, 1)
  } else {
    subset <- rep(1, length(dates))
  }
  Z <- as.data.frame(cbind(x = x, y1 = y1, y2 = y2, s = subset))
  binwidth <- (max(Z[, 2], na.rm = T) - min(Z[, 2], na.rm = T)) / bins
  p2 <- ggplot(Z) +
    # raw areas
    geom_point(
      data = Z, mapping = aes(x = x, y = y1),
      color = lowlight_color,
      size = 1,
      alpha = 0.5,
      shape = 20
    )
  # drift + bias compensated areas
  if (show.compensated.areas) {
    p2 <- p2 + geom_point(
      data = Z, mapping = aes(x = x, y = y2),
      color = highlight_color, # factor(s)),
      size = 1,
      shape = 20,
      alpha = 0.5
    )
  }
  p2 <- p2 +
    geom_hline(aes(yintercept = 0), color = line_color, linetype = "dashed", linewidth = 1) +
    labs(title = element_blank(), subtitle = element_blank(), x = xlab, y = ylab) +
    ylim(l.y.limit, h.y.limit) +
    residuals_plot_theme(size = size)
  # histogram of residuals
  p3 <- ggplot(data = na.omit(Z), mapping = aes(x = y1)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = binwidth, na.rm = T, position = "identity", fill = lowlight_color, color = lowlight_color)
  if (show.compensated.areas) {
    p3 <- p3 +
      geom_histogram(aes(x = y2, y = after_stat(density)), binwidth = binwidth, na.rm = T, position = "identity", fill = highlight_color, alpha = 0.4)
  }
  p3 <- p3 +
    geom_vline(aes(xintercept = 0), color = line_color, linetype = "dashed", linewidth = 1) +
    geom_function(color = line_color, linewidth = 1, fun = dnorm, args = list(mean = mean(y1, na.rm = T), sd = sd(y1, na.rm = T))) +
    histogram_plot_theme(size = size) +
    xlim(l.y.limit, h.y.limit) +
    coord_flip()

  p <- ggpubr::ggarrange(p1, p2 + p3 + patchwork::plot_layout(widths = c(3, 1)),
    ncol = 1, nrow = 2, heights = c(6, 3)
  )
  return(p)
}

plot_control_chart <- function(data, sample, peak, colors, size = 4, show.compensated.areas = TRUE, plot.residuals = TRUE, bins = 50) {
  if (any(is.na(colors))) {
    stop("Colors missing in arguments")
  }
  if (is.na(size)) size <- 4
  # plot settings
  lowlight_color <- unname(unlist(match.arg(colors$lowlight_color, colors)))
  highlight_color <- unname(unlist(match.arg(colors$highlight_color, colors)))
  line_color <- unname(unlist(match.arg(colors$line_color, colors)))
  xlab <- "date"
  ylab <- "peak area"
  title <- "Control Chart"
  subtitle <- paste(sample, "-", peak)
  # plot data
  x <- dates(object = data@correction, sample = sample)
  y1 <- rawAreas(object = data@correction, sample = sample)[, peak]
  y2 <- compensatedCorrectedAreas(object = data@correction, sample = sample)[, peak]
  # plot margins
  area.ref <- data@correction@local.fits[[sample]]$area.ref[data@correction@local.fits[[sample]]$peak == peak]
  span <- 0.025 * area.ref
  l.y.limit <- area.ref - span
  h.y.limit <- area.ref + span
  brk <- c(1, 3, 6, 12)
  nb <- as.numeric((max(x) - min(x)) / 300)
  y.breaks <- brk[which(abs(brk - nb) == min(abs(brk - nb)))]
  # outliers
  outliers <- NULL
  if (!is.null(outliers)) {
    subset <- ifelse(outliers, 2, 1)
  } else {
    subset <- rep(1, length(dates))
  }
  # plot data
  X <- data.frame(x = x, y1 = y1, y2 = y2, s = subset)
  X[, 1] <- as.Date(X[, 1], origin = "1970-01-01")
  std.dev <- sd(X[, 3])
  # control chart
  p1 <- ggplot(data = na.omit(X)) +
    geom_point(
      mapping = aes(x = x, y = y1),
      color = lowlight_color, size = 2,
      alpha = 1,
      shape = 20
    ) +
    scale_x_date(date_breaks = paste0(y.breaks, " months"), date_labels = "%m/%y")
  if (show.compensated.areas) {
    p1 <- p1 + geom_point(
      data = X, mapping = aes(x = x, y = y2),
      color = highlight_color,
      size = 2,
      shape = 20,
      alpha = 1
    )
  }
  p1 <- p1 +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    ylim(l.y.limit, h.y.limit) +
    geom_hline(yintercept = area.ref, linetype = "solid", color = line_color, linewidth = 1) +
    geom_hline(yintercept = area.ref + std.dev, linetype = "dotted", color = line_color, linewidth = 1) +
    geom_hline(yintercept = area.ref - std.dev, linetype = "dotted", color = line_color, linewidth = 1) +
    geom_hline(yintercept = area.ref + 2 * std.dev, linetype = "dashed", color = line_color, linewidth = 1) +
    geom_hline(yintercept = area.ref - 2 * std.dev, linetype = "dashed", color = line_color, linewidth = 1) +
    main_plot_theme(size = size)
  # sapply(bpts, function(a) geom_vline(xintercept=a, linetype = "dotted", color = line_color, linewidth = 1))
  #
  # residuals
  ylab <- "residuals"
  # plot data
  x <- dates(object = data@correction, sample = sample)
  y1 <- residualAreas(object = data@correction, sample = sample, type = "compensated.corrected.area")[, peak]
  y2 <- compensatedCorrectedAreas(object = data@correction, sample = sample)[, peak] - area.ref
  l.y.limit <- -span / 2.5
  h.y.limit <- span / 2.5
  # outliers
  outliers <- NULL
  if (!is.null(outliers)) {
    subset <- ifelse(outliers, 2, 1)
  } else {
    subset <- rep(1, length(x))
  }
  Z <- as.data.frame(cbind(x = x, y1 = y1, y2 = y2, s = subset))
  Z[, 1] <- as.Date(Z[, 1], origin = "1970-01-01") # as.Date(Z[, 1])
  binwidth <- (max(Z[, 2], na.rm = T) - min(Z[, 2], na.rm = T)) / bins
  p2 <- ggplot(Z) +
    # raw areas
    geom_point(
      data = Z, mapping = aes(x = x, y = y1),
      color = lowlight_color,
      size = 1,
      alpha = 0.5,
      shape = 20
    ) +
    labs(title = element_blank(), subtitle = element_blank(), x = xlab, y = ylab) +
    ylim(l.y.limit, h.y.limit) +
    scale_x_date(date_breaks = paste0(y.breaks, " months"), , date_labels = "%m/%y")
  # drift + bias compensated areas
  if (show.compensated.areas) {
    p2 <- p2 + geom_point(
      data = Z, mapping = aes(x = x, y = y2),
      color = highlight_color,
      size = 1,
      shape = 20,
      alpha = 0.5
    )
  }
  p2 <- p2 + geom_hline(aes(yintercept = 0), color = line_color, linetype = "dashed", linewidth = 1) +
    residuals_plot_theme(size = size)
  # histogram of residuals
  p3 <- ggplot(data = na.omit(Z), mapping = aes(x = y1)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = binwidth, na.rm = T, position = "identity", fill = lowlight_color, color = lowlight_color)
  if (show.compensated.areas) {
    p3 <- p3 +
      geom_histogram(aes(x = y2, y = after_stat(density)), binwidth = binwidth, na.rm = T, position = "identity", fill = highlight_color, alpha = 0.4)
  }
  p3 <- p3 +
    geom_vline(aes(xintercept = 0), color = line_color, linetype = "dashed", linewidth = 1) +
    geom_function(color = line_color, linewidth = 1, fun = dnorm, args = list(mean = mean(y1, na.rm = T), sd = sd(y1, na.rm = T))) +
    histogram_plot_theme(size = size) +
    xlim(l.y.limit, h.y.limit) +
    coord_flip()

  p <- ggpubr::ggarrange(p1, p2 + p3 + patchwork::plot_layout(widths = c(3, 1)),
    ncol = 1, nrow = 2, heights = c(6, 3)
  )
  return(p)
}
