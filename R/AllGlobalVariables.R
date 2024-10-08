
column_names <- c("Sample_Name", "Peak_Name", "Injection_Date",
                  "Air_Pressure", "Raw_Area")

default_palette <- list(
  highlight_color = "cadetblue2",
  lowlight_color = "honeydew4",
  line_color = "darkgoldenrod4",
  fill_color = "darkgoldenrod1"
)

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

x <- NULL
y <- NULL

