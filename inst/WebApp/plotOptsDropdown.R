# Code for the 'Plot Options' dropdown in the 'Model Playground' Panel
plotOptsDropdown <- dropdown(
  inputId = "plot_dropdown",
  label = "Plot Options",
  style = "default",
  icon = icon("chart-area"),
  status = "default btn-controls",
  size = "lg",
  up = TRUE,
  #animate = TRUE,
  width = "100%",
  selectInput(
    inputId = "xaxis",
    label = "X-axis variable",
    choices = plot_options_params,
    selected = "time",
    multiple = FALSE,
    selectize = TRUE
  ),
  selectInput(
    inputId = "yaxis",
    label = "Y-axis variable",
    choices = plot_options_params[-1],
    selected = "Plv",
    multiple = TRUE,
    selectize = TRUE
  )
)

