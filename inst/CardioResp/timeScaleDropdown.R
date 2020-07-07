# Code for the 'Time Scale' dropdown in the 'Model Playground' Panel
timeScaleDropdown <- dropdown(
  inputId = "timescale_dropdown",
  label = "Time Scale",
  icon = HTML("<i class=\"far fa-hourglass\">", "</i>"),
  style = "default",
  status = "default btn-controls",
  size = "lg",
  up = TRUE,
  #animate = TRUE,
  width = "100%",
  inputPanel(
    prettyRadioButtons(
      inputId = "timescale",
      label = "Time scale", 
      thick = TRUE,
      shape = "curve",
      choices = c("sec", "min"),
      animation = "pulse", 
      status = "primary",
      selected = "sec",
      inline = TRUE,
      fill = TRUE
      #icon = icon("clock-o")
    )
  ),
  inputPanel(uiOutput("solver_parms"))
)
