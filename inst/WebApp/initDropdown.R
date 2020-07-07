# Code for the 'Variables' dropdown in the 'Model Playground' Panel
initDropdown <- dropdown(
  inputId = "init_dropdown",
  label = "Variables",
  icon = icon("arrows-alt-v"),
  style = "default",
  status = "default btn-controls",
  size = "lg",
  up = TRUE,
  #animate = TRUE,
  width = "100%",
  shiny::tags$div(
    style = "height: 200px; overflow-y: scroll;",
    lapply(1:length(state_gadget), FUN = function(i) {
      numericInput(
        inputId = paste0(names(state_gadget)[[i]], "-init"),
        label = names(state_gadget)[[i]],
        value = state_gadget[[i]],
        min = 0
      )
    })
  )
)
