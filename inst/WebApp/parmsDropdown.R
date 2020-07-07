# Code for the 'Parametres' dropdown in the 'Model Playground' Panel --> deleted
# parmsDropdown <- dropdown(
#   inputId = "parms_dropdown",
#   label = "Parameters",
#   icon = icon("sliders"),
#   style = "default",
#   status = "default btn-controls",
#   size = "lg",
#   up = TRUE,
#   #animate = TRUE,
#   width = "100%",
#   shiny::tags$div(
#     style = "height: 200px; overflow-y: scroll; overflow-x: hidden;",
#     lapply(1:length(pars_gadget), FUN = function(i) {
#       sliderInput(
#         inputId = paste0(names(pars_gadget)[[i]], "-value"),
#         label = names(pars_gadget)[[i]],
#         value = pars_gadget[[i]],
#         min = 0,
#         max = pars_gadget[[i]] * 10,
#         width = "80%"
#         
#         
#           # sliderInput(
#           #   inputId = "vtid",
#           #   label = "Vtid",
#           #   value = pars_gadget[["Vtid"]],
#           #   min = 0,
#           #   max = pars_gadget[["Vtid"]] * 10,
#           #   width = "80%"
# 
#     # lapply(1:length(pars_dropdown_names), FUN = function(i) {
#     #   sliderInput(
#     #     inputId = paste0(names(pars_dropdown_names)[[i]], "-value"),
#     #     label = names(pars_dropdown_names)[[i]],
#     #     value = pars_dropdown_names[[i]],
#     #     min = 0,
#     #     max = pars_dropdown_names[[i]] * 10,
#     #     width = "80%"
# 
#       ) 
#       # %>%
#         # shinyInput_label_embed(
#         #   icon("undo") %>%
#         #     actionBttn(
#         #       inputId = paste("reset", names(pars_gadget)[[i]], sep = "_"),
#         #       #inputId = paste("reset", names(pars_dropdown_names)[[i]], sep = "_"),
#         #       label = "",
#         #       size = "xs"
#         #     )
#         # )
#     })
#   )
# )
