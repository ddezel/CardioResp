# #' Debug UI Function
# #'
# #' @param id, character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
# #'
# #' @return a \code{shiny::\link[shiny]{tag}} containing UI elements
# #' 
# #' @export
# deBugUI <- function(id) {
#   ns <- shiny::NS(id)
#   bs4Dash::bs4Card(
#     title = "Network Debug",
#     elevation = 3,
#     width = 12,
#     closable = FALSE,
#     collapsible = FALSE,
#     shiny::verbatimTextOutput(ns("debug"))
#   )
# }


# #' Debug Server Function
# #'
# #' @param input Shiny inputs.
# #' @param output Shiny Outputs.
# #' @param session Session object.
# #' @param positions Node coordinates and camera position returned by \link{patientNetwork}
# #'
# #' @return list with following components
# #' \describe{
# #'   \item{xvar}{reactive character string indicating x variable selection}
# #'   \item{yvar}{reactive character string indicating y variable selection}
# #' }
# #' 
# #' @export
# deBug <- function(input, output, session, positions) {
#   #ns <- session$ns
#   output$debug <- shiny::renderPrint(positions())
#   #output$debugKidney <- renderPrint(positionsKidney())
   
# }
