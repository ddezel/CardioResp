#' Model Playground (Gadget) UI Function
#'
#' @param id, character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements
#' 
#' @export

patientGraphUI <- function(id) {
  ns <- shiny::NS(id)
  bs4Dash::bs4Card(
    title = "Gadget Playground",
    elevation = 3,
    width = 12,
    closable = FALSE,
    collapsible = FALSE,
    headerBorder = FALSE,
    style = 'padding_0px'
  )
}

#' Graph Output Server Function
#'
#' @param input Shiny inputs.
#' @param output Shiny Outputs.
#' @param session Session object.
#'
#' @return list with following components
#' \describe{
#'   \item{xvar}{reactive character string indicating x variable selection}
#'   \item{yvar}{reactive character string indicating y variable selection}
#' }
#' 
#' @export

patientGraph <- function(input, output, session) {
  ns <- session$ns
  output$patient_graph <- shiny::renderUI({
     # fluidRow(
     #   column(
     #     width = 6,
     #     style = 'padding:0px;'
     #     #uiOutput("graph_box")
     #   )
     # )
  
    
    ## Add Davids gadget
    #tags$iframe(src="https://cardiomodel.shinyapps.io/gadget/", height=600, width=535)

    #tagList(
     # bs4Card(
        
        # withSpinner(
        #   plotlyOutput(
        #     "plot_node", 
        #     height = "500px", 
        #     width = "100%"
        #   ),
        #   size = 2,
        #   type = 8,
        #   color = "#000000"
        # )
    #  )
    #)
  })
}
