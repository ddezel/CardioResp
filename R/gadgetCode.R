#' Patient Network UI Function
#'
#' @param id, character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
#'
#' @return a \code{shiny::\link[shiny]{tag}} containing UI elements
#' 
#' @export

gadgetCodeUI <- function(id) {
  ns <- shiny::NS(id)
  bs4Dash::bs4Card(
    
    # Allow js interactions
    useShinyjs(),
    title = fluidRow(
      column( 
        width = 10,
        tagList(
        HTML("<b>Model Playground</b>"), br(), h6("Choose your timescale, change parameters 
        and plot it all. Have fun!"), h6(tags$i("FYI: If you play with the numbers of the variables, 
        increase the timescale to approx. 5 or 10 minutes to discover the effects."))
        )
      ),
      column( 
        width = 2,
        actionBttn(
          inputId = "resetAll",
          label = "Reset",
          style = "unite",
          size = "sm",
          color = "default",
          icon = icon("undo")
        )
      )
    ),
    width = 12,
    elevation = 4,
    solidHeader = TRUE,
    collapsible = FALSE,
    closable = FALSE,
    headerBorder = FALSE,
    style = 'padding_0px',
    height = "600px",
        fluidRow(
          column(
            width = 12,
            withSpinner(
              plotOutput("plot", height = "380px"), # plot of the gadget
              size = 2,
              type = 8,
              color = "darkred",
              color.background = "white"
            )
          ) 
          ), 
          # organizing all the dropdowns for the gadget
          fluidRow(
            column( 
              width = 4, 
              timeScaleDropdown
              ),
            # Model parameters --> too complicated for user to handle
            # column(
            #   width = 3,
            #   parmsDropdown
            #   ),
            # Initial conditions
            column(
              width = 4, 
              initDropdown
              ),
            # plot options
            column(
              width = 4,
              plotOptsDropdown 
              
              )
          ),
    
        
        tags$div(
          class = "btn-group-charter btn-group-justified-charter",
          # Timescale setup
          
          tags$script("$('.sw-dropdown').addClass('btn-group-charter');"),
          tags$script(HTML("$('.sw-dropdown > .btn').addClass('btn-charter');")),
          tags$script("$('#sw-content-timescale_dropdown').addClass('sw-show');"),
          tags$script("$('#sw-content-plot_dropdown').addClass('sw-show');"),
          #tags$script("$('#sw-content-filterdrop').click(function (e) {e.stopPropagation();});"),
          singleton(tags$head(includeScript(system.file("www", "shiny-utils.js", package = "esquisse"))))
        )
        
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

gadgetCode <- function(input, output, session) {
  ns <- session$ns
  output$gadget_code <- shiny::renderUI({
    # 
    # fluidRow(
    #   column(
    #     
    #   )
    # )

  })
}