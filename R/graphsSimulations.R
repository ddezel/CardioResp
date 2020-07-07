#' Graph Simulations (Hyper/Hypocapnia) Output UI Function
#'
#' @param id, character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
#'
#' @return a \code{shiny::\link[shiny]{tag}} containing UI elements
#' 
#' @export

graphSimulationsUI <- function(id) {
  ns <- shiny::NS(id)
  bs4Dash::bs4Card(
    id = "simulation_tabs",
    title = 
      fluidRow(
        tagList(
        HTML("<b>Simulation-Box</b>"), br(), 
        h6("Click on the buttons and observe the influence of different CO2 levels 
           on the body. For some background information, check out the 'Lexicon' tab 
           on the left sidebar."), 
        h6(tags$i("Hypercapnia was achieved by decreasing the breathing rate from 16 breaths/min to
           ~8 breaths/min, to model hypocapnia the breathing frequency was increased to ~43 breaths/min."))
        )
      ),
      fluidRow(
        column(
        width = 4,
        actionBttn(
          inputId = "iso",
          size = "m",
          label = "Isocapnia",
          style = "jelly",
          color = "danger",
          icon = icon("arrows-alt-h")
          )
        ),
        column(
        width = 4,
        actionBttn(
          inputId = "hypercap",
          size = "m",
          label = "Hypercapnia",
          style = "jelly",
          color = "danger",
          icon = icon("long-arrow-alt-up")
        )
      ),
      column(
        width = 4,
        actionBttn(
          inputId = "hypocap",
          size = "m",
          label = "Hypocapnia",
          style = "jelly",
          color = "danger",
          icon = icon("long-arrow-alt-down")
        )
      )
    ),
    elevation = 3,
    width = 12,
    closable = FALSE,
    collapsible = FALSE,
    headerBorder = FALSE,
    status = "transparent",
    tabStatus = "white",
       bs4TabPanel(
        tabName = " ",
        # Output of the Simulation Box bottom right
        fluidRow(
          column(
            width = 12,
            align = "center",
            uiOutput("graph_simulations")
            # bs4TabSetPanel(
            #   status = "transparent",
            #   tabStatus = "light",
            #    id = "tabset1",
            #    side = "left", 
            #     bs4TabPanel(
            #       tabName = "Resistances and SN-Activity",
            #       active = TRUE,
            #       withSpinner(
            #           plotlyOutput("simulations_plot_1"),
            #         size = 2,
            #         type = 8,
            #         color = "darkred",
            #         color.background = "white"
            #       )
            #      ),
            #      bs4TabPanel(
            #        tabName = "Flow and CO",
            #        active = FALSE,
            #        withSpinner(
            #          plotlyOutput("simulations_plot_2"),
            #          size = 2,
            #          type = 8,
            #          color = "darkred",
            #          color.background = "white"
            #        )
            #      ),
            #      bs4TabPanel(
            #        tabName = "MAP and HR",
            #        active = FALSE,
            #        withSpinner(
            #          plotlyOutput("simulations_plot_3"),
            #          size = 2,
            #          type = 8,
            #          color = "darkred",
            #          color.background = "white"
            #        )
            #      )
            #  )
          )
        )
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

graphSimulations <- function(input, output, session) {
  
  ns <- session$ns
  
 # output$graph_simulations <- shiny::renderUI({
    # 
    # bs4TabSetPanel(
    #   status = "transparent",
    #   tabStatus = "light",
    #   id = "tabset1",
    #   side = "left", 
    #   bs4TabPanel(
    #     tabName = "Resistances and SN-Activity",
    #     active = TRUE,
    #     withSpinner(
    #       # plotlyOutput(
    #       #   if (is.null(input)) "simulations_plot_1" 
    #       #   else if (input$hypercapnia) "plot_hypercapnia_1"
    #       #   else "plot_hypocapnia_1"),
    #       size = 2,
    #       type = 8,
    #       color = "darkred",
    #       color.background = "white"
    #     )
    #   ),
    #   bs4TabPanel(
    #     tabName = "Flow and CO",
    #     active = FALSE,
    #     withSpinner(
    #       plotlyOutput("simulations_plot_2"),
    #       size = 2,
    #       type = 8,
    #       color = "darkred",
    #       color.background = "white"
    #     )
    #   ),
    #   bs4TabPanel(
    #     tabName = "MAP and HR",
    #     active = FALSE,
    #     withSpinner(
    #       plotlyOutput("simulations_plot_3"),
    #       size = 2,
    #       type = 8,
    #       color = "darkred",
    #       color.background = "white"
    #     )
    #   )
    # )
    
    # fluidRow(
    #   column(
    #     width = 6,
    #     style = 'padding:0px;'
    #     #uiOutput("graph_box")
    #   )
    # )
    
  
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
  #})
}
