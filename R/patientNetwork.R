#' Patient Network UI Function (Human Overview)
#'
#' @param id, character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements
#' 
#' @export

patientNetworkUI <- function(id) {
  # Create a namespace function using the provided id
  ns <- shiny::NS(id)
  shiny::column(
    width = 12,
    style = 'padding:0px',
    # the guidlinequestions in the top left box
    bs4Dash::bs4Card(
      title = HTML("<b>Guideline Questions</b>"),
      elevation = 3,
      width = 12,
      height = "100%",
      closable = FALSE,
      collapsible = FALSE,
      headerBorder = FALSE,
      fluidRow(
        tagList(
          HTML(
            "<ol>
          <li>Try to replicate the pressure-volume loop of the left ventricle with the help
          of the 'Model Playground'. 
          Compare the model output to textbook illustrations or search it on the internet.
          Is the human simulator behaving as it should? 
          </li>
          <li>With this human simulator, you can replicate states of blood hypocapnia by increasing the 
          breathing frequency, but it cannot be called hyperventilation (although this also results in hypocapnia). 
          Do you know what functional unit is missing? <i>Hint: Look up what pulmonary stretch
          receptors are for.</i></li>
          <li>How do the resistances in the systemic and cerebral vascular bed change 
          with altering CO2 -levels? </li>
          </ol>"
          )
        )
      ),
      fluidRow(
        column(
          offset = 2,
          width = 4,
          actionBttn(
            inputId = "solution1",
            label = "1",
            style = "jelly",
            size = "sm",
            color = "default",
            icon = icon("lightbulb")
          )
        ),
        column(
          width = 4,
          actionBttn(
            inputId = "solution2",
            label = "2",
            style = "jelly",
            size = "sm",
            color = "default",
            icon = icon("lightbulb")
          )
        ),
        column(
          width = 4,
          actionBttn(
            inputId = "solution3",
            label = "3",
            style = "jelly",
            size = "sm",
            color = "default",
            icon = icon("lightbulb")
          )
        )
      )
    ),
    
    # patient simulator on the right
    bs4Dash::bs4Card(
      title = fluidRow(
        width = 10,
        tagList(
          HTML("<b>Human Overview</b>"), br(),
          h6("The images represent the scope of the human simulator.
                   Doubleclick on the organs and see what pops up."))
      ),
      elevation = 3,
      width = 12,
      closable = FALSE,
      collapsible = FALSE,
      headerBorder = FALSE,
      shiny::div(
        id = "Patient_Overview",
        style =
          "background: url('human_network/human_wholebody.svg');
        background-repeat: no-repeat;background-size: cover;
        background-position: center; width: 100%; height: 70%;",
        shinycssloaders::withSpinner(
          visNetwork::visNetworkOutput(
            ns("patient_overview"), 
            height = "800px"
          ),
          size = 2,
          type = 8,
          color = "darkred",
          color.background = "white"
        )
      )
    )
  )
}


#' Patient Network Server Function
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

patientNetwork <- function(input, output, session) {
  
  ns <- session$ns
  
  #------------------------------------------------------------------------- 
  #  Generate the patient overview network
  #-------------------------------------------------------------------------
  
  nodes_patient <- shiny::reactive({
    generate_nodes_pat()
  })
  
  output$patient_overview <- visNetwork::renderVisNetwork({
    
    nodes_patient <- nodes_patient()
    generate_network(
      nodes = nodes_patient,
      edges = NULL,
      usephysics = TRUE
    )%>%
      # simple click event to select a node
      visNetwork::visEvents(selectNode = "function(nodes) {
      Shiny.onInputChange('current_node_id', nodes.nodes);}"
      ) %>% 
      # very important: change the whole graph position after drawing
      visNetwork::visEvents(
        type = "on",
        stabilized = "function() { 
      this.moveTo({position: {x:-12.96973, y:10.16007}, offset: {x: 0, y:0} })
      }"
      ) %>%
      # add doubleclcick for nodes
      visNetwork::visEvents(
        doubleClick = "function(nodes) {
      Shiny.onInputChange('current_node_bis_id', nodes.nodes);
      }")
  })
  
  #------------------------------------------------------------------------- 
  # Get node position
  #-------------------------------------------------------------------------
  
  vals <- shiny::reactiveValues(
    coords = NULL, 
    viewposition = NULL,
    connected_edges = NULL
  )
  
  # Node position
  # useful to set a proper layout
  shiny::observe({
    shiny::invalidateLater(1000)
    visNetwork::visNetworkProxy(ns("patient_overview")) %>% visNetwork::visGetPositions()
    vals$coords <- if (!is.null(input$patient_overview_positions)) 
      do.call(rbind, input$patient_overview_positions)
  })
  
  # view position (of the camera)
  # useful to set a proper view
  shiny::observe({
    shiny::invalidateLater(1000)
    visNetwork::visNetworkProxy(ns("patient_overview")) %>% visNetwork::visGetViewPosition()
    vals$viewposition <- if (!is.null(input$patient_overview_viewPosition))
      do.call(rbind, input$patient_overview_viewPosition)
  })
  
  positions_patient <- shiny::reactive({
    list(
      position = vals$coords,
      view = vals$viewposition
    )
  })
  
  return(positions_patient)
  
}