#' Function to create a visNetwork object in the 'Human Overview' Panel
#'
#' @param nodes Network nodes. A dataframe.
#' @param edges Network edges. A dataframe.
#' @param usephysics Whether to use visNetwork physics. FALSE by default.
#'
#' @return a visNetwork object.
#' 
#' @export
generate_network <- function(nodes, edges, usephysics = FALSE) {
  
  visNetwork::visNetwork(
    nodes, 
    edges, 
    width = "100%", 
    height = "100%"
  ) %>%
    visNetwork::visNodes(
      shapeProperties = list(
        interpolation = TRUE
      )
    ) %>%
    # put shadow on false
    visNetwork::visEdges(
      shadow = FALSE, 
      smooth = TRUE,
      font = list(align = "horizontal")
    ) %>%
    # add group selection option
    visNetwork::visOptions(
      highlightNearest = FALSE, 
      clickToUse = FALSE, 
      manipulation = FALSE, 
      collapse = FALSE,
      autoResize = TRUE
    ) %>% 
    # prevent edge from being selected when a node is selected
    visNetwork::visInteraction(
      hover = TRUE, 
      hoverConnectedEdges = FALSE, 
      selectConnectedEdges = FALSE, 
      multiselect = FALSE, 
      dragNodes = TRUE,
      dragView = FALSE, 
      zoomView = FALSE,
      navigationButtons = FALSE,
      selectable = TRUE,
      tooltipStyle = 
        'position: fixed;
         visibility:hidden;
         padding: 5px;
         padding-right: 5px;
         padding-bottom: 5px;
         white-space: nowrap;
         font-family: Calibri;
         font-size:12px;
         font-color:#9fb2c6;
         background-color: #FFFFFF;
         -moz-border-radius: 3px;
         -webkit-border-radius: 3px;
         border-radius: 3px;
         border: 2px solid #e6f2ff;
         box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);
         z-index: 100;',
      tooltipStay = 300
    ) %>% 
    # stabilization prevents arrows from bouncing
    visNetwork::visPhysics(
      stabilization = TRUE, 
      enabled = usephysics
    )
}

#' Function to create nodes for a patient
#' 
#' @return a dataframe containing all nodes
#' @export
generate_nodes_pat <- function() {
  data.frame(
    id = c("ki", "he", "lu", "ur", "br"),
    title = c(rep("doubleclick!",3), NA, "doubleclick!" ),
    shape = rep("image", 5),
    image = c("human_network/kidney.svg",
              "human_network/heart.svg",
              "human_network/lungs.svg",
              "human_network/urine.svg",
              "human_network/brain.svg"
    ),
    x = c(-40, 30, -5, -12, -10),
    y = c(132, -40, -142, 250, -285),
    size = 43,
    physics = rep(FALSE, 5),
    stringsAsFactors = FALSE,
    fixed = list("x" = TRUE, "y" = TRUE)
  )
}


#' Function to create nodes for the kidney network
#' 
#' @return a dataframe containing all nodes
#'
# nodes_kidney_diagram <- function() {
#   data.frame(
#     id = c("ADH1", "ADH2", "Ald1", "Ald2", "ANP", "AT2_1", "AT2_2", "Na", "brain"),
#     shape = rep("image", 9),
#     image = c(
#       rep("human_network/Kidney_zoom/ADH.svg", 2), 
#       rep("human_network/Kidney_zoom/Ald.svg", 2),
#       "human_network/Kidney_zoom/ANP.svg",
#       rep("human_network/Kidney_zoom/AT2.svg", 2),
#       "human_network/Kidney_zoom/Naplus.svg",
#       "human_network/Kidney_zoom/brain.svg"
#     ),
#     size = c(rep(200, 8), 40),
#     x = c(135, 124, 130, 135, 202, 58, -57, 124, -186),
#     y = c(266, -13, 74, 359, 236, 7, 178, -106, -245),
#     physics = rep(FALSE, 9),
#     stringsAsFactors = FALSE,
#     fixed = list("x" = TRUE, "y" = TRUE)
#   )
# }


#' Function to create nodes for the heart network
#' 
#' @return a dataframe containing all nodes
#'
# nodes_heart_diagram <- function() {
#   data.frame(
#     id = c("ADH", "art", "AT2", "brain", "heart", "Na", "veins"),
#     shape = rep("image", 7),
#     image = c(
#       "human_network/Heart_zoom/ADH.svg", 
#       "human_network/Heart_zoom/arteries.svg",
#       "human_network/Heart_zoom/AT2.svg",
#       "human_network/Heart_zoom/brain.svg",
#       "human_network/Heart_zoom/heart_big.svg",
#       "human_network/Heart_zoom/NAplus.svg",
#       "human_network/Heart_zoom/veins.svg"
#     ),
#     size = c(rep(200, 3), 50, 150, rep(200, 2)),
#     physics = rep(FALSE, 7),
#     x = c(60, 109, 260, -100, -149, 79, 140),
#     y = c(143, -227, -6, 120, -212, 208, -77),
#     stringsAsFactors = FALSE,
#     fixed = list("x" = TRUE, "y" = TRUE)
#   )
# }


#' Function to create nodes for the brain network
#' 
#' @return a dataframe containing all nodes
#'
# nodes_brain_diagram <- function() {
#   data.frame(
#     id = c("Ald", "AT2", "brain", "Kplus", "NAplus", "renin"),
#     shape = rep("image", 6),
#     image = c(
#       "human_network/Brain_zoom/Ald.svg",
#       "human_network/Brain_zoom/AT2.svg",
#       "human_network/Brain_zoom/brain.svg",
#       "human_network/Brain_zoom/Kplus.svg",
#       "human_network/Brain_zoom/NAplus.svg",
#       "human_network/Brain_zoom/renin.svg"
#     ),
#     size = c(rep(250, 2), 60, rep(250, 3)),
#     physics = rep(FALSE, 6),
#     x = c(173, -0, -115, 177, 64, -24),
#     y = c(169, 169, -58, -63, -63, -63),
#     fixed = list("x" = TRUE, "y" = TRUE),
#     stringsAsFactors = FALSE
#   )
# }

#------------------------------------------------------------------------- 
#  Generate the zoom information if "+" button is doubleclicked
#-------------------------------------------------------------------------

# load svg to zoom in the kidney
# kidney_zoom <- shiny::HTML("<div id=\"kidneyzoom\">",
#                     "<img height=\"520px\" src=\"human_network/kidney_diagram.svg\" 
#                     width=\"450px\" align=\"center\">", "</div>")


# load svg to zoom in the brain
# brain_zoom <- shiny::HTML("<div id=\"brainzoom\">",
#                     "<img height=\"250px\" src=\"human_network/brain_zoom.svg\" 
#                     width=\"450px\" align=\"center\">", "</div>")
#   
#   shiny::HTML(
#   paste("<div class=\"row\">", "<div class=\"col-sm-6\">",
#         "<a href=\"human_network/brain_diagram.svg\" target=\"_blank\">
#         <img id = \"zoom_image\" width=\"440\" height=\"300\" border=\"0\"
#         align=\"center\"  src=\"human_network/brain_diagram.svg\"/></a>",
#         "</div>", "<div class=\"col-sm-6\">"
#   )
# )

# load svg to zoom in the heart
# heart_zoom <- shiny::HTML("<div id=\"heartzoom\">",
#                    "<img height=\"480px\" src=\"human_network/heart_zoom.svg\" 
#                     width=\"450px\" align=\"center\">", "</div>")
  
  # paste("<div class=\"row\">", "<div class=\"col-sm-6\">",
  #       "<a href=\"human_network/heart_diagram.svg\" target=\"_blank\">
  #       <img id = \"zoom_image\" width=\"440\" height=\"480\" border=\"0\"
  #       align=\"center\"  src=\"human_network/heart_diagram.svg\"/></a>",
  #       "</div>", "<div class=\"col-sm-6\">"
  # )
# )


# lung_zoom <- shiny::HTML(
#   paste("Some information about the lungs."
#   )
# )

