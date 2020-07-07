#' Team Card UI module 
#'
#' @param id character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
#' @param width Card width. 4 by default.
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements
#' 
#' @export
teamCardUI <- function(id, width = 4) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("staffCard"), container = shiny::column, width = width)
}

#' Team Card Server module
#'
#' @param input Shiny inputs.
#' @param output Shiny Outputs.
#' @param session Session object.
#' @param src Path to the profile image.
#' @param name Staff member name.
#' @param position Staff member current position
#'
#' @return a card containing people's information
#' 
#' @export
teamCard <- function(input, output, session, src, name, position) {
  output$staffCard <- shiny::renderUI({
    bs4Dash::bs4Card(
      width = 12,
      solidHeader = TRUE,
      status = "danger",
      title = NULL,
      gradientColor = NULL,
      elevation = 4,
      collapsible = FALSE,
      closable = FALSE,
      bs4Dash::cardProfile(src = src, title = name, subtitle = position)
    )
  })
}