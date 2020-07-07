#' Lexicon Page
#'
#' @param id, character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
#'
#' @return a \code{shiny::\link[shiny]{tag}} containing UI elements
#' 
#' @export
#' 
lexiconUI <- function(id) {
  ns <- shiny::NS(id)
  fluidPage(
    fluidRow(
      column(
        width = 8,
        # top left card displaying the model graphically
        bs4Dash::bs4Card(
          title = fluidRow(
            column(
              width = 9,
              h5(tags$b("The Underlying Model")) ),
            column(
              width = 3,
              actionButton(
                inputId = "Ellweinpaper",
                label = "Ellwein Paper",
                onclick = "window.open('https://pdfs.semanticscholar.org/e427/f0107ac74652cd6e4dbec4d89d94ee009841.pdf', '_blank')",
                icon = icon("sticky-note")
              )
            )
            ,
            fluidRow(
              style = 'padding:10px',
              h6("The plots and simulations in the main page are all based on the output of this mathematical model, 
               which calculates all of the parameters in the background of this application and then reveals the 
               not-so-easy-to-interpret output in the form of understandable plots. 
               The raw structure is taken from Ellwein and colleagues.")
            )
          ),
          elevation = 3,
          width = 12,
          closable = FALSE,
          collapsible = FALSE,
          headerBorder = FALSE,
          fluidPage(
            div(
              img(src = "lexicon/model_structure.pdf", height = "100%", width = "100%")
            ),
            HTML("<i>The red arrows represent oxygenated arteries, blue arrows the 
          oxygen-poor veins and in green by the SNS affected variables.</i>")
          )
        ),
        # bottom left card displaying the model abbreviations 
        bs4Card(
          title = HTML("<b>Model Abbreviations</b>"),
          width = 12,
          elevation = 3,
          closable = FALSE,
          collapsible = FALSE,
          headerBorder = FALSE,
          div(
            img(src = "lexicon/abbreviations.png", height = "100%", width = "100%")
          )
        ),
        bs4Card(
          title = "If you want to read more on physiological modeling..",
          width = 12,
          elevation = 3,
          closable = FALSE,
          collapsible = FALSE,
          headerBorder = FALSE,
          div(
            h6( 
              a(href = "https://physoc.onlinelibrary.wiley.com/doi/epdf/10.1113/jphysiol.2005.087320", 
                target = "_blank", "Ainslie et al. 2005, Differential responses to CO2 and sympathetic stimulation [..]"
              )
            ),
            h6(
              a(href = "https://repository.lib.ncsu.edu/bitstream/handle/1840.16/3576/etd.pdf?sequence=1&isAllowed=y", 
                target = "_blank", "Ellwein. 2008, Cardiovascular and Respiratory Regulation, Modeling and Parameter Estimation"
              )
            ),
            h6(
              a(href = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4896934/pdf/fphys-07-00189.pdf", 
                target = "_blank", "Fresiello et al. 2016, A Model of the Cardiorespiratory Response [..]"
              )
            )
          )
        )
      ),
      column(
        width = 4,
        # code for the sortable info cards on the right panel
        bs4Sortable(
          bs4Card(
            width = 12,
            title = HTML("<b>Hypercapnia</b>"),
            closable = FALSE,
            collapsible = FALSE,
            headerBorder = FALSE,
            overflow = TRUE,
            elevation = 3,
            "Hypercapnia is an abnormal high concentration of carbon dioxide in the 
          blood. CO2 is a gaseous product of the metabolism and normally quits the 
          body through the lungs. Hypercapnia is usually caused by acute respiratory
          failures (extremely slow breathing)."
          ),
          bs4Card(
            width = 12,
            title = HTML("<b>Hypocapnia</b>"),
            closable = FALSE,
            collapsible = FALSE,
            headerBorder = FALSE,
            overflow = TRUE,
            elevation = 3,
            "Hypocapnia is a state of reduced carbon dioxide in the blood. The decline
          of CO2 concentration usually results from deep or rapid breathing, 
          also known as hyperventilation."
          ),
          bs4Card(
            width = 12,
            title = HTML("<b>Additional information about the model</b>"),
            closable = FALSE,
            collapsible = FALSE,
            headerBorder = FALSE,
            overflow = TRUE,
            elevation = 3,
            ## DDZ height = "100%",
            tagList(
              div(
                class = "card-body",
                style="overflow-y: auto; max-height: 600px;")
              ,
              HTML(
                "Briefly, the cardiovascular system was represented using a limited number of compartments, 
          to minimize the number of parameters and unknowns in the system but still 
          allow for the adequate hemodynamics to be captured. For the <b>heart</b>, we model 
          the right and left ventricles and the four valves but neglect the two atria. 
          On the vascular side, we distinguish between cerebral vasculature and a compartment 
          accounting for the rest of the systemic vasculature. Each <b>vascular bed</b> is modeled 
          using a two element Windkessel model capturing the resistance and compliance 
          provided by the depicted vessels. The <b>lungs</b> are illustrated as one dead space 
          and a gas exchange compartment in the alveoli. <b>Chemosensing</b> is only considered 
          in the cerebral vasculature, chemoreceptors in the aorta are disregarded. 
          Signals concerning the carbon dioxide concentration are transmitted to the 
          medulla oblongata, from where the <b>sympathetic nervous system</b> acts on various 
          players of the model: we included regulations triggered by sympathetic nervous 
          activity in the ventricular pressures, the cerebral and systemic resistances 
          and the vascular compliances. "
              )
            )
          ),
          bs4Card(
            width = 12,
            title = HTML("<b>What is mathematical modeling?</b>"),
            closable = FALSE,
            collapsible = FALSE,
            headerBorder = FALSE,
            overflow = TRUE,
            elevation = 3,
            "A numerical model helps explaining a system, enables to study the effects of the 
          various components and make predictions about the behaviour. A physiological model
          usually contains fixed parameters (cannot be changed, e.g. the weight), initial 
          conditions (starting values for the system; they fluctuate over time, e.g. volumes).
          Both of them are then used to calculate the remaining quantities (e.g. pressures), 
          which are permanently updated."
          )
        )
      )
    )
  )
}


#' Lexicon Server Function
#'
#' @param input Shiny inputs.
#' @param output Shiny Outputs.
#' @param session Session object.
#'
#' @return list with following components
#' @export

lexicon <- function(input, output, session) {
  ns <- session$ns
  output$lexicon 
  
}