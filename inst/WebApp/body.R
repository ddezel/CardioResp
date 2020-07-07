# staff infos
staff_data <- list(
  src = c(
    "http://interfacegroup.ch/content/uploads/2013/02/Diane_larger_200_200.jpg",
    "http://interfacegroup.ch/content/uploads/2017/10/David-Grandjon2_200_200.jpg",
    "http://interfacegroup.ch/content/uploads/2017/09/Isabelle-Rudolf_200_200.jpg"
  ),
  names = c("Diane", "David", "Isabelle"),
  positions = c("Project manager", "Lead Developer", "Developer")
)

# Body functions
body <- bs4DashBody(
  # include css file 
  includeCSS(path = "www/css/cardio_renal_app.css"),
  # javascript file
  includeScript(path = "www/js/screensize.js"),
  
  bs4TabItems(
    # MAIN SECTION
    bs4TabItem(
      tabName = "main",
      fluidRow(
        column(
          width = 5,
          # link to patientNetwork.R module
          patientNetworkUI("patient_network")
        ),
        column(
          width = 7,
          fluidRow(
            gadgetCodeUI("gadget_code")
          ),
          # link to graphs_simulations.R module
          fluidRow(
            graphSimulationsUI("graph_simulations")
            #uiOutput("graph_simulations")
          )
        )
      )
    ),
    # ITEM SECTION
    bs4TabItem(
      tabName = "lexicon",
      
        lexiconUI("lexicon")
    ),
    # ABOUT SECTION
    bs4TabItem(
      tabName = "about",
      fluidRow(column(width = 4, align = "center", br(), h1("Our Team"), br())),
      fluidRow(lapply(seq_along(staff_data), FUN = function(i) teamCardUI(staff_data$names[i]))),
      fluidRow(
        column(
          width = 12,
          br(),
          align = "center",
          div(
            h5( 
              a(
                href = "http://interfacegroup.ch/people/", 
                target = "_blank",
                img(src = "logos/interface_logo.svg", height = "40px")
              ), 
              "The Interface Group, University of Zurich"
            ) 
          )
        )
      )
    )
  )
)
