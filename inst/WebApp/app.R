#------------------------------------------------------------------------- 
#  This codes loads all packages needed by the application.
#  Moreover, it contains all mandatory UI and server elements. 
#-------------------------------------------------------------------------

# load packages
library(shiny)
library(plotly)
library(deSolve)
require(visNetwork)
library(shinyjs)
library(shinycssloaders)
library(shinyjqui)
library(bsplus)
library(purrr)
library(stringr)
library(shinyFeedback)
library(shinyWidgets)
library(bs4Dash)
library(dplyr)
library(shinyEffects)
library(latex2exp)
library(ggplot2)
library(readr)
library(magrittr)
library(FME)
library(htmltools)
library(colorspace)

library(CardioResp)

# define color palette for plots
colors <- c("royalblue1", "red4", "brown2", "gray40", "navy", "indianred1", "black")

# load the modules
source("simulationPlotUtils.R")


# UI components for Gadget
source("parmsDropdown.R")
source("timeScaleDropdown.R")
source("initDropdown.R")
source("plotOptsDropdown.R")
#source("parmsDropdown.R")    # removed

# Load the UI components for the Web App
source("header.R")
source("sidebar.R")
source("body.R")


# set the current time zone to Zurich (for shiny server)
Sys.setenv(TZ = "Europe/Zurich")

# compile the C code containing equations
# system("R CMD SHLIB Ellwein_circulation_model.c")
# dyn.load(paste("Ellwein_circulation_model", .Platform$dynlib.ext, sep = ""))

#------------------------------------------------------------------------- 
#  The Gadget Skeleton
#-------------------------------------------------------------------------

# model directory
coreFolder <- system.file("Ellwein_Gadget/model_core", package = "CardioResp")
#con <- file(paste0(coreFolder, "/equations.c"), open = "r")
#con <- file("/Users/isabellerudolf\ 1/Desktop/VirtualPatient/inst/Ellwein_Gadget/model_core/equations.c", open = "r")
#temp_model <- readLines(con)

# A function which extracts all inpuId of all sliders only for parameters 
# that are not events. It is needed in the observeEvent to reset slider.
# I used this function to avoid passing big vector to the event handler.
extract_slider_input <- function(input, parameters) {
  res <- lapply(1:length(parameters), function(i) {
    paste("input$reset", names(parameters)[[i]], sep = "_")
  })
  res_vec <- unlist(res)
  res_eval <- lapply(1:length(res_vec), function(i) {
    eval(parse(text = res_vec[i]))
  })
  unlist(res_eval)
}

# Initialize the parameter change dataframe
#pars_change <- data.frame(pars_name = NULL, percent_change = NULL)


#------------------------------------------------------------------------- 
#  The App Skeleton
#-------------------------------------------------------------------------

# Define UI
ui <- bs4DashPage(
  sidebar_collapsed = TRUE,
  title = "Cardio-Respiratory Human Simulator",
  navbar = navbar,
  sidebar = sidebar,
  body = body,
  footer = NULL
  # controlbar = dashboardControlbar
)

  #####################
  #    Server code    # 
  #####################

server <- function(input, output, session) {
  
  #------------------------------------------------------------------------- 
  # call Modules
  # callModule(module server function, "id" = same as the id for UI)
  #-------------------------------------------------------------------------
  callModule(patientNetwork, "patient_network")
  callModule(patientGraph, "patient_graph")
  # callModule(lexicon, "lexicon")
  # callModule(deBug, "debug", positions)
  # callModule(graphSimulations, "graph_simulations")
   callModule(gadgetCode, "gadget_code")
  
  lapply(seq_along(staff_data), FUN = function(i) {
    callModule(
      module = teamCard, 
      id = staff_data$names[i], 
      src = staff_data$src[i], 
      name = staff_data$names[i], 
      position = staff_data$positions[i]
    )
  })
  
  # display wholebody human as background
  observe({
    addClass(id = "Patient_Overview", class = "patient_overview")
  })
  
  # define reactive values to modify them later
  events <- reactiveValues(
    logged = FALSE
  )
  
  # When WebApp starts, show sweetalert message for guidline questions
  # Delay of 0.2 seconds 
  observe({
    if (!events$logged) {
      shinyjs::delay(
        100,
        confirmSweetAlert(
          session,
          inputId = "start_app",
          title = "Hello you!",
          text = tagList(
            # author credit for intro_heart.svg: Icons smashicons, flaticons 
            img(src = "human_network/intro_lungs.svg", width = "70px", height = "70px"), HTML("&nbsp;"),
            img(src = "human_network/intro_exchange.svg", width = "30px", height = "30px"), HTML("&nbsp;"),
            img(src = "human_network/intro_heart.svg", width = "70px", height = "70px"),
            br(),
            br(),
            HTML(
              "This web application is a virtual human simulator with the focus on the <b>cardiovascular</b>
              and the <b>respiratory system</b>. You will learn how the breathing is connected to other body functions
              and what role the sympathetic nervous system plays in all of this.", "<br>", "<br>",
              "With the help of various tools you'll be able to <b>answer the three questions</b> asked in the box on the top left of the window."
            ),
            hr(),
            column(
              align = "center",
              width = 12
            )
          ),
          btn_labels = c(NULL, "Alright"),
          type = "info",
          html = TRUE
        )
      )
    }
  })
  
   
  #------------------------------------------------------------------------- 
  # generate plots for simulation box (bottom right in webapp) 
  #-------------------------------------------------------------------------
  
  # render isocapnic simulation plots
  output$plot_isocapnia_1 <- renderPlotly({
  plotly::subplot(
    resistances_plot, symp_act_plot,
    titleX = TRUE,
    titleY = TRUE,
    nrows = 2,
    margin = c(0.1, 0.07,0.15, 0.15),
    heights = rep(0.5, 2)
  ) %>% 
    layout(showlegend = TRUE, plot_bgcolor='rgb(245, 245, 245)')
  })
  
  
  output$plot_isocapnia_2 <- renderPlotly({
    plotly::subplot(
      flow_plot, CO_plot,
      titleX = TRUE,
      titleY = TRUE,
      nrows = 2,
      margin = c(0.1, 0.07,0.15, 0.15),
      heights = rep(0.5, 2)
    ) %>% layout(showlegend = TRUE, plot_bgcolor='rgb(245, 245, 245)')
  })
  
  output$plot_isocapnia_3 <- renderPlotly({
    plotly::subplot(
      HR_plot, MAP_plot, 
      titleX = TRUE,
      titleY = TRUE,
      nrows = 2,
      margin = c(0.1, 0.07,0.15, 0.15),
      heights = rep(0.5, 2)
    ) %>% layout(showlegend = FALSE, plot_bgcolor='rgb(245, 245, 245)')
  })
  
  # render hypercapnic simulation plots
  output$plot_hypercapnia_1 <- renderPlotly({
    plotly::subplot(
      resistances_plot_hyper, symp_act_plot_hyper,
      titleX = TRUE,
      titleY = TRUE,
      nrows = 2,
      margin = c(0.1, 0.07,0.15, 0.15),
      heights = rep(0.5, 2)
    ) %>% layout(showlegend = TRUE, plot_bgcolor='rgb(245, 245, 245)')
    })
  
  output$plot_hypercapnia_2 <- renderPlotly({
  plotly::subplot(
    flow_plot_hyper, CO_plot_hyper,
    titleX = TRUE,
    titleY = TRUE,
    nrows = 2,
    margin = c(0.1, 0.07,0.15, 0.15),
    heights = rep(0.5, 2)
  ) %>% layout(showlegend = TRUE, plot_bgcolor='rgb(245, 245, 245)')
  })
  
  output$plot_hypercapnia_3 <- renderPlotly({
  plotly::subplot(
    HR_plot_hyper, MAP_plot_hyper, 
    titleX = TRUE,
    titleY = TRUE,
    nrows = 2,
    margin = c(0.1, 0.07,0.15, 0.15),
    heights = rep(0.5, 2)
  ) %>% layout(showlegend = FALSE, plot_bgcolor='rgb(245, 245, 245)')
  })
    
  # render hypocapnic simulation plots
  output$plot_hypocapnia_1 <- renderPlotly({
   plotly::subplot(
      resistances_plot_hypo, symp_act_plot_hypo,
      titleX = TRUE,
      titleY = TRUE,
      nrows = 2,
      margin = c(0.1, 0.07,0.15, 0.15),
      heights = rep(0.5, 2)
    ) %>% 
      layout(showlegend = TRUE, plot_bgcolor='rgb(245, 245, 245)')
  })
  
  output$plot_hypocapnia_2 <- renderPlotly({
   plotly::subplot(
      flow_plot_hypo, CO_plot_hypo,
      titleX = TRUE,
      titleY = TRUE,
      nrows = 2,
      margin = c(0.1, 0.07,0.15, 0.15),
      heights = rep(0.5, 2)
    ) %>% layout(showlegend = TRUE, plot_bgcolor='rgb(245, 245, 245)')
    
  })
  
  output$plot_hypocapnia_3 <- renderPlotly({
   plotly::subplot(
      HR_plot_hypo, MAP_plot_hypo, 
      titleX = TRUE,
      titleY = TRUE,
      nrows = 2,
      margin = c(0.1, 0.07,0.15, 0.15),
      heights = rep(0.5, 2)
    ) %>% layout(showlegend = FALSE, plot_bgcolor='rgb(245, 245, 245)')
    
  })
  
    # render the plots for the 'Isocapnia' actionbutton
    observeEvent(input$iso, {
      if(input$iso) {
      output$graph_simulations <- shiny::renderUI({
      if (input$iso) {
        bs4TabSetPanel(
        status = "transparent",
        tabStatus = "light",
        id = "tabset1",
        side = "left", 
        bs4TabPanel(
          tabName = "Resistances & SN-Activity",
          active = TRUE,
          withSpinner(
            plotlyOutput("plot_isocapnia_1"),
            size = 2,
            type = 8,
            color = "darkred",
            color.background = "white"
          )
        ),
        bs4TabPanel(
          tabName = "Flow & Cardiac Output",
          active = FALSE,
          withSpinner(
            plotlyOutput("plot_isocapnia_2"),
            size = 2,
            type = 8,
            color = "darkred",
            color.background = "white"
          )
        ),
        bs4TabPanel(
          tabName = "MAP & Heart Rate",
          active = FALSE,
          withSpinner(
            plotlyOutput("plot_isocapnia_3"),
            size = 2,
            type = 8,
            color = "darkred",
            color.background = "white"
          )
        )
       )
      }
    })
   }
  })
    
    # render the plots for the 'Hypercapnia' actionbutton
    observeEvent(input$hypercap, {
      if(input$hypercap) {
     output$graph_simulations <- shiny::renderUI({
      bs4TabSetPanel(
        status = "transparent",
        tabStatus = "light",
        id = "tabset1",
        side = "left",
        bs4TabPanel(
          tabName = "Resistances & SN-Activity",
          active = TRUE,
           withSpinner(
            plotlyOutput("plot_hypercapnia_1",
                         height = "400px",
                         width = "100%"),
            size = 2,
            type = 8,
            color = "blue",
            color.background = "white"
          )
        ),
        bs4TabPanel(
          tabName = "Flow & Cardiac Output",
          active = FALSE,
          withSpinner(
            plotlyOutput("plot_hypercapnia_2",
                         height = "400px",
                         width = "100%"),
            size = 2,
            type = 8,
            color = "blue",
            color.background = "white"
          )
        ),
        bs4TabPanel(
          tabName = "MAP & Heart Rate",
          active = FALSE,
          withSpinner(
            plotlyOutput("plot_hypercapnia_3",
                         height = "400px",
                         width = "100%"),
            size = 2,
            type = 8,
            color = "blue",
            color.background = "white"
          )
        )
      )
    })
    }
  })
    
    # render the plots for the 'Hypocapnia' actionbutton
    observeEvent(input$hypocap, {
      if(input$hypocap) {
      output$graph_simulations <- shiny::renderUI({

        bs4TabSetPanel(
          status = "transparent",
          tabStatus = "light",
          id = "tabset1",
          side = "left",
          bs4TabPanel(
            tabName = "Resistances & SN-Activity",
            active = TRUE,
            withSpinner(
              plotlyOutput("plot_hypocapnia_1",
                           height = "400px",
                           width = "100%"),
              size = 2,
              type = 8,
              color = "dodgerblue",
              color.background = "white"
            )
          ),
          bs4TabPanel(
            tabName = "Flow & Cardiac Output",
            active = FALSE,
            withSpinner(
              plotlyOutput("plot_hypocapnia_2",
                           height = "400px",
                           width = "100%"),
              size = 2,
              type = 8,
              color = "dodgerblue",
              color.background = "white"
            )
          ),
          bs4TabPanel(
            tabName = "MAP & Heart Rate",
            active = FALSE,
            withSpinner(
              plotlyOutput("plot_hypocapnia_3",
                           height = "400px",
                           width = "100%"),
              size = 2,
              type = 8,
              color = "dodgerblue",
              color.background = "white"
            )
          )
        )
      })
      }
    })
  # 
  # # delete compiled files right after session is closed...
  # session$onSessionEnded(function() {
  #   if (.Platform$OS.type == "unix") {
  #     file.remove("Ellwein_circulation_model.o")
  #     file.remove("Ellwein_circulation_model.so")
  #   } else if (.Platform$OS.type == "windows") {
  #     file.remove("Ellwein_circulation_model.dll")
  #   }
  # })
  
  #------------------------------------------------------------------------- 
  # Sweet alert functions are activated if organs are double-clicked
  #-------------------------------------------------------------------------
  
  observeEvent(input$current_node_bis_id, {
     # print(input$current_node_bis_id)

     # show detailed diagram if brain is doubleclicked
     if(input$current_node_bis_id == "br") {
       sendSweetAlert(
         session = session,
         title = "The pathway of CO2-sensing in the brain (medulla oblongata) 
         and the periphery",
         html = TRUE,
         btn_labels = FALSE,
         type = "info",
         closeOnClickOutside = TRUE,
         text = shiny::HTML("In green are the mechanisms that happens if the concentration of CO2 is 
         too high in the blood, in red if the concentration is too low.",
           "<div id=\"brainzoom\">",
                            "<img height=\"100%\" src=\"human_network/brain_zoom.pdf\" 
                    width=\"100%\" align=\"center\">", "</div>")
         #text = brain_zoom    # details in networks.R
       )
     }

     # show detailed diagram if heart is doubleclicked
     if(input$current_node_bis_id == "he") {
       sendSweetAlert(
         session = session,
         title = "The heart",
         html = TRUE,
         btn_labels = FALSE,
         type = "info",
         closeOnClickOutside = TRUE,
         text = shiny::HTML("This image is a  representation of the heart in this simulator.
         Note that it is partly limited: there are no atriums and no baroreceptors in this 
         model (see gray shaded areas). Latter are responsible for the maintenance of a healthy 
         blood pressure via <i>Baroreflex</i>. Their absence indicates that there is
         no negative feedback loop to regulate the blood pressure.",
                            "<div id=\"heartzoom\">",
                            "<img height=\"480px\" src=\"human_network/heart_model.svg\" 
            width=\"450px\" align=\"center\">", "</div>")
       )
     }
     # show kidney infobox if ki node is doubleclicked
    if(input$current_node_bis_id == "ki") {
      sendSweetAlert(
        session = session,
        title = "The kidney is is not part of this model.",
        html = TRUE,
        btn_labels = FALSE,
        type = "info",
        closeOnClickOutside = TRUE,
        text = shiny::HTML("Nevertheless, here is a nice schematic
            of a nephron (the smallest functional unit of the kidney) and its hormonal interactions.",
                           "<div id=\"kidneyzoom\">",
                           "<img height=\"520px\" src=\"human_network/kidney_diagram.svg\" 
           width=\"450px\" align=\"center\">", "</div>")
      )
    }
     # show lung infobox if lu node is doubleclicked
    if(input$current_node_bis_id == "lu") {
      sendSweetAlert(
        session = session,
        title = "How the lungs work",
        html = TRUE,
        btn_labels = FALSE,
        type = "info",
        closeOnClickOutside = TRUE,
        text = shiny::HTML("The respiratory system is responsible for the gas exchange with the blood, 
         where oxygen is transported via atmospheric air into the blood stream while carbon dioxide is removed from
         it and breathed out.",
                           "<div id=\"kidneyzoom\">",
                           "<img height=\"520px\" src=\"human_network/lungs_model.svg\" 
           width=\"450px\" align=\"center\">", "</div>") 
      )
    }
    
    if(input$current_node_bis_id == "ur") {
      sendSweetAlert(
        session = session,
        title = "This model does not account for any form of excretion.",
        html = TRUE,
        btn_labels = FALSE,
        type = "info",
        closeOnClickOutside = TRUE,
        text = shiny::HTML(
          paste("But here is a fun fact about our waste product management in 
          the body:", br(), "The kidney produces <b>daily around 170 litres of primary urine</b> by 
                ultrafiltration of the plasma. But no worries, approximately 99% 
                of it is reabsorbed.")
        ) 
      )
    }
  })
  
  #------------------------------------------------------------------------- 
  # Sweet alert functions are activated if solution buttons are klicked
  #-------------------------------------------------------------------------
  
    observeEvent(input$solution1, {
      sendSweetAlert(
        session = session,
        title = "Solution to Question 1",
        html = TRUE,
        btn_labels = FALSE,
        type = "info",
        closeOnClickOutside = TRUE,
        text = shiny::HTML("<div id=\"solution2\">",
                           "<img height=\"480px\" src=\"solutions/solution2.png\" 
                    width=\"400px\" align=\"center\">", "</div>", "The top 
                    figure is from a textbook, the bottom is the result of the model. 
                    Due to different scaling, they might look quite different, but after
                    a closer look it gets obvious that both loops resemble in terms of quantities
                    and shape. It can therefore be concluded that <b>the model reveals feasible 
                    results and behaves properly.</b>") 
      )
    })
    
    observeEvent(input$solution2, {
      sendSweetAlert(
        session = session,
        title = "Solution to Question 2",
        html = TRUE,
        btn_labels = FALSE,
        type = "info",
        closeOnClickOutside = TRUE,
        text = shiny::HTML("The human simulator does <b>not account for pulmonary stretch 
                         receptors</b> in the bronchia, which belong to the family of mechanoreceptors. 
                         Those would be needed to respond to physical stretching 
                         of the airways (=breathing) by signaling the state of the 
                         lungs to the respiratory control center in the brain. 
                         If the <b>receptors are excessively stretched during large 
                         inspirations</b>, as it happens during hyperventilating, the brain sends <b>signals to inhibit inspiration</b>, 
                         allowing expiration to occur. Thereby the breathing is brought 
                         back to a normal frequency. This reflex is called
                         <i>Hering-Breuer reflex</i>.", "<br><br/>",
                           "<div id=\"solution1\">",
                           "<img height=\"100%\" src=\"solutions/stretchreceptors.pdf\" 
                    width=\"100%\" align=\"center\">", "</div>"
        )
      )
    })
    
    observeEvent(input$solution3, {
      sendSweetAlert(
        session = session,
        title = "Solution to Question 3",
        html = TRUE,
        btn_labels = FALSE,
        type = "info",
        closeOnClickOutside = TRUE,
        text = shiny::HTML("Voilà, the three states in comparison:",
                           "<div id=\"solution3.1\">",
                           "<img height=\"430px\" src=\"solutions/solution3.png\" 
                          width=\"440px\" align=\"center\">", "</div>",
                           "The model reveals that an increase in blood carbon dioxide (<i>hypercapnia</i>) promotes the
                         <b>systemic resistance to increase</b> whereas the <b>cerebral resistance decreases</b>.
                         In the case of low carbon dioxide levels (<i>hypocapnia</i>), the responses are <b>reversed</b>.", "<br>",
                           "Have a look at the 'Brain' node in the human overview box to see the whole pathway."
                           
        )
      )
    })
    
    
  #------------------------------------------------------------------------- 
  # Server code for gadget
  #-------------------------------------------------------------------------
    
  #-------------------------------------------------------------------------
  #
  #  Integrate equations using deSolve package to generate table
  #  out is a reactive intermediate component that is called by
  #  to make plots or other stuffs. We used the compiled version of
  #  the code, to make computations faster
  #
  #-------------------------------------------------------------------------
  
    
  
    pars_gadget <- reactive({
     c(
      context = 0,
      # Parameters rates, tissue volumes, and gas dissociation constants
      M_CO2  = 259.98, #4.333 * 60,  #  [mLSTPD/min] Systemic tissue metabolic rate of CO2 4.333
      M_O2   = 310.02, #5.167 * 60,  #  [mLSTPD/min] Systemic tissue metabolic rate of O2 5.167
      MB_CO2 = 42, #0.7 * 60,        #  [mLSTPD/min] Cerebral tissue metabolic rate of CO2 0.875
      MB_O2  = 52.5, #0.875 * 60,    #  [mLSTPD/min] Cerebral tissue metabolic rate of O2 = MC_CO2
      MS_CO2 = 217.98, #3.633 * 60,  #  [mLSTPD/min] Systemic (body) tissue metabolic rate of CO2 (M_CO2 - MB_CO2)
      MS_O2  = 259.98,       #  [mLSTPD/min] Systemic (body) tissue metabolic rate of O2 (M_O2 - MB_O2)
      VT_CO2 = 15000,        #  [mLSTPD] Systemic tissue volume of CO2 15000
      VT_O2  = 6000,         #  [mLSTPD] Systemic tissue volume of O2 6000
      VB_CO2 = 900,          #  [mLSTPD] Cerebral tissue volume of CO2 900
      VB_O2  = 1000,         #  [mLSTPD] Cerebral tissue volume of O2 1000
      VS_CO2 = 14100,        #  [mLSTPD] Systemic tissue volume of CO2 VT_CO2 - VB_CO2
      VS_O2  = 5000,         #  [mLSTPD] Systemic tissue volume of O2 VT_O2  - VB_O2
      VA_CO2 = 3200,         #  [mLBTPS] Alveolar tissue volume of CO2 3200
      VA_O2  = 2500,         #  [mLBTPS] Alveolar tissue volume of O2 2500
      VD     = 151,          #  [mLBTPS] Total dead space volume 181.1, optimized value
      K1     = 0.2,          #  [mLSTPD mL] Dissociation coefficient for O2 0.2
      K2     = 0.046,        #  [mmHg−1 ] Dissociation coefficient for O2 0.046
      KCO2   = 0.0065,       #  [mLSTPDmmHg mL] Dissociation coefficient for CO2 0.0065
      kCO2   = 0.244,        #  [mLSTPD mL] Dissociation coefficient for CO2 0.244
      qp     = 4920, #82*60, #  [mL/min] Mean pulmonary flow 82*60 ml/min
      pv_CO2 = 45, #46,      #  [mmHg] venous CO2 partial pressure [h]
      pv_O2  = 35, # 40,     #  [mmHg] venous O2 partial pressure [h]
      pi_O2  = 159,          #  [mmHg] air O2 partial pressure at sea level
      pi_CO2 = 0.3,          #  [mmHg] air CO2 partial pressure at sea level
      
      tau_CO = 0.25,          # Cardiac output time constant min
      #Caorta = 1.35,          # [mL/mmHg] Aortic compliance 
      #Raorta = 0.01,          # [mmHg min/mL] Aortic resistance
      d_aorta_0 = 40,         # [mm] Nominal aortic diameter (thorax)
      L_aorta = 47,           # [mm] Nominal aortic lenght (portion of thoracic aorta)
      tau_strain = 4.2,       # [min] Aortic strain time constant
      f_0 = 15,               # [min-1] Baroreceptor gain parameter
      delta_strain_0 = 0.5,   # [without units] Baroreceptor saturation constant
      a = 3.9,                # [min-1] Baroreceptor activation rate
      b = 12.02,              # [min-1] Baroreceptor deactivation rate
      f_SN = 0.525,           # [min-1] Baroreflex arc parameter
      
      tau_autoreg = 6.77,      # [min] Autoregulation time constant
      CO_0 = 2889.769,         # [ml.min-1] Autoregulation param
      CO_1 = 2942.761,         # [ml.min-1] Autoregulation param
      
      P1 = 20.12,            # [mmHg] Steady-state renin-angiotensin system tone
      P2 = 24.98,            # [mmHg] Steady-state renin-angiotensin system tone
      tau_renin = 12.61,     # [min] Time constant for renin production
      tau_at = 1.117,        # [min] Time constant for AT2 production
      tau_al = 30,           # [min] Time constant for aldosterone production
      tau_adh = 6,           # [min] Time constant for ADH production
      tau_map = 0.25,        # [min] Time constant for mean pressure calculation
      
      alpha1 = 0.319,   # Arterial and venous compliance parameter
      alpha2 = 14.18,   # Arterial resistance parameter
      
      Vtid = 600,         # [ml]
      resp2CC_ratio = 0.2 # ratio of the respiratory cycle to cardiac cycle duration: slow breathing: 0.08; //hyperventilation: 0.8; 
    )
  })
    
    
  # Allow js interactions
  useShinyjs()
  
  # model file and first compilation
  # if( dir.exists("temp_core") ) {system("rm -r temp_core")}
  # dir.create("temp_core")
  # file.copy(file.path(coreFolder, list.files(coreFolder)), "temp_core")
  system("R CMD SHLIB temp_core/equations.c")
  dyn.load(paste("temp_core/equations", .Platform$dynlib.ext, sep = ""))
  
  # update the model when the user change equations via the gadget and re-compile
  observeEvent(input$confirm_model, {
    writeLines(input$ellwein_model, "temp_core/equations.c")
    # model compilation
    system("R CMD SHLIB temp_core/equations.c")
    dyn.load(paste("temp_core/equations", .Platform$dynlib.ext, sep = ""))
  })
  
  # model equations
  output$model_eq <- renderUI({
    aceEditor(
      outputId = "ellwein_model", 
      mode = "c_cpp", 
      theme = "tomorrow_night_bright",
      autoComplete = "live",
      height = if (input$screenSize$height <= 700) {
        "400px"
      } else if (input$screenSize$height <= 800) {
        "500px"
      } else {
        "700px"
      }, 
      value = paste(temp_model, sep = "", collapse = "\n")
    )
  })
  
  # define initial time
  times_gadget <- reactive(seq(0, input$tmax, by = 0.0001))  # by = input$dt
  
  # define parameters
  # parameters <- reactive({
  #   temp_parms <- lapply(1:length(pars_gadget), FUN = function(i) {
  #     input[[paste0(names(pars_gadget)[[i]], "-value")]]
  #   })
  #   names(temp_parms) <- names(pars_gadget)
  #   return(temp_parms)
  # })

  
  # define initial conditions
  init <- reactive({
    temp_state <- lapply(1:length(state_gadget), FUN = function(i) {
      input[[paste0(names(state_gadget)[[i]], "-init")]]
    })
    names(temp_state) <- names(state_gadget)
    return(temp_state)
  })
  
  # find solutions
  out_gadget <- reactive({
    pars_gadget <- pars_gadget()
    #state_gadget <- state_gadget()
    input$confirm_model
    req(input$tmax, 0.0005) #0.0001)  # input$dt
    as.data.frame(
      ode(
        y = unlist(init()),
        #y = state_gadget,
        times = times_gadget(),
        func = "derivs",
        #method = "impAdams_d", # impAdams_d is 3 times faster than lsoda here but less precise
        #parms  = unlist(parameters()),
        parms  = pars_gadget,
        dllname = "equations",
        initfunc = "initmod",
        nout = 53,
        outnames = c(
          "Qtv", "Qpv", "Qmv", "Qav", "Qaa", "Qp", "Qs",
          "Qca", "Qc", "Qcv",
          "Ppa", "Ppv", "Paorta", "Psa", "Psv", "Pca", "Pcv",
          "Plv", "Prv",
          "Vtot",
          "A_aorta_0", "A_aorta", "Aortic_strain",
          "delta_strain", "f_Baro", "HR", "P_ra", 
          "Pb_0","Pa_0", "Rc_pB_CO2", "Rs_pa_CO2", "Rs", "Rbrain", "Raorta", "Rc", "Rtot",
          "Csa", "Csv", "ca_CO2", "ca_O2", "cv_CO2", "cv_O2", "cS_CO2", "cS_O2",
          "cB_CO2", "cB_O2", "Va", "d_Va", "Rav", "Rmv", "Rpv", "Rtv", "Raa"
        )
        
        # nout = 40,
        # outnames = c(
        #   "Qtv", "Qpv", "Qmv", "Qav", "Qaa", "Qp", "Qs",
        #   "Qca", "Qc", "Qcv",
        #   "Ppa", "Ppv", "Paorta", "Psa", "Psv", "Pca", "Pcv",
        #   "Plv", "Prv", "Vtot", "HR",
        #   "Rc_pB_CO2", "Rs_pa_CO2", "Rs", "Rbrain", "Raorta", "Rtot",
        #   "Csa", "Csv", "ca_CO2", "ca_O2", "cv_CO2", "cv_O2", "ConcSystemicCO2", "ConcSystemicO2",
        #   "ConcCerebralCO2", "ConcCerebralO2", "VolAlveolar", "d_VolAlveolar", "Raa"
        # )
      )
    )
  })
  
  # Solver UI
  output$solver_parms <- renderUI({
    req(input$timescale)
    # Solver parameters
    tagList(
      numericInput(
        label = "Max. time",
        inputId = "tmax",
        width = "45%",
        value = if (input$timescale == "sec" ) 0.5 else 1,
        min = 0
      )
      # ,
      # numericInput(
      #   label = "Step size",
      #   inputId = "dt",
      #   width = "45%",
      #   value = if (input$timescale == "sec" ) 1e-005 else 0.001,
      #   min = 0
      # )
    )
  })
  
  # returns a table of changed parameters
  # parms_change <- reactive({
  #   for(i in 1:length(pars)) {
  #     if (pars[i] != unlist(parameters()[i])) {
  #       change <- (unlist(parameters()[i]) - pars[i]) / pars[i] * 100
  #       # significant change
  #       if (abs(change) > 1) {
  #         pars_change_temp <- data.frame(
  #           pars_name = names(pars)[i],
  #           percent_change = paste0(change, " %")
  #         )
  #         pars_change <- rbind(pars_change, pars_change_temp)
  #       }
  #     }
  #   }
  #   pars_change
  # })
  # 
  # output$pars_change <- renderPrint({ if (nrow(parms_change()) > 0) parms_change() })
  
  
   output$plot <- renderPlot({
     out_gadget <- out_gadget()[-c(1:2000), ]
     #out_gadget <- out_gadget()
     req(input$xaxis, input$yaxis)
    # check if multiple output variables are selected
     if (length(input$yaxis) == 1) {
       plot(
         x = if (input$timescale == "sec" && input$xaxis == "time" ) {
           out_gadget[[input$xaxis]] * 60
         } else {
           out_gadget[[input$xaxis]]
         },
         y = out_gadget[[input$yaxis]],
         type = "l",
         col = "skyblue1",
         lwd = 4,
         xlab = if (input$timescale == "sec" && input$xaxis == "time" ) {
           paste0(input$xaxis, " [sec]")
         } else if (input$timescale == "min" && input$xaxis == "time" ) {
           paste0(input$xaxis, " [min]" )
         } else {
           input$xaxis
         },
         ylab = input$yaxis,
         ylim = c(min(out_gadget[[input$yaxis]]) * 0.9, max(out_gadget[[input$yaxis]]) * 1.1)
       )
    } else {
      # initialize the plot
        ymin = min(out_gadget[[input$yaxis[1]]]) * 0.9
        ymax = max(out_gadget[[input$yaxis[1]]]) * 1.1
        for (i in 2:length(input$yaxis)) {
          ymin = min(ymin,min(out_gadget[[input$yaxis[i]]]) * 0.9)
          ymax = max(ymax,max(out_gadget[[input$yaxis[i]]]) * 1.1)
        }

        plot(
          x = if (input$timescale == "sec" && input$xaxis == "time" ) {
            out_gadget[[input$xaxis]] * 60
          } else {
            out_gadget[[input$xaxis]]
          },
          y = out_gadget[[input$yaxis[1]]],
          type = "l",
          col = "skyblue1",
          lwd = 4,
          xlab = if (input$timescale == "sec" && input$xaxis == "time" ) {
            paste0(input$xaxis, " [sec]")
          } else if (input$timescale == "min" && input$xaxis == "time" ) {
            paste0(input$xaxis, " [min]" )
          } else {
            input$xaxis
          },
          ylab = input$yaxis,
          ylim = c(ymin, ymax)
        )

        lapply(2:length(input$yaxis), FUN = function(i) {
          lines(
            x = if (input$timescale == "sec" ) {
              out_gadget[[input$xaxis]] * 60
            } else {
              out_gadget[[input$xaxis]]
            },
            y = out_gadget[[input$yaxis[i]]],
            type = "l",
            col = colors[i],
            lwd = 4,
            ylab = input$yaxis[i],
            ylim = c(ymin, ymax)
          )
        })
        }
    })
  
  # returns table
  output$summary <- DT::renderDataTable( 
    out_gadget(), 
    options = list(
      scrollX = TRUE, 
      pageLength = 10, 
      scrollY = FALSE
    )
  )
  
  # reset all inputs
  observeEvent(input$resetAll, {
    all_inputs <- names(reactiveValuesToList(input))
    lapply(1:length(all_inputs), FUN = function(i) {
      shinyjs::reset(all_inputs[[i]])
    })
    
    # show an alert after 2 seconds
    shinyjs::delay(1000, 
      sendSweetAlert(
        session,
        title = "Yay!",
        text = "All inputs are reset",
        type = "success",
        closeOnClickOutside = TRUE
      )
    )
    
  })

  
  # #reset sliders individually, DOES NOT WORK YET
  # rv <- reactiveValues(lastBtn = character())
  # # get the last selected reset button
  # lapply(
  #   X = 1:length(pars_gadget),
  #   FUN = function(i) {
  #     observeEvent(input[[paste0("reset_", names(pars_gadget)[[i]])]], {
  #       if (input[[paste0("reset_", names(pars_gadget)[[i]])]] > 0) {
  #         rv$lastBtn <- paste0("reset_", names(pars_gadget)[[i]])
  #       }
  #     })
  #   })
  # # use shinyjs to reset buttons
  # observeEvent(extract_slider_input(input, pars_gadget),{
  #   # reset only the slider which corresponds to the selected reset button
  #   slider_name <- unlist(str_split(rv$lastBtn, "_"))[[2]]
  #   shinyjs::reset(paste0(slider_name, "-value"))
  # })
  
  # close the connection when app is shut down
  # onSessionEnded(function() {
  #   close(con)
  # })
  
  # close the gadget
  observeEvent(input$done, {
    stopApp(TRUE)
    # remove old dll files
    if (.Platform$OS.type == "unix") {
      file.remove("temp_core/equations.o")
      file.remove("temp_core/equations.so")
    } else if (.Platform$OS.type == "windows") {
      file.remove("temp_core/equations.dll")
    }
  })

}


shinyApp(ui = ui, server = server)
