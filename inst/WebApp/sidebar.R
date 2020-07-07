#  This code contains the sidebar of shinydashboard. 

sidebar <- bs4DashSidebar(
  width = 300,
  title = "Human Simulator",
  skin = "light",
  brandColor = NULL,
  url = "http://physiol-seafile.uzh.ch/", 
  src = "logos/heartbeat1.jpeg",
  elevation = 0,
  opacity = 1,
  status = "danger",

  # user panel info
   uiOutput("user_panel"),
  
  # sidebar menu with 3 tabs
  bs4SidebarMenu(
    bs4SidebarMenuItem(
      "App", 
      tabName = "main", 
      icon = "home"
    ),
    bs4SidebarMenuItem(
      "Lexicon", 
      tabName = "lexicon", 
      icon = "book"
    ),
    bs4SidebarMenuItem(
      "About", 
      tabName = "about", 
      icon = "info"
    )
  )
)

