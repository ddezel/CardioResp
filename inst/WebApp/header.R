#-------------------------------------------------------------------------
#  This code contains the header of shinydashboard. It is modified compared
#  to classic header. Indeed, some buttons to save, load, reset, download are
#  inserted in the header bar. Moreover, users can change the global theme
#  clicking on the theme selector.
#
#  David Granjon, the Interface Group, Zurich
#  December 4th, 2017
#-------------------------------------------------------------------------

navbar <- bs4DashNavbar(
  skin = "light", 
  status = "gray-light", 
  border = FALSE,
  sidebarIcon = "bars", 
  controlbarIcon = NULL, 
  leftUi = NULL,
  rightUi = NULL, 
  fixed = FALSE,
  id = "navbar"
)