#----------------------------------#
# Simulation Box helper functions #
#----------------------------------#

library(plotly)
library(magrittr)
library(data.table)

# define plot axis -----------------------------------------------------------
#'
#' @export
c(
  yaxisSymp <- list(
                gridcolor = 'rgb(255,255,255)',
                zeroline = FALSE
  ),
  yaxisResist <- list(title = "[mmHg·min/mL]", 
                      gridcolor = 'rgb(255,255,255)',
                      # showline = FALSE,
                      # showgrid = FALSE,
                       zeroline = FALSE
  ),
  yaxisFlow <- list(title = "[mL/min]", 
                    gridcolor = 'rgb(255,255,255)',
                    zeroline = FALSE
  ),
  yaxisPress <- list(title = "[mmHg]", 
                     gridcolor = 'rgb(255,255,255)',
                     zeroline = FALSE
  ),
  yaxisRate <- list(title = "[beats/min]", 
                    gridcolor = 'rgb(255,255,255)',
                    zeroline = FALSE
  ),
  
  xaxis_iso <- list(title = "time [sec]", zeroline = FALSE, gridcolor = 'rgb(255,255,255)',
                titlefont = list(size = 12), range = c(1000, 6000)),
  
  xaxis <- list(title = "time [sec]", zeroline = FALSE, gridcolor = 'rgb(255,255,255)',
                    titlefont = list(size = 12)),
  
  # define plot annotations ----------------------------------------------------
  # annotations are necessary to keep the plot titles in the gathered subplots()
  aResist <- list(
    text = "Resistances",
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 1,
    showarrow = FALSE
  ),
  aSN <- list(
    text = "Sympathetic Activity",
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 1,
    showarrow = FALSE
  ),
  
  aResist <- list(
    text = "Resistances",
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 1,
    showarrow = FALSE
  ),
  aSN <- list(
    text = "Sympathetic Activity",
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 1,
    showarrow = FALSE
  ),
  aCO <- list(
    text = "Cardiac Output",
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 1,
    showarrow = FALSE
  ),
  aFlows <- list(
    text = "Flows",
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 1,
    showarrow = FALSE
  ),
  aMAP <- list(
    text = "Mean Arterial Pressure",
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 1,
    showarrow = FALSE
  ),
  aHR <- list(
    text = "Heart Rate",
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 1,
    showarrow = FALSE
  )
)

#------------------------------------------------------------------------- 
# PLOTS FOR ISOCAPNIC STATE  
#-------------------------------------------------------------------------
# path to isocapnia csv file
iso_table <- fread("csvfiles/isocapnia.csv", header = T, sep = ",") # freads is WAY faster than read.csv

#' @export
resistances_plot <- plot_ly(
  iso_table, x = iso_table$time*60) %>%
  add_lines(y = iso_table$Rtot,
            name = "Total Resistance", line = list(color = "darkred", width = 1,
                                                   dash = "dot")) %>%
  add_lines(y = iso_table$Rs, ymin = 0.5 * min(iso_table$Rbrain), ymax = 1.5 * max(iso_table$Rbrain),
            name = "Systemic Resistance", line = list(color = "darkred", width = 1,
                                                      dash = "dash")) %>%
  add_lines(y = iso_table$Rbrain,
            name = "Cerebral Resistance", line = list(color = "darkred", width = 1)) %>%
  layout(paper_bgcolor='rgb(255,255,255)', annotations = aResist, 
         xaxis = xaxis_iso, yaxis = yaxisResist)

#' @export
symp_act_plot <- plotly::plot_ly(
  iso_table, x = iso_table$time*60) %>%
  plotly::add_lines(y = iso_table$SN_activity, ymin = 0.5 * min(iso_table$SN_activity), ymax = 1.5 * max(iso_table$SN_activity),
                    name = "Sympathetic Activity", line = list(color = "darkred", width = 1)) %>%
  plotly::layout(annotations = aSN, paper_bgcolor='rgb(255,255,255)', showlegend = FALSE,
                 xaxis = xaxis_iso, yaxis = yaxisSymp)


# plot ISOCAPNIC cardiac output and the aortic and cerebral flows 
#' @export
  CO_plot <- plotly::plot_ly(
    iso_table, x = iso_table$time*60) %>%
    plotly::add_lines(y = iso_table$CO, ymin = 0.5 * min(iso_table$CO), ymax = 1.5 * max(iso_table$CO),
                      name = "Cardiac Output", line = list(color = "darkred", width = 1), 
    ) %>%
    plotly::layout( paper_bgcolor='rgb(255,255,255)', xaxis = xaxis_iso, yaxis = yaxisFlow,
                    annotations = aCO)
  
  #' @export
  flow_plot <- plotly::plot_ly(
    iso_table, x = iso_table$time*60) %>%
    plotly::add_lines(y = iso_table$Qaa, ymin = 0.5 * min(iso_table$Qaa), ymax = 1.5 * max(iso_table$Qaa),
                      name = "Arterial Flow", line = list(color = "darkred", width = 1), 
    ) %>%
    plotly::add_lines(y = iso_table$Qc, 
                      name = "Cerebral Flow", line = list(color = "red", width = 1), 
    ) %>%
    plotly::layout(annotations = aFlows, paper_bgcolor='rgb(255,255,255)', 
                   xaxis = xaxis_iso, yaxis = yaxisFlow)


# plot ISOCAPNIC MAP and the heart rate
  #' @export
  MAP_plot <- plotly::plot_ly(
    iso_table, x = iso_table$time*60) %>%
    plotly::add_lines(y = iso_table$MAP, ymin = 0.5 * min(iso_table$MAP), ymax = 1.5 * max(iso_table$MAP),
                      name = "MAP", line = list(color = "darkred", width = 1), ) %>%
    plotly::layout(annotations = aMAP, paper_bgcolor='rgb(255,255,255)', 
                   xaxis = xaxis_iso, yaxis = yaxisPress)
  
  #' @export
  HR_plot <- plotly::plot_ly(
    iso_table, x = iso_table$time*60) %>%
    plotly::add_lines(y = iso_table$HR, ymin = 0.5 * min(iso_table$HR), ymax = 1.5 * max(iso_table$HR),
                      name = "Heart Rate", line = list(color = "darkred", width = 1), 
    ) %>%
    plotly::layout(annotations = aHR, paper_bgcolor='rgb(255,255,255)', 
                   xaxis = xaxis_iso, yaxis = yaxisRate)
  

#------------------------------------------------------------------------- 
# PLOTS FOR HYPERCAPNIC STATE  
#-------------------------------------------------------------------------

# path to hypercapnia csv file
hyper_table <- fread("csvfiles/hypercapnia.csv", header = T, sep = ",") # freads is WAY faster than read.csv

#' @export
  resistances_plot_hyper <- plot_ly(
    hyper_table, x = hyper_table$time*60) %>%
    add_lines(y = hyper_table$Rtot,
              name = "Total Resistance", line = list(color = "darkred", width = 1,
                                                     dash = "dot")) %>%
    add_lines(y = hyper_table$Rs, ymin = 0.5 * min(hyper_table$Rbrain), ymax = 1.5 * max(hyper_table$Rbrain),
              name = "Systemic Resistance", line = list(color = "darkred", width = 1,
                                                        dash = "dash")) %>%
    add_lines(y = hyper_table$Rc,
              name = "Cerebral Resistance", line = list(color = "darkred", width = 1)) %>%
    layout(paper_bgcolor='rgb(255,255,255)', annotations = aResist,
           xaxis = xaxis, yaxis = yaxisResist)

  #' @export
  symp_act_plot_hyper <- plotly::plot_ly(
    hyper_table, x = hyper_table$time*60) %>%
    plotly::add_lines(y = hyper_table$SN_activity, ymin = 0.5 * min(hyper_table$SN_activity), ymax = 1.5 * max(hyper_table$SN_activity),
                      name = "Sympathetic Activity", line = list(color = "darkred", width = 1)) %>%
    plotly::layout(annotations = aSN, paper_bgcolor='rgb(255,255,255)', showlegend = FALSE,
                   xaxis = xaxis, yaxis = yaxisSymp)

#' @export
  CO_plot_hyper <- plotly::plot_ly(
    hyper_table, x = hyper_table$time*60) %>%
    plotly::add_lines(y = hyper_table$CO, ymin = 0.5 * min(hyper_table$CO), ymax = 1.5 * max(hyper_table$CO),
                      name = "Cardiac Output", line = list(color = "darkred", width = 1), 
    ) %>%
    plotly::layout( paper_bgcolor='rgb(255,255,255)', xaxis = xaxis, yaxis = yaxisFlow,
                    annotations = aCO)

#' @export
  flow_plot_hyper <- plotly::plot_ly(
    hyper_table, x = hyper_table$time*60) %>%
    plotly::add_lines(y = hyper_table$Qaa, ymin = 0.5 * min(hyper_table$Qaa), ymax = 1.5 * max(hyper_table$Qaa),
                      name = "Arterial Flow", line = list(color = "darkred", width = 1), 
    ) %>%
    plotly::add_lines(y = hyper_table$Qc, 
                      name = "Cerebral Flow", line = list(color = "red", width = 1), 
    ) %>%
    plotly::layout(annotations = aFlows, paper_bgcolor='rgb(255,255,255)', 
                   xaxis = xaxis, yaxis = yaxisFlow)

#' @export
  MAP_plot_hyper <- plotly::plot_ly(
    hyper_table, x = hyper_table$time*60) %>%
    plotly::add_lines(y = hyper_table$MAP, ymin = 0.5 * min(hyper_table$MAP), ymax = 1.5 * max(hyper_table$MAP),
                      name = "MAP", line = list(color = "darkred", width = 1), ) %>%
    plotly::layout(annotations = aMAP, paper_bgcolor='rgb(255,255,255)', 
                   xaxis = xaxis, yaxis = yaxisPress)
#' @export
  HR_plot_hyper <- plotly::plot_ly(
    hyper_table, x = hyper_table$time*60) %>%
    plotly::add_lines(y = hyper_table$HR, ymin = 0.5 * min(hyper_table$HR), ymax = 1.5 * max(hyper_table$HR),
                      name = "Heart Rate", line = list(color = "darkred", width = 1), 
    ) %>%
    plotly::layout(annotations = aHR, paper_bgcolor='rgb(255,255,255)', 
                   xaxis = xaxis, yaxis = yaxisRate)
  
#------------------------------------------------------------------------- 
# PLOTS FOR HYPOCAPNIC STATE  
#-------------------------------------------------------------------------

# path to hypocapnia csv file
  hypo_table <- fread("csvfiles/hypocapnia.csv", header = T, sep = ",")

#' @export
  resistances_plot_hypo <- plot_ly(
    hypo_table, x = hypo_table$time*60) %>%
    add_lines(y = hypo_table$Rtot,
              name = "Total Resistance", line = list(color = "darkred", width = 1,
                                                     dash = "dot")) %>%
    add_lines(y = hypo_table$Rs, ymin = 0.5 * min(hypo_table$Rc), ymax = 1.5 * max(hypo_table$Rc),
              name = "Systemic Resistance", line = list(color = "darkred", width = 1,
                                                        dash = "dash")) %>%
    add_lines(y = hypo_table$Rc,
              name = "Cerebral Resistance", line = list(color = "darkred", width = 1)) %>%
    layout(paper_bgcolor='rgb(255,255,255)', annotations = aResist, 
           xaxis = xaxis, yaxis = yaxisResist)
#' @export
  symp_act_plot_hypo <- plotly::plot_ly(
    hypo_table, x = hypo_table$time*60) %>%
    plotly::add_lines(y = hypo_table$SN_activity, ymin = 0.5 * min(hypo_table$SN_activity), ymax = 1.5 * max(hypo_table$SN_activity),
                      name = "Sympathetic Activity", line = list(color = "darkred", width = 1)) %>%
    plotly::layout(annotations = aSN, paper_bgcolor='rgb(255,255,255)', showlegend = FALSE,
                   xaxis = xaxis, yaxis = yaxisSymp)
  
 #' @export
  CO_plot_hypo <- plotly::plot_ly(
    hypo_table, x = hypo_table$time*60) %>%
    plotly::add_lines(y = hypo_table$CO, ymin = 0.5 * min(hypo_table$CO), ymax = 1.5 * max(hypo_table$CO),
                      name = "Cardiac Output", line = list(color = "darkred", width = 1), 
    ) %>%
    plotly::layout( paper_bgcolor='rgb(255,255,255)', xaxis = xaxis, yaxis = yaxisFlow,
                    annotations = aCO)
  #' @export
  flow_plot_hypo <- plotly::plot_ly(
    hypo_table, x = hypo_table$time*60) %>%
    plotly::add_lines(y = hypo_table$Qaa, ymin = 0.5 * min(hypo_table$Qaa), ymax = 1.5 * max(hypo_table$Qaa),
                      name = "Arterial Flow", line = list(color = "darkred", width = 1), 
    ) %>%
    plotly::add_lines(y = hypo_table$Qc, 
                      name = "Cerebral Flow", line = list(color = "red", width = 1), 
    ) %>%
    plotly::layout(annotations = aFlows, paper_bgcolor='rgb(255,255,255)', 
                   xaxis = xaxis, yaxis = yaxisFlow)

#' @export
MAP_plot_hypo <- plotly::plot_ly(
    hypo_table, x = hypo_table$time*60) %>%
    plotly::add_lines(y = hypo_table$MAP, ymin = 0.5 * min(hypo_table$MAP), ymax = 1.5 * max(hypo_table$MAP),
                      name = "MAP", line = list(color = "darkred", width = 1), ) %>%
    plotly::layout(annotations = aMAP, paper_bgcolor='rgb(255,255,255)', 
                   xaxis = xaxis, yaxis = yaxisPress)
#' @export
  HR_plot_hypo <- plotly::plot_ly(
    hypo_table, x = hypo_table$time*60) %>%
    plotly::add_lines(y = hypo_table$HR, ymin = 0.5 * min(hypo_table$HR), ymax = 1.5 * max(hypo_table$HR),
                      name = "Heart Rate", line = list(color = "darkred", width = 1), 
    ) %>%
    plotly::layout(annotations = aHR, paper_bgcolor='rgb(255,255,255)', 
                   xaxis = xaxis, yaxis = yaxisRate)
  
