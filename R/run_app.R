##*********************************************
##*
##* @file: run_app.R
##*
##* Run ORBayesTtest.App to run Shiny App
##*
##* Author:
##* Israel Almodovar-Rivera PhD
##* Department of Mathematical Sciences
##* University of Puerto Rico at Mayaguez
##* israel.almodovar@upr.edu
##* Copyright July 2025
##*********************************************


ORBayesianTtest.app <- function() {
  shinyApp(ui = ui(), server = server)
}