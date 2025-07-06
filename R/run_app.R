##*********************************************
##*
##* @file: run_app.R
##*
##* Run ORBayesTtest as a R function
##*
##* Author:
##* Israel Almodovar-Rivera PhD
##* Department of Mathematical Sciences
##* University of Puerto Rico at Mayaguez
##* israel.almodovar@upr.edu
##* Copyright July 2025
##*********************************************


ORBayesTtest <- function() {
  shiny::shinyApp(ui = ui(), server = server)
}