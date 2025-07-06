library(shiny)
library(ORBayesTtest)

shinyApp(ui = ORBayesTtest::ui(), server = ORBayesTtest::server)