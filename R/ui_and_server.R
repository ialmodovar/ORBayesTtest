
##*********************************************
##*
##* @file: ui_and_server.R
##*
##* UI and server function
##*
##* Author:
##* Israel Almodovar-Rivera PhD
##* Department of Mathematical Sciences
##* University of Puerto Rico at Mayaguez
##* israel.almodovar@upr.edu
##* Copyright July 2025
##*********************************************


ui <- function(){
  fluidPage(
  withMathJax(),
  titlePanel("Objective Bayesian \\(t\\)-test for one and two mean comparison"),
  p(withMathJax("Student's \\(t\\) has been over 100 years since the discovery of one of the most fundamental statistical tests: the Student's \\(t\\) test. In this shiny app, employs the objective and robust Bayesian approach for hypothesis testing for one-sample and two-sample mean comparisons for the assumption of equal variances.")),
  tabsetPanel(
    tabPanel("Data",
             sidebarLayout(
               sidebarPanel(
                 radioButtons("upload_method", "Select Data Input Method:",
                              choices = c("Upload File" = "file", 
                                          "Manual Entry" = "manual")),
                 conditionalPanel(
                   condition = "input.upload_method == 'file'",
                   fileInput("file1", "Choose file", accept = c(".csv", ".xlsx", ".xls", ".ods"))
                 ),
                 conditionalPanel(
                   condition = "input.upload_method == 'manual'",
                   textAreaInput("manual_data", "Input data (including header)", "", rows = 8),
                   downloadButton("download", "Download .csv"),
                   br(),
                   downloadButton("download2", "Download .xlsx/.xls"),
                   br(),
                   downloadButton("download3", "Download .ods")
                 )
               ),
               mainPanel(tableOutput("contents"))
             )
    ),
    tabPanel("One Sample Inference",
             sidebarLayout(
               sidebarPanel(
                 radioButtons(inputId = "one_sample_options",label = "One Sample options", choices=c("Data" = "upload_input_data",
                                                                                                     "Summaries" = "one_summary")),
                 conditionalPanel(condition="input.one_sample_options=='upload_input_data'",
                                  uiOutput("var_select_one_sample_ui")),
                 conditionalPanel(condition="input.one_sample_options=='one_summary'",withMathJax(),
                                  numericInput(inputId = "one.sample.n",label = "Sample size \\( n\\)",value = 1),
                                  numericInput(inputId = "one.sample.mean",label = "Sample mean \\( \\bar{x}\\)",value = 0),
                                  numericInput(inputId = "one.sample.std",label = "Sample standard deviation \\( s\\)",value = 1),
                                  actionButton("run_one_sm", "Run one sample summary")
                                  
                 ),
               ),
               mainPanel(
                 conditionalPanel(
                   condition = "input.one_sample_options == 'upload_input_data'",
                   withMathJax(),
                   uiOutput("UIone_sample")
                 ),
                 
                 conditionalPanel(
                   condition = "input.one_sample_options == 'one_summary'",
                   withMathJax(),
                   uiOutput("UIone_sm")
                 )
               )
             )
    ),
    tabPanel("Two Sample Inference",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("var_select_two_sample_ui"),
                 helpText("Parameters for Gonen (2005) Bayes Factor"),
                 numericInput("lambda", label = withMathJax("\\(\\lambda\\): Effect size under \\( H_1 \\)"), value = 0, step = 0.1),
                 numericInput("sigma2d", label = withMathJax("\\(\\sigma^2_d\\): Prior variance"), value = 1/3, step = 0.1),
                 
                 actionButton("run_two", "Run Two Sample Test")
               ),
               mainPanel(
                 conditionalPanel(condition = "input.two_sample_options == 'two_samples_data'",
                                  withMathJax(),
                                  uiOutput("UItwo_sample")
                 ),
                 conditionalPanel(
                   condition = "input.two_sample_options == 'two_summary'",
                   withMathJax(),
                   uiOutput("UItwo_sm")
                 )
               )
             )
    ),
    tabPanel("Bayes Factor Info",
             withMathJax(),
             h3("What is a Bayes Factor?"),
             p("The Bayes Factor is tool for model comparison and allow us to quantifies how well a dataset support one statistical hypothesis over another."),
             h4("Why Use Bayes Factors?"),
             tags$ul(
               tags$li("Unlike \\( p\\)-values, Bayes Factors  quantify evidence for or against both \\( H_0 \\) and \\( H_1 \\)."),
               tags$li("They incorporate prior knowledge into the comparison.")
             ),
             tags$ul(
               tags$li("Given two hypotheses: \\( H_0 \\): Null hypothesis and \\( H_1:\\) Alternative"
               ),
               p("The Bayes factor \\( B_{01} \\) is defined as:"),
               withMathJax("$$
  B_{01} = \\frac{P(\\text{data} | H_0)}{P(\\text{data} | H_1)}
  $$")),
             tags$li("They relate to posterior probability"),  
             withMathJax("$$
  P(H_0|data) = \\frac{1}{1+\\frac{P(H_0)}{P(H_1)} \\frac{1}{B_{01}} }.
  $$"),
             p("In here \\(P(H_0)\\) and \\(P(H_1)\\) are the priors probabilities of the hypotheses.Also \\( B_{10} = 1/B_{01}. \\)" ),
             br(),
             h3("Interpretation of \\( B_{01} \\)"),
             withMathJax(),
             p("Two commonly used interpretations of Bayes factors are Jeffreys (1961) and Kass and Raftery (1995). Jeffrey's suggested evaluating them on the base-10 logarithm and Kass and Raftery using twice the natural logarithm."),
             fluidRow(
               column(6,
                      h4("Jeffreys Scale"),
                      uiOutput("jeffreys_table")
               ),
               column(6,
                      h4("Kass & Raftery Scale"),
                      uiOutput("kr_table")
               )
             )
    )
    
  )

)
}

server <- function(input,output,session){

# ---- Data input ----
mydata <- reactiveVal(NULL)

observeEvent({
  input$file1
  input$manual_data
  input$upload_method
}, {
  req(input$upload_method)
  
  if (input$upload_method == "file" && !is.null(input$file1)) {
    ext <- tools::file_ext(input$file1$name)
    df <- switch(ext,
                 csv = read.csv(input$file1$datapath),
                 xlsx = read_excel(input$file1$datapath),
                 xls  = read_excel(input$file1$datapath),
                 ods  = read_ods(input$file1$datapath),
                 NULL)
    mydata(df)
  } else if (input$upload_method == "manual") {
    tryCatch({
      df <- read.csv(text = input$manual_data)
      mydata(df)
    }, error = function(e) {
      showNotification("Error reading manual input.", type = "error")
    })
  }
})

output$contents <- renderTable({
  head(mydata(), 15)
})

output$download <- downloadHandler(
  filename = function() "manual_input.csv",
  content = function(file) {
    write.csv(mydata(), file, row.names = FALSE)
  }
)

output$download2 <- downloadHandler(
  filename = function() "manual_input.xlsx",
  content = function(file) {
    wb <- createWorkbook()
    addWorksheet(wb, "Sheet 1")
    writeData(wb, sheet = "Sheet 1", mydata(), rowNames = FALSE)
    saveWorkbook(wb, file, overwrite = TRUE)
  }
)

output$download3 <- downloadHandler(
  filename = function() "manual_input.ods",
  content = function(file) {
    write_ods(mydata(), file, row_names = FALSE)
  }
)

# ---- One sample inference ----

output$var_select_one_sample_ui <- renderUI({
  df <- mydata()
  req(df)
  num_vars <- names(df)[sapply(df, is.numeric)]
  selectInput("num.vars", "Select numerical variables:", choices = num_vars, multiple = TRUE)
})

output$UIone_sample <- renderUI({
  req(input$num.vars)
  df <- mydata()
  var_list <- input$num.vars
  
  lapply(var_list, function(varname) {
    safe_varname <- gsub("[^a-zA-Z0-9]", "_", varname)
    t_output_id <- paste0("mean_", safe_varname)
    mu_input_id <- paste0("mu0_", safe_varname)
    p0_id <- paste0("p0_", safe_varname)
    stats_id <- paste0("stats_", safe_varname)
    
    output[[t_output_id]] <- renderTable({
      req(df[[varname]], input[[mu_input_id]], input[[p0_id]])
      mu0 <- as.numeric(input[[mu_input_id]])
      p0 <- as.numeric(input[[p0_id]])
      x <- df[[varname]]
      one.sample.t.test.bf(x = x, mu0 = mu0, p0 = p0)
    }, rownames = TRUE, digits = 5)
    
    output[[stats_id]] <- renderUI({
      x <- df[[varname]]
      m <- mean(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      n <- sum(!is.na(x))
      
      withMathJax(HTML(
        paste0(
          "<table class='table'>",
          "<tr><td>\\( \\bar{x} = \\)</td><td>", round(m, 3), "</td></tr>",
          "<tr><td>\\( s = \\)</td><td>", round(s, 3), "</td></tr>",
          "<tr><td>\\( n = \\)</td><td>", n, "</td></tr>",
          "</table>"
        )
      ))
    })
    
    tagList(
      tags$h4(paste("Variable:", varname)),
      fluidRow(
        column(6, textInput(mu_input_id, label = withMathJax(paste0("\\( H_0: \\mu = \\mu_0 \\) for ", varname)), value = "0")),
        column(6, sliderInput(p0_id, label = withMathJax("Prior probability \\(P(H_0)\\)"), min = 0, max = 1, value = 0.5, step = 0.01))
      ),
      uiOutput(stats_id),
      tableOutput(t_output_id),
      tags$hr()
    )
  })
})

##******
##* summary one sample
##*******

output$UIone_sm <- renderUI({
  tagList(
    fluidRow(
      column(6, numericInput("mu0_sm", label = withMathJax("\\( H_0: \\mu = \\mu_0 \\)"), value = 0)),
      column(6, sliderInput("p0_sm", label = withMathJax("Prior probability \\(P(H_0)\\)"), min = 0, max = 1, value = 0.5, step = 0.01))
    ),
    tableOutput("one_summary_result"),
    tags$hr()
  )
})

output$one_summary_result <- renderTable({
  req(input$one.sample.mean, input$one.sample.std, input$one.sample.n, input$mu0_sm, input$p0_sm)
  
  one.sample.t.test.bf.sm(
    xbar = input$one.sample.mean,
    sx = input$one.sample.std,
    n = input$one.sample.n,
    mu0 = input$mu0_sm,
    p0 = input$p0_sm
  )
}, rownames = TRUE, digits = 5)

# ---- Two sample inference ----

output$var_select_two_sample_ui <- renderUI({
  tagList(
    radioButtons("two_sample_options", "Select option:",
                 choices = c("Data Input" = "two_samples_data", "Summaries" = "two_summary"),
                 selected = "two_samples_data"),
    
    conditionalPanel(condition = "input.two_sample_options == 'two_samples_data'",
                     
                     radioButtons("two_sample_mode", "Select input mode:",
                                  choices = c("Two separate numeric variables" = "vars",
                                              "One numeric variable + grouping variable" = "group"),
                                  selected = "vars"),
                     checkboxInput("paired", "Paired samples", FALSE),
                     
                     uiOutput("two_sample_mode_ui")  # this sub-UI depends on df()
    ),
    
    conditionalPanel(condition = "input.two_sample_options == 'two_summary'",
                     withMathJax(),
                     br(),
                     h4("Group 1"),
                     numericInput(inputId = "two.sample.n1", label = "Sample size \\( n_1 \\)", value = 2),
                     numericInput(inputId = "two.sample.mean1", label = "Sample mean \\( \\bar{x}_1 \\)", value = 0),
                     numericInput(inputId = "two.sample.std1", label = "Sample standard deviation \\( s_1 \\)", value = 1),
                     
                     h4("Group 2"),
                     numericInput(inputId = "two.sample.n2", label = "Sample size \\( n_2 \\)", value = 2),
                     numericInput(inputId = "two.sample.mean2", label = "Sample mean \\( \\bar{x}_2 \\)", value = 0),
                     numericInput(inputId = "two.sample.std2", label = "Sample standard deviation \\( s_2 \\)", value = 1),
                    )
  )
})

output$two_sample_mode_ui <- renderUI({
  req(mydata())
  df <- mydata()
  num_vars <- names(df)[sapply(df, is.numeric)]
  group_vars <- names(df)[sapply(df, function(x) is.factor(x) || is.character(x))]
  
  if (input$two_sample_mode == "vars") {
    tagList(
      selectInput("var1", "First Numeric Variable:", choices = num_vars),
      selectInput("var2", "Second Numeric Variable:", choices = num_vars)
    )
  } else if (input$two_sample_mode == "group") {
    tagList(
      selectInput("numvar_grouped", "Select Numeric Variable:", choices = num_vars),
      selectInput("groupvar", "Select Grouping Variable (2 levels):", choices = group_vars)
    )
  }
})

df <- reactive({ req(mydata()); mydata() })

two.vars <- eventReactive(input$run_two, {
  data <- df()
  mode <- req(input$two_sample_mode)
  
  if (mode == "group") {
    req(input$numvar_grouped, input$groupvar)
    data <- na.omit(data[, c(input$numvar_grouped, input$groupvar)])
    data[[input$groupvar]] <- as.factor(data[[input$groupvar]])
    data[[input$groupvar]] <- droplevels(data[[input$groupvar]])
    levels_group <- levels(data[[input$groupvar]])
    validate(need(length(levels_group) == 2, "Grouping variable must have exactly 2 levels"))
    
    x <- data[data[[input$groupvar]] == levels_group[1], input$numvar_grouped]
    y <- data[data[[input$groupvar]] == levels_group[2], input$numvar_grouped]
  } else {
    req(input$var1, input$var2)
    x <- data[[input$var1]]
    y <- data[[input$var2]]
    levels_group <- c("Group 1", "Group 2")
  }
  
  list(x = as.numeric(x), y = as.numeric(y), group_labels = levels_group)
})


output$UItwo_sample <- renderUI({
  req(two.vars())
  x <- two.vars()$x
  y <- two.vars()$y
  
  tagList(
    fluidRow(
      column(6, textInput("mu0_two", label = withMathJax("\\( H_0: \\mu_1 - \\mu_2 =  \\)"), value = "0")),
      column(6, sliderInput("p0_two", label = withMathJax("Prior probability \\(P(H_0)\\)"), min = 0, max = 1, value = 0.5, step = 0.01))
    ),
    
    if (!input$paired) {
      m1 <- mean(x, na.rm = TRUE)
      m2 <- mean(y, na.rm = TRUE)
      s1 <- sd(x, na.rm = TRUE)
      s2 <- sd(y, na.rm = TRUE)
      n1 <- sum(!is.na(x))
      n2 <- sum(!is.na(y))
      
      withMathJax(HTML(
        paste0(
          "<table class='table'>",
          "<tr><th></th><th>Group 1</th><th>Group 2</th></tr>",
          "<tr><td>\\( \\bar{x} \\)</td><td>", round(m1, 3), "</td><td>", round(m2, 3), "</td></tr>",
          "<tr><td>\\( s \\)</td><td>", round(s1, 3), "</td><td>", round(s2, 3), "</td></tr>",
          "<tr><td>\\( n \\)</td><td>", n1, "</td><td>", n2, "</td></tr>",
          "</table>"
        )
      ))
    } else {
      d <- x - y
      m <- mean(d, na.rm = TRUE)
      s <- sd(d, na.rm = TRUE)
      n <- sum(!is.na(d))
      
      withMathJax(HTML(
        paste0(
          "<table class='table'>",
          "<tr><td>\\( \\bar{d} = \\)</td><td>", round(m, 3), "</td></tr>",
          "<tr><td>\\( s_d = \\)</td><td>", round(s, 3), "</td></tr>",
          "<tr><td>\\( n = \\)</td><td>", n, "</td></tr>",
          "</table>"
        )
      ))
    },
    
    tableOutput("mean")
  )
})

output$mean <- renderTable({
  req(two.vars())
  two.sample.t.test.bf(two.vars()$x, two.vars()$y,
                       mu0 = as.numeric(input$mu0_two),
                       p0 = as.numeric(input$p0_two),
                       paired = input$paired,
                       lambda = input$lambda,
                       sigma2d = input$sigma2d)
}, rownames = TRUE, digits = 5)


output$UItwo_sm <- renderUI({
  tagList(
    fluidRow(
      column(6, numericInput("mu0_sm2", label = withMathJax("\\( H_0: \\mu_1-\\mu_2 =  \\)"), value = 0)),
      column(6, sliderInput("p0_sm2", label = withMathJax("Prior probability \\(P(H_0)\\)"), min = 0, max = 1, value = 0.5, step = 0.01))
    ),
    tableOutput("two_summary_result"),
    tags$hr()
  )
})

output$two_summary_result <- renderTable({
  req(input$two.sample.mean1, input$two.sample.std1, input$two.sample.n1,input$two.sample.mean2, input$two.sample.std2, input$two.sample.n2, input$mu0_sm2, input$p0_sm2)
  
  two.sample.t.test.bf.sm(
    xbar1 = input$two.sample.mean1,
    s21 = (input$two.sample.std1^2),
    n1 = input$two.sample.n1,
    xbar2 = input$two.sample.mean2,
    s22 = (input$two.sample.std2^2),
    n2 = input$two.sample.n2,
    mu0 = input$mu0_sm2,
    p0 = input$p0_sm2)
}, rownames = TRUE, digits = 5)



output$jeffreys_table <- renderUI({
  HTML('
    <table class="table table-striped table-sm">
      <thead>
        <tr>
          <th>\\( B_{01} \\)</th>
          <th>\\( \\log_{10}(B_{01}) \\)</th>
          <th>Evidence</th>
        </tr>
      </thead>
      <tbody>
        <tr><td>\\( [0, 1) \\)</td><td>\\( (-\\infty, 0) \\)</td><td>Negative (Evidence favors \\( H_1 \\))</td></tr>
        <tr><td>\\( [1, 3) \\)</td><td>\\( [0, 0.48) \\)</td><td>Weak evidence for \\( H_0 \\)</td></tr>
        <tr><td>\\( [3, 10) \\)</td><td>\\( [0.48, 1) \\)</td><td>Moderate evidence for \\( H_0 \\)</td></tr>
        <tr><td>\\( [10, 30) \\)</td><td>\\( [1, 1.48) \\)</td><td>Strong evidence for \\( H_0 \\)</td></tr>
        <tr><td>\\( [30, 100) \\)</td><td>\\( [1.48, 2) \\)</td><td>Very strong evidence for \\( H_0 \\)</td></tr>
        <tr><td>\\( [100, \\infty) \\)</td><td>\\( [2, \\infty) \\)</td><td>Decisive evidence for \\( H_0 \\)</td></tr>
      </tbody>
    </table>
  ')
})

output$kr_table <- renderUI({
  HTML('
    <table class="table table-striped table-sm">
      <thead>
        <tr>
          <th>\\( B_{01} \\)</th>
          <th>\\( 2 \\log(B_{01}) \\)</th>
          <th>Evidence</th>
        </tr>
      </thead>
      <tbody>
        <tr><td>\\( [0, 1) \\)</td><td>\\( (-\\infty, 0) \\)</td><td>Negative (Evidence favors \\( H_1 \\))</td></tr>
        <tr><td>\\( [1, 3) \\)</td><td>\\( [0, 2.20) \\)</td><td>Weak evidence for \\( H_0 \\)</td></tr>
        <tr><td>\\( [3, 20) \\)</td><td>\\( [2.20, 6) \\)</td><td>Positive (Moderate evidence for \\( H_0 \\))</td></tr>
        <tr><td>\\( [20, 150) \\)</td><td>\\( [6, 10) \\)</td><td>Strong evidence for \\( H_0 \\)</td></tr>
        <tr><td>\\( [150, \\infty) \\)</td><td>\\( [10, \\infty) \\)</td><td>Very strong evidence for \\( H_0 \\)</td></tr>
      </tbody>
    </table>
  ')
})

}




