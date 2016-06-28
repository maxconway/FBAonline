
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("FBA online"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      actionButton("reevaluate", "Reevaluate Model"),
      textInput('model_url_1','model_url_1'),
      textInput('model_url_2','model_url_2'),
      selectizeInput('reaction_filter', 'Filter reactions',NULL,multiple=TRUE)
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel(title = 'table',
          dataTableOutput('model_table')
        ),
        tabPanel(title = 'bar comparison',
                 plotOutput('bar_chart')
                 ),
        tabPanel(title = 'heatmap',
                 selectInput('contrast_1', 'contrast_1', 'Group1', 'Group1'),
                 selectInput('contrast_2', 'contrast_2', 'Group2', 'Group2'),
                 plotOutput('heatmap')
                 )
      )
    )
  )
))
