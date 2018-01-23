
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(visNetwork)

shinyUI(fluidPage(

  # Application title
  titlePanel("FBA online"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      actionButton("reevaluate", "Evaluate Model"),
      textInput('model_url_1','Model URL'),
      selectizeInput('reaction_filter', 'Filter reactions',NULL,multiple=TRUE),
      h3('Instructions'),
      p('Paste a link to a google sheet containing a metabolic model into one of the "model url" boxes.'),
      p('An example of a suitable sheet for an E. Coli model is here: ', a('iJO1366', href='https://docs.google.com/spreadsheets/d/1XdpAKFyEpGjPmI3UZYYK4zQw-RUzKWD2_GJnp3LO_Pk/edit?usp=sharing')),
      p('To compare models, enter both into the "model url" boxes.'),
      p('Use the tabs to view different visualizations of the models'),
      p('You can filter to certain reactions using the Filter reactions box')
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel(title = 'Flux Table',
          dataTableOutput('model_table')
        ),
        tabPanel(title = 'Bar Comparison',
                 inputPanel(textInput('model_url_2','URL for Model 2')),
                 plotOutput('bar_chart')
                 ),
        tabPanel(title = 'Network Visualization',
                 p('Select a model and filter some reactions to see their network.'),
                 visNetworkOutput("network")
                 ),
        tabPanel(title = 'Metabolite Table',
                 dataTableOutput('metabolite_table')
                ),
        tabPanel(title = 'Reaction Clusters',
                 inputPanel(sliderInput('cut_h', 'Cut Height', min=0, max=1000, value=100),
                            sliderInput('min_members', 'Minimum Members', min=1, max=50, value=4, round = TRUE)),
                 plotOutput('cluster_chart')
        )
        # tabPanel(title = 'heatmap',
        #          selectInput('contrast_1', 'contrast_1', 'Group1', 'Group1'),
        #          selectInput('contrast_2', 'contrast_2', 'Group2', 'Group2'),
        #          plotOutput('heatmap')
        #          ),
      )
    )
  )
))
