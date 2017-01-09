
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
      actionButton("reevaluate", "Evaluate Model"),
      textInput('model_url_1','model url 1'),
      textInput('model_url_2','model url 2'),
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
        tabPanel(title = 'table',
          dataTableOutput('model_table')
        ),
        tabPanel(title = 'bar comparison',
                 plotOutput('bar_chart')
                 ),
        # tabPanel(title = 'heatmap',
        #          selectInput('contrast_1', 'contrast_1', 'Group1', 'Group1'),
        #          selectInput('contrast_2', 'contrast_2', 'Group2', 'Group2'),
        #          plotOutput('heatmap')
        #          ),
        # tabPanel(title = 'metabolites',
        #          dataTableOutput('metabolite_table')
        #          ),
        tabPanel(title = 'settings',
                 textInput(inputId = 'pattern_arrow', 
                           label = 'Regex for arrow in equations',
                           value = '<?[-=]+>')
                 )
      )
    )
  )
))
