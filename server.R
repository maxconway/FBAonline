
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(gsheet)
library(dplyr)
library(purrr)
library(logging)
library(magrittr)
library(ggplot2)
library(tidyr)

basicConfig('FINEST')
addHandler(writeToConsole)

shinyServer(function(input, output, session) {
  
  model1 <- reactive({
    logfine('Started evaluation: model1')
    input$reevaluate
    model_url_1 <- input$model_url_1
    if(nchar(model_url_1)<10){
      loginfo('url too short')
      data_frame()
    }else{
      gsheet::gsheet2tbl(model_url_1)
    }
  })
  
  model1_parsed <- reactive({
    logfine('Started evaluation: model1_parsed')
    model1 <- model1()
    if(!is.null(model1) & nrow(model1) > 0){
      FluxBalanceAnalyzeR::parse_reaction_table(model1)
    }else{
      'model missing'
    }
  })
  
  model1_evaluated <- reactive({
    logfine('Started evaluation: model1_evaluated')
    model_parsed <- model1_parsed()
    if(!is.null(model_parsed) & is.list(model_parsed)){
      gurobi::gurobi(model_parsed)
    }else{
      NULL
    }
  })
  
  model2 <- reactive({
    logfine('Started evaluation: model2')
    input$reevaluate
    model_url_2 <- input$model_url_2
    if(nchar(model_url_2)<10){
      loginfo('url too short')
      data_frame()
    }else{
      gsheet::gsheet2tbl(model_url_2)
    }
  })
  
  model2_parsed <- reactive({
    logfine('Started evaluation: model2_parsed')
    model2 <- model2()
    if(!is.null(model2) & nrow(model2) > 0){
      FluxBalanceAnalyzeR::parse_reaction_table(model2)
    }else{
      'model missing'
    }
  })
  
  model2_evaluated <- reactive({
    logfine('Started evaluation: model1_evaluated')
    model_parsed <- model2_parsed()
    if(!is.null(model_parsed) & is.list(model_parsed)){
      gurobi::gurobi(model_parsed)
    }else{
      NULL
    }
  })
  
  full_list_results <- reactive({
    logfine('Started evaluation: full_list_results')
    model1 <- model1()
    model1_parsed <- model1_parsed()
    model1_evaluated <- model1_evaluated()
    model2 <- model2()
    model2_parsed <- model2_parsed()
    model2_evaluated <- model2_evaluated()
    
    logfine('Retrieved all reactives: full_list_results')
    
    data_frame(group = c('Group1', 'Group2'),
               fba_mod = list(model1, model2),
               gurobi_mod = list(model1_parsed, model2_parsed),
               gurobi_result = list(model1_evaluated, model2_evaluated)
    )%>%
      filter(!map_lgl(group, is_null), 
             !map_lgl(fba_mod, is_null), 
             !map_lgl(gurobi_result, is_null)
      ) %>% 
      mutate(fluxmat = map2(gurobi_mod, gurobi_result, function(mod, res){
        mod$A * matrix(res$x, nrow=nrow(mod$A), ncol=ncol(mod$A), byrow=TRUE)
      }))
  })
  
  reaction_metab_results <- reactive({
    logfine('Started evaluation: reaction_metab_results')
    full_list_results <- full_list_results()
    
    if(nrow(full_list_results)>0){  
      full_list_results %>%
        select(group, fluxmat) %>%
        mutate(fluxdf = map(fluxmat, function(x){
          mat <- as(x, 'TsparseMatrix')
          data_frame(reaction = colnames(x)[mat@j+1], metabolite = rownames(x)[mat@i+1], flux = mat@x)
        })) %>%
        select(-fluxmat) %>%
        unnest(fluxdf)
    }else{
      data_frame()
    }
  })
  
  ##########################################################################################################################
  
  observe({
    logfine('Started evaluation: observe')
    reaction_metab_results <- reaction_metab_results()
    
    if(is.data.frame(reaction_metab_results) && ('reaction' %in% names(reaction_metab_results))){
      updateSelectizeInput(session, 'reaction_filter', choices = unique(reaction_metab_results$reaction))
    }
  })
  
  
  ##########################################################################################################################
  
  filtered_reactions_with_fluxes <- reactive({
    logfine('Started evaluation: filtered_reactions_with_fluxes')
    full_list_results <- full_list_results()
    reaction_filter <- input$reaction_filter
    logfine('Retrieved all reactives: filtered_reactions_with_fluxes')
    if(nrow(full_list_results)==0){
      data_frame(abbreviation = character(0),
                 flux = double(0),
                 group=character(0))
    }else{
      full_list_results %>%
        select(-fluxmat, -gurobi_mod) %>%
        mutate(fluxes = map(gurobi_result, function(y){
          data_frame(flux=y$x)
        })) %>%
        unnest(fba_mod, fluxes) %>%
        filter(abbreviation %in% reaction_filter | length(reaction_filter)==0)
    }
  })
  
  
  
  
  #########################################################################################################################
  
  output$model_table <- renderDataTable({
    logfine('Started evaluation: model_table')
    filtered_reactions_with_fluxes <- filtered_reactions_with_fluxes()
    logfine('Retrieved all reactives: model_table')
    filtered_reactions_with_fluxes %>%
      spread(group, flux)
  })
  
  output$bar_chart <- renderPlot({
    logfine('Started evaluation: bar_chart')
    filtered_reactions_with_fluxes <- filtered_reactions_with_fluxes()
    logfine('Retrieved all reactives: bar_chart')
    filtered_reactions_with_fluxes %>%
      ggplot(aes(x=abbreviation, y=flux, colour=group)) + 
      geom_point() + 
      coord_flip()
  })
  
  output$heatmap <- renderPlot({
    logfine('Started evaluation: heatmap')
    reaction_metab_results <- reaction_metab_results()
    contrast_1 <- input$contrast_1
    contrast_2 <- input$contrast_2
    logfine('Retrieved all reactives: heatmap')
    
    reaction_metab_results %>%
      inner_join(filter(., group==contrast_1),
                 filter(., group==contrast_2),
                 by=c('reaction', 'metabolite')) %T>% logfiner('joined',.) %>%
      select(-group.x, -group.y) %>%
      mutate(diff = flux.x-flux.y) %>%
      select(-flux.x, -flux.y)  %>%
      spread(metabolite, diff, fill=0) %T>% logfiner('spread',.)  %>%
      (function(x){
        res <- as.matrix(x %>% select(-reaction))
        rownames(res) <- x$reaction
        return(res)
      }) %T>% logfiner('matricised',.) %>%
      heatmap.plus::heatmap.plus(distfun = (function(x){dist(abs(x))}), 
                                 col=RColorBrewer::brewer.pal(11,'RdYlGn'), 
                                 main=contrast
      )
  })
  
})
