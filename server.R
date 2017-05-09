library(shiny)
library(dplyr)
library(purrr)
library(ggplot2)
library(tibble)
library(gsheet)
library(fbar)
library(ROI)
library(magrittr)
library(logging)
library(ROI.plugin.ecos)
#safely(library)('ROI.plugin.glpk')
library(visNetwork)
library(stringr)

source('R/methelpers.R')

basicConfig('FINEST')
addHandler(writeToConsole)

shinyServer(function(input, output, session) {
  #### Inputs
  
  model1 <- reactive({
    logfine('Started evaluation: model1')
    input$reevaluate
    model_url_1 <- input$model_url_1
    logfine('Finished loading: model1')
    if(nchar(model_url_1)<10){
      loginfo('url too short')
      result <- tibble::tribble(~abbreviation, ~equation, ~lowbnd, ~uppbnd, ~obj_coef)
    }else{
      result <- tryCatch(gsheet::gsheet2tbl(model_url_1),
                         error = function(e){
                           validate(need(FALSE, label = 'URL for model 1', message='URL for model 1 is not a google spreadsheet'))
                           return(tibble::tribble(~abbreviation, ~equation, ~lowbnd, ~uppbnd, ~obj_coef))
                         })
    }
    result
  })
  
  model2 <- reactive({
    logfine('Started evaluation: model2')
    input$reevaluate
    model_url_2 <- input$model_url_2
    logfine('Finished loading: model2')
    if(nchar(model_url_2)<10){
      loginfo('url too short')
      result <- tibble::tribble(~abbreviation, ~equation, ~lowbnd, ~uppbnd, ~obj_coef)
    }else{
      result <- tryCatch(gsheet::gsheet2tbl(model_url_2),
                         error = function(e){
                           validate(need(FALSE, label = 'URL for model 2', message='URL for model 2 is not a google spreadsheet'))
                           return(tibble::tribble(~abbreviation, ~equation, ~lowbnd, ~uppbnd, ~obj_coef))
                         })
    }
    result
  })
  
  ### Processing
  
  model1_parsed <- reactive({
    logfine('Started evaluation: model1_parsed')
    model1 <- model1()
    logfine('Finished loading: model1_parsed')
    if(!is.null(model1) & nrow(model1) > 0){
      result <- fbar::reactiontbl_to_expanded(model1)
    }else{
      result <- NULL
    }
    result
  })
  
  model2_parsed <- reactive({
    logfine('Started evaluation: model2_parsed')
    model2 <- model2()
    logfine('Finished loading: model2_parsed')
    if(!is.null(model2) & nrow(model2) > 0){
      result <- fbar::reactiontbl_to_expanded(model2)
    }else{
      result <- NULL
    }
    result
  })
  
  model1_evaluated <- reactive({
    logfine('Started evaluation: model1_evaluated')
    model_parsed <- model1_parsed()
    logfine('Finished loading: model1_evaluated')
    if(!is.null(model_parsed) & is.list(model_parsed)){
      result <- ROI_solve(expanded_to_ROI(model_parsed))
      logfine(paste('mod1 LP status: ',result$status$code))
    }else{
      result <- NULL
    }
    result
  })
  
  model2_evaluated <- reactive({
    logfine('Started evaluation: model2_evaluated')
    model_parsed <- model2_parsed()
    logfine('Finished loading: model2_evaluated')
    if(!is.null(model_parsed) & is.list(model_parsed)){
      result <- ROI_solve(expanded_to_ROI(model_parsed))
  logfine(paste('mod2 LP status: ',result$status$code))
    }else{
      result <- NULL
    }
    result
  })
  
  filtered_reactions_with_fluxes <- reactive({
    logfine('Started evaluation: filtered_reactions_with_fluxes')
    model1_parsed <- model1_parsed()
    model2_parsed <- model2_parsed()
    model1_evaluated <- model1_evaluated()
    model2_evaluated <- model2_evaluated()
    reaction_filter <- input$reaction_filter
    logfine('Finished loading: filtered_reactions_with_fluxes')
    
    list(
      model1 = list(parsed = model1_parsed, evaluated = model1_evaluated, name='Model 1'),
      model2 = list(parsed = model2_parsed, evaluated = model2_evaluated, name='Model 2')
    ) %>%
      keep(~is_list(.[['parsed']])) %>%
      map_df(function(x){
        x$parsed %>% 
          fbar::expanded_to_reactiontbl() %>% 
          inner_join(data_frame(abbreviation = x$parsed$rxns$abbreviation, flux = ROI::solution(x$evaluated)))
        }, .id='name') %>%
      bind_rows(tribble(~abbreviation, ~equation, ~lowbnd, ~uppbnd, ~obj_coef, ~flux),.) %>% # normalizes empty tables
      filter((abbreviation %in% reaction_filter) | is_empty(reaction_filter)) %>%
      select(flux, everything())
  })

  ### Observations

  observe({
    logfine('Started evaluation: observe')
    model1_parsed <- model1_parsed()
    model2_parsed <- model2_parsed()
    logfine('Finished loading: observe')
    
    abbreviations <- list(model1_parsed, model2_parsed) %>%
      keep(is_list) %>%
      map_df(extract2, 'rxns') %>%
      extract2('abbreviation') %>%
      unique()
      
    updateSelectizeInput(session, 'reaction_filter', choices = abbreviations)
  })
  
  #### Outputs
  
  output$model_table <- renderDataTable({
    logfine('Started evaluation: model_table')
    filtered_reactions_with_fluxes <- filtered_reactions_with_fluxes()
    logfine('Finished loading: model_table')

    filtered_reactions_with_fluxes %>%
      mutate(flux = signif(flux, 3)) %>%
      arrange(desc(abs(obj_coef)), abbreviation)
  })
  
  output$bar_chart <- renderPlot({
    logfine('Started evaluation: bar_chart')
    filtered_reactions_with_fluxes <- filtered_reactions_with_fluxes()
    logfine('Finished loading: bar_chart')
    filtered_reactions_with_fluxes %>%
      ggplot(aes(x=abbreviation, y=flux, colour=name)) + 
      geom_point() + 
      coord_flip()
  })
  
  output$metabolite_table <- renderDataTable({
    logfine('Started evaluation: metabolite_table')
    model1_parsed <- model1_parsed()
    if(!is.null(model1_parsed)){
      result <- model1_parsed %>%
        getElement('mets') %>%
        expand_metabolites() %>%
        mutate(present = 'present') %>%
        tidyr::spread(compartment, present, fill='absent')
    }else{
      result <- NULL
    }
    return(result)
  })
  
  output$network <- renderVisNetwork({
    logfine('Started evaluation: network')
    model1_parsed <- model1_parsed()
    reaction_filter <- input$reaction_filter
    logfine('Finished loading: network')
    
    if('model missing' %in% model1_parsed){
      model1_parsed <- list(rxns = tibble::tribble(~abbreviation, ~equation, ~lowbnd, ~uppbnd, ~obj_coef),
                            mets = tibble::tribble(~met),
                            stoich = tibble::tribble(~stoich, ~abbreviation, ~met))
    }
      with(model1_parsed,{
        rxns <- rxns %>% filter(abbreviation %in% reaction_filter)
        stoich <- stoich %>% filter(abbreviation %in% reaction_filter)
        mets <- mets %>% filter(met %in% stoich$met)
        
        if(nrow(rxns)>50){
          shiny::validate('Too many reactions. Filter down to 50 or less.')
        }
        
        visNetwork(nodes=bind_rows(rxns %>% 
                                     transmute(id=as.character(abbreviation),
                                               type='rxn'),
                                   mets %>% 
                                     transmute(id=met, type='met')) %>%
                     mutate(title=id,
                            group=type),
                   edges=stoich %>%
                     mutate(from = ifelse(stoich>0, as.character(abbreviation), met),
                            to = ifelse(stoich>0, met, as.character(abbreviation)),
                            value = abs(stoich),
                            title = abs(stoich),
                            arrows='middle')) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visEdges(smooth = FALSE) %>%
        visPhysics(stabilization = FALSE)
      })
    
  })
  
})