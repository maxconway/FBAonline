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
library(ROI)
library(ROI.plugin.ecos)
safely(library)('ROI.plugin.glpk')
library(visNetwork)
library(stringr)
library(igraph)
library(ggdendro)

options(shiny.sanitize.errors = FALSE)

basicConfig('FINEST')
addHandler(writeToConsole)

dummytable <- tibble(lowbnd=double(), uppbnd=double(), abbreviation=character(), name=character(), equation=character(), obj_coef=double())
dummyexpanded <- dummytable %>% fbar::reactiontbl_to_expanded()
dummytableflux <- dummytable %>% mutate(flux=double())

shinyServer(function(input, output, session) {
  #### Inputs
  
  model1 <- reactive({
    logfine('Started evaluation: model1')
    input$reevaluate
    model_url_1 <- input$model_url_1
    logfine('Finished loading: model1')
    possibly(gsheet::gsheet2tbl, otherwise = dummytable)(model_url_1)
  })
  
  model2 <- reactive({
    logfine('Started evaluation: model2')
    input$reevaluate
    model_url_2 <- input$model_url_2
    logfine('Finished loading: model2')
    possibly(gsheet::gsheet2tbl, otherwise = dummytable)(model_url_2)
  })
  
  ### Processing
  
  model1_parsed <- reactive({
    logfine('Started evaluation: model1_parsed')
    model <- model1()
    logfine('Finished loading: model1_parsed')
    possibly(fbar::reactiontbl_to_expanded, dummyexpanded)(model)
  })
  
  model2_parsed <- reactive({
    logfine('Started evaluation: model2_parsed')
    model <- model2()
    logfine('Finished loading: model2_parsed')
    possibly(fbar::reactiontbl_to_expanded, dummyexpanded)(model)
  })
  
  model1_reactions_with_fluxes <- reactive(label = 'model1_reactions_with_fluxes', x={
    logfine('Started evaluation: model1_reactions_with_fluxes')
    model <- model1()
    logfine('Finished loading: model1_reactions_with_fluxes')
    
    model %>%
      possibly(find_fluxes_df, dummytableflux)(do_minimization = FALSE) %>%
      bind_rows(dummytableflux,.)
  })
  
  model2_reactions_with_fluxes <- reactive(label = 'model2_reactions_with_fluxes', x={
    logfine('Started evaluation: model2_reactions_with_fluxes')
    model <- model2()
    logfine('Finished loading: model2_reactions_with_fluxes')
    model %>%
      possibly(find_fluxes_df, dummytableflux)(do_minimization = FALSE) %>%
      bind_rows(dummytableflux,.)
  })
  
  model1_filtered_reactions_with_fluxes <- reactive(label = 'model1_filtered_reactions_with_fluxes', x={
    logfine('Started evaluation: model1_filtered_reactions_with_fluxes')
    reactions_with_fluxes <- model1_reactions_with_fluxes()
    reaction_filter <- input$reaction_filter
    logfine('Finished loading: model1_filtered_reactions_with_fluxes')
    
    reactions_with_fluxes %>%
      filter((abbreviation %in% reaction_filter) | is_empty(reaction_filter))
  })
  
  model2_filtered_reactions_with_fluxes <- reactive(label = 'model2_filtered_reactions_with_fluxes', x={
    logfine('Started evaluation: model2_filtered_reactions_with_fluxes')
    reactions_with_fluxes <- model2_reactions_with_fluxes()
    reaction_filter <- input$reaction_filter
    logfine('Finished loading: model2_filtered_reactions_with_fluxes')
    
    reactions_with_fluxes %>%
      filter((abbreviation %in% reaction_filter) | is_empty(reaction_filter))
  })
  
  clusters <- reactive(label = 'clusters', x = {
    logfine('Started evaluation: clusters')
    model1_parsed <- model1_parsed()
    logfine('Finished loading: clusters')
    
    if(!all(map_int(model1_parsed, nrow)>0)){
      return(NULL)
    }
    
    model1_parsed %>%
      `$`(stoich) %>%
      select(-stoich) %>%
      inner_join(.,.,by='met') %>%
      select(from=abbreviation.x, to=abbreviation.y) %>%
      unique() %>%
      filter(from!=to) %>%
      graph_from_data_frame(vertices = model1_parsed$rxns) %>%
      as.undirected(mode='collapse') %>%
      cluster_walktrap() %>%
      as.dendrogram()
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
  
  observe(priority=-100, x={
    clusters <- clusters()
    if(is.null(clusters)){
      return(NULL)
    }
    updateSliderInput(session, 'cut_h', max = attr(clusters,'height'))
    updateSliderInput(session, 'min_members', max = attr(clusters,'members'))
  })
  
  #### Outputs
  
  output$model_status_1 <- renderText({
    logfine('Started evaluation: model_status_1')
    model1 <- model1()
    model1_parsed <- model1_parsed()
    model1_reactions_with_fluxes <- model1_reactions_with_fluxes()
    logfine('Finished loading: model_status_1')
    if(nrow(model1)==0){
      return('Waiting for model')
    }
    if(nrow(model1)>0 & is_null(model1_parsed)){
      return('Model failed to parse')
    }
    if(nrow(model1)>0 & !is_null(model1_parsed) & !has_name(model1_reactions_with_fluxes, 'flux')){
      return('Model failed to evaluate')
    }
    if(nrow(model1)>0 & !is_null(model1_parsed) & has_name(model1_reactions_with_fluxes, 'flux')){
      return('Model OK')
    }
    
  })
  
  output$model_table <- renderDataTable({
    logfine('Started evaluation: model_table')
    filtered_reactions_with_fluxes <- model1_filtered_reactions_with_fluxes()
    logfine('Finished loading: model_table')
    
    filtered_reactions_with_fluxes %>%
      mutate(flux = signif(flux, 3)) %>%
      #mutate_if(is_character, str_trunc, width=40) %>% # Produces weird error
      arrange(desc(abs(obj_coef)), abbreviation)
  })
  
  output$bar_chart <- renderPlot({
    logfine('Started evaluation: bar_chart')
    model1_filtered_reactions_with_fluxes <- model1_filtered_reactions_with_fluxes()
    model2_filtered_reactions_with_fluxes <- model2_filtered_reactions_with_fluxes()
    logfine('Finished loading: bar_chart')
    
    if((nrow(model1_filtered_reactions_with_fluxes) + nrow(model1_filtered_reactions_with_fluxes)) == 0){
      return('Waiting for models...')
    }
    
    bind_rows(model1_filtered_reactions_with_fluxes %>% select(abbreviation, flux) %>% mutate(name='model1'),
              model2_filtered_reactions_with_fluxes %>% select(abbreviation, flux) %>% mutate(name='model2'),
    ) %>%
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
        decompose_metabolites() %>%
        mutate(present = 'present') %>%
        tidyr::spread(compartment, present, fill='absent')
    }else{
      result <- NULL
    }
    return(result)
  })
  
  output$cluster_chart <- renderPlot({
    logfine('Started evaluation: cluster_chart')
    model1_parsed <- model1_parsed()
    clusters <- clusters()
    logfine('Finished loading: cluster_chart')
    validate(need(all(map_int(model1_parsed, nrow)>0), 'model not yet parsed'),
             need(clusters, 'clusters not yet calculated')
    )
    
    b <- clusters %>% 
      cut(h = input$cut_h) %>%
      `$`(lower) %>%
      discard(is.leaf) %>%
      keep(~attr(.x,'members') >= input$min_members) %>%
      map(dendro_data) %>%
      transpose %>%
      `[`(c('labels','segments')) %>%
      map(function(x){quietly(bind_rows)(x, .id='tree')$result}) 
    
    ggplot(b$segments, aes(x=x,y=y)) + 
      geom_segment(aes(xend=xend, yend=yend), data=b$segments, alpha=0.5) + 
      geom_text(aes(label=label), data=b$labels, hjust='left', nudge_x = 0.15) + 
      facet_wrap(~tree, scales='free_y') + 
      ylab('Height') +
      coord_flip() +
      theme_bw() +
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y = element_blank(),
            strip.text = element_blank(),
            strip.background = element_blank())
  })
  
  output$network <- renderVisNetwork({
    logfine('Started evaluation: network')
    model1_parsed <- model1_parsed()
    reaction_filter <- input$reaction_filter
    logfine('Finished loading: network')
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