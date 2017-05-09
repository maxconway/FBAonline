# Metabolite helper functions

library(stringr)

expand_metabolites <- function(met_table){
  met_table %>%
    transmute(chemical = str_replace(met, '\\[[a-zA-Z0-9]+]$', ''),
              compartment = str_extract(met, '\\[[a-zA-Z0-9]+]$') %>% 
                coalesce('') %>%
                str_replace_all('\\[|]', ''))
}

condense_metabolites <- function(expanded_metabolites){
  expanded_metabolites %>%
    transmute(met = str_c(chemical, '[', compartment, ']'))
}