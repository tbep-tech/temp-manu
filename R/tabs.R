library(tidyverse)
library(here)
library(flextable)
library(mgcv)
library(broom)

source(here('R/funcs.R'))

# mixef mod summary table ---------------------------------------------------------------------

load(file = here('data/mixmods.RData'))

totab <- mixmods %>% 
  select(-yrstr, -yrstrse, -yrend, -yrendse) %>% 
  mutate(
    pvl = p_ast(pvl), 
    slo = as.character(round(slo, 2)), 
    salithr = gsub('^.*_', '', salithr),
    tempthr = gsub('^.*_', '', tempthr), 
    thrtyp = factor(thrtyp, levels = c('tempcnt', 'salicnt', 'bothcnt'), 
                     labels = c('Temperature', 'Salinity', 'Both')),
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) %>% 
  unite('slo', slo, pvl, sep = '') %>% 
  mutate(
    slo = ifelse(slo == 'NANA', '-', slo)
  ) %>%
  arrange(bay_segment, thrtyp) %>% 
  pivot_wider(names_from = 'thrtyp', values_from = 'slo')

mixtab <- totab %>% 
  as_grouped_data(groups = 'bay_segment') %>% 
  flextable() %>% 
  add_header_row(values = c('', 'Thresholds', 'Thresholds', rep('Slopes', 3))) %>% 
  merge_at(i = 1, j = c(2:3), part = 'header') %>% 
  merge_at(i = 1, j = c(4:6), part = 'header') %>% 
  set_header_labels(i = 2, values = c('Bay Segment', 'Temperature', 'Salinity', 'Temperature', 'Salinity', 'Both')) %>% 
  padding(padding = 0, part = 'all') %>% 
  font(part = 'all', fontname = 'Times New Roman') %>% 
  fontsize(size = 9, part = 'body')

save(mixtab, file = here('tabs/mixtab.RData'))

# mixef summary start and end number of days --------------------------------------------------

load(file = here('data/mixmods.RData'))

totab <- mixmods %>% 
  select(-slo, -pvl) %>% 
  mutate(
    yrstr = round(yrstr, 0), 
    yrstrse = paste0('(', round(yrstrse, 1), ')'),
    yrend = round(yrend, 0), 
    yrendse = paste0('(', round(yrendse, 1), ')'),
    salithr = gsub('^.*_', '', salithr),
    tempthr = gsub('^.*_', '', tempthr), 
    thrtyp = factor(thrtyp, levels = c('tempcnt', 'salicnt', 'bothcnt'), 
                    labels = c('Temperature', 'Salinity', 'Both')),
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) %>% 
  unite('yrstr', yrstr, yrstrse, sep = ' ') %>%
  unite('yrend', yrend, yrendse, sep = ' ') %>%
  pivot_longer(names_to = 'def', values_to = 'est', matches('yr')) %>% 
  unite('thrtyp', thrtyp, def) %>% 
  mutate(
    thrtyp = factor(thrtyp, 
                    levels = c('Temperature_yrstr', 'Temperature_yrend', 
                               'Salinity_yrstr', 'Salinity_yrend', 'Both_yrstr', 'Both_yrend')
    ), 
    est = ifelse(grepl('NA', est), '-', est)
  ) %>%
  arrange(bay_segment, thrtyp) %>% 
  pivot_wider(names_from = 'thrtyp', values_from = 'est')

daytab <- totab %>% 
  as_grouped_data(groups = 'bay_segment') %>% 
  flextable() %>% 
  add_header_row(values = c('', 'Thresholds', 'Thresholds', rep('Temperature', 2), 
                            rep('Salinity', 2), rep('Both', 2))) %>% 
  merge_at(i = 1, j = c(2:3), part = 'header') %>% 
  merge_at(i = 1, j = c(4:5), part = 'header') %>% 
  merge_at(i = 1, j = c(6:7), part = 'header') %>% 
  merge_at(i = 1, j = c(8:9), part = 'header') %>% 
  set_header_labels(i = 2, values = c('Bay Segment', 'Temperature', 'Salinity', rep(c('Start', 'End'), 3))) %>% 
  padding(padding = 0, part = 'all') %>% 
  font(part = 'all', fontname = 'Times New Roman') %>% 
  fontsize(size = 9, part = 'body')

save(daytab, file = here('tabs/daytab.RData'))

# GAM summary ---------------------------------------------------------------------------------

# load(file = here('data/cmbmod.RData'))
# 
# totab <- cmbmod %>% 
#   tidy() %>% 
#   mutate(
#     p.value = p_ast(p.value), 
#     term = gsub('^ti\\(|\\)$', '', term),
#     Effect = case_when(
#       grepl('\\,', term) ~ 'Interaction',
#       T ~ 'Main'
#     )
#   ) %>% 
#   mutate_if(is.numeric, round, 2) %>% 
#   unite('statistic', statistic, p.value, sep = '') %>% 
#   select(
#     Effect, 
#     Term = term, 
#     edf, 
#     ref.df, 
#     `F` = statistic
#   )
# 
# gamtab <- totab %>% 
#   as_grouped_data(groups = 'Effect') %>% 
#   flextable() %>% 
#   padding(padding = 0, part = 'all') %>% 
#   font(part = 'all', fontname = 'Times New Roman') %>% 
#   align(j = 3:5, align = 'left', part = 'all')
# 
# save(gamtab, file = here('tabs/gamtab.RData'))

