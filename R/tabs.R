library(tidyverse)
library(here)
library(flextable)
library(mgcv)
library(broom)

source(here('R/funcs.R'))

# mixef mod summary table and start/end -------------------------------------------------------

load(file = here('data/mixmods.RData'))

salthr <- '25'
tmpthr <- '30'

totab <- mixmods %>% 
  filter(salithr == paste0('sali_', salthr)) %>% 
  filter(tempthr == paste0('temp_', tmpthr)) %>%  
  select(-tempthr, -salithr) %>% 
  mutate(
    pvl = p_ast(pvl), 
    slo = as.character(round(slo, 2)), 
    yrstr = round(yrstr, 0), 
    yrstrse = paste0('(', round(yrstrse, 1), ')'),
    yrend = round(yrend, 0), 
    yrendse = paste0('(', round(yrendse, 1), ')'),
    thrtyp = factor(thrtyp, levels = c('tempcnt', 'salicnt', 'bothcnt'), 
                    labels = c(paste('Temperature >', tmpthr), paste('Salinity <', salthr), 'Both')),
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) %>% 
  unite('slo', slo, pvl, sep = '') %>% 
  unite('yrstr', yrstr, yrstrse, sep = ' ') %>%
  unite('yrend', yrend, yrendse, sep = ' ') %>% 
  mutate(
    slo = ifelse(slo == 'NANA', '-', slo),
    yrstr = ifelse(grepl('NA', yrstr), '-', yrstr),
    yrend = ifelse(grepl('NA', yrend), '-', yrend)
  ) %>%
  arrange(bay_segment, thrtyp)

mixdaytab <- totab %>% 
  as_grouped_data(groups = 'bay_segment') %>% 
  flextable() %>% 
  set_header_labels(i = 1, values = c('Bay Segment', 'Threshold', 'Slope', 'Start', 'End')) %>% 
  padding(padding = 0, part = 'all') %>% 
  width(j = 2, width = 1.5) %>% 
  font(part = 'all', fontname = 'Times New Roman')

save(mixdaytab, file = here('tabs/mixdaytab.RData'))

# supp gam performance, temp/sal only ---------------------------------------------------------

load(file = here('data/gamfit.RData'))

fittab <- gamfit %>% 
  mutate(
    bay_segment = as.character(bay_segment),
    bay_segment = ifelse(duplicated(bay_segment), '', bay_segment), 
    .by = 'parameter'
  ) %>% 
  group_nest(parameter) %>% 
  mutate(
    data = purrr::map(data, function(x){
      
      x %>% 
        knitr::kable(
          col.names = c('Bay Segment', 'Station', 'Lon', 'Lat', 'AIC', 'GCV', 'R$^2$')
        ) 
      
    })
  )

supptempfittab <- fittab %>% 
  filter(parameter == 'temp') %>% 
  pull('data') %>% 
  .[[1]]
suppsalifittab <- fittab %>% 
  filter(parameter == 'sal') %>% 
  pull('data') %>% 
  .[[1]]

save(supptempfittab, file = here('tabs/supptempfittab.RData'))
save(suppsalifittab, file = here('tabs/suppsalifittab.RData'))

# supp1 mixef mod summary table ----------------------------------------------------------------

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

supp1mixtab <- totab %>% 
  as_grouped_data(groups = 'bay_segment') %>% 
  flextable() %>% 
  add_header_row(values = c('', 'Thresholds', 'Thresholds', rep('Slopes', 3))) %>% 
  merge_at(i = 1, j = c(2:3), part = 'header') %>% 
  merge_at(i = 1, j = c(4:6), part = 'header') %>% 
  set_header_labels(i = 2, values = c('Bay Segment', 'Temperature', 'Salinity', 'Temperature', 'Salinity', 'Both')) %>% 
  padding(padding = 0, part = 'all') %>% 
  font(part = 'all', fontname = 'Times New Roman') %>% 
  fontsize(size = 9, part = 'body')

save(supp1mixtab, file = here('tabs/supp1mixtab.RData'))

# supp2 mixef mod summary table ----------------------------------------------------------------

load(file = here('data/suppmixmods.RData'))

totab <- suppmixmods %>% 
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

supp2mixtab <- totab %>% 
  as_grouped_data(groups = 'bay_segment') %>% 
  flextable() %>% 
  add_header_row(values = c('', 'Thresholds', 'Thresholds', rep('Slopes', 3))) %>% 
  merge_at(i = 1, j = c(2:3), part = 'header') %>% 
  merge_at(i = 1, j = c(4:6), part = 'header') %>% 
  set_header_labels(i = 2, values = c('Bay Segment', 'Temperature', 'Salinity', 'Temperature', 'Salinity', 'Both')) %>% 
  padding(padding = 0, part = 'all') %>% 
  font(part = 'all', fontname = 'Times New Roman') %>% 
  fontsize(size = 9, part = 'body')

save(supp2mixtab, file = here('tabs/supp2mixtab.RData'))

# supp1 mixef summary start and end number of days --------------------------------------------

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

supp1daytab <- totab %>% 
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

save(supp1daytab, file = here('tabs/supp1daytab.RData'))

# supp2 mixef summary start and end number of days ---------------------------------------------

load(file = here('data/suppmixmods.RData'))

totab <- suppmixmods %>% 
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

supp2daytab <- totab %>% 
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

save(supp2daytab, file = here('tabs/supp2daytab.RData'))

