library(tidyverse)
library(here)
library(mgcv)
library(broom)
library(tbeptools)
library(flextable)
library(knitr)

source(here('R/funcs.R'))

# use mid salinity for 1975, bottom missing
data(epcdata)
epcdat <- epcdata %>% 
  mutate(
    Sal_Bottom_ppth = case_when(
      yr == 1975 & is.na(Sal_Bottom_ppth) ~ Sal_Mid_ppth, 
      T ~ Sal_Bottom_ppth
    )
  ) %>% 
  filter(yr < 2023)

# dataset summaries ---------------------------------------------------------------------------

totab <- tibble(
  Dataset = c('SWFWMD aerial maps', 'Transect data', 'EPC', 'FIM', 'PDEM', 'Tampa International Airport', 'SWFWMD precipitation'),
  Description = c('Seagrass coverage in acres', 'Seagrass frequency occurrence by species', 'Water quality monitoring samples', 'Nearshore temperature and salinity, seagrass species and cover', 'Water quality and seagrass presence/absence', 'Air temperature', 'Area-weighted precipitation for the wet season (June-September) for the Tampa Bay watershed'),
  Temporal = c('1988 - 2022, biennial', '1999 - 2022, annual', '1975 - 2022, monthly', '1996 to 2022, monthly', '2003 - 2022, monthly', '1975 - 2022, annual', '1975-2022, annual'),
  Spatial = c('Whole bay', 'Whole bay, 62 transects', 'Whole bay, fixed sites', 'Whole bay shallow, stratified random sites', 'Old Tampa Bay, stratified random sites', '27.979$^\\circ$N, 82.535$^\\circ$W', 'Whole watershed'),
  Analysis = c('Biennial trends by bay segment, visual only', 'Annual trends by bay segment and species, comparison with temperature, salinity, and light attenuation as stressor metrics or observed data at annual scale', 'Trends in annual change and seasonal Kendall tests, estimate of stressor metrics as number of days above/below threshold', 'Trends in annual observed temperature, salinity, comparison to annual seagrass % cover', 'Trends in annual observed temperature, salinity, comparison to annual seagrass frequency occurrence', 'Annual trend', 'Annual trend')
)

dattab <- knitr::kable(totab)

save(dattab, file = here('tabs/dattab.RData'))

# linear trend summaries ----------------------------------------------------------------------

load(file = here('data/lintrnds.RData'))

totab <- lintrnds %>% 
  mutate_at(vars(strvest, endvest), round, 1) %>% 
  mutate(
    chng = endvest - strvest,
    slo = round(slo, 2), 
    slose = paste0('(', round(slose, 2), ')'),
    pval = p_ast(pval), 
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')), 
    org = factor(org, levels = c('EPC', 'FIM', 'PDEM'))
  ) %>% 
  unite('slo', slo, slose, sep = ' ') %>% 
  select(
    var, 
    `Start year` = yrstr,
    `Bay segment` = bay_segment, 
    Dataset = org, 
    `Change / year` = slo,
    `Start value` = strvest, 
    `End value` = endvest, 
    `Total change` = chng
  ) %>% 
  arrange(var, `Start year`, `Bay segment`, Dataset) %>% 
  filter(`Change / year` != 'NA (NA)') %>% 
  mutate(
    `Bay segment` = ifelse(duplicated(`Bay segment`), '', as.character(`Bay segment`)), 
    .by = c(var, `Start year`)
  ) %>% 
  mutate(
    `Start year` = ifelse(duplicated(`Start year`), '', as.character(`Start year`)),
    .by = var
  )

# save temp and sali trend tables  
tab <- totab %>% 
  group_nest(var) %>% 
  deframe() %>% 
  iwalk(function(data, ind){
    
    outnm <- paste0(ind, 'trndtab')
    
    out <- knitr::kable(data)
      
    assign(outnm, out)
    fl <- here(paste0('tabs/', outnm, '.RData'))
    save(list = outnm, file = fl, compress = 'xz')
    
  })

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
  select(-pvl) %>% 
  unite('yrstr', yrstr, yrstrse, sep = ' ') %>%
  unite('yrend', yrend, yrendse, sep = ' ') %>% 
  mutate(
    slo = ifelse(is.na(slo), '-', slo),
    yrstr = ifelse(grepl('NA', yrstr), '-', yrstr),
    yrend = ifelse(grepl('NA', yrend), '-', yrend)
  ) %>%
  arrange(bay_segment, thrtyp)

mixdaytab <- totab %>% 
  mutate(
    bay_segment = ifelse(duplicated(bay_segment), '', as.character(bay_segment))
  ) %>%
  knitr::kable(col.names = c('Bay Segment', 'Threshold', 'Slope', 'Start', 'End'))

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
  select(-yrstr, -yrstrse, -yrend, -yrendse, -pvl) %>% 
  mutate(
    slo = as.character(round(slo, 2)), 
    salithr = gsub('^.*_', '', salithr),
    tempthr = gsub('^.*_', '', tempthr), 
    thrtyp = factor(thrtyp, levels = c('tempcnt', 'salicnt', 'bothcnt'), 
                     labels = c('Temperature', 'Salinity', 'Both')),
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) %>% 
  mutate(
    slo = ifelse(is.na(slo), '-', slo)
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
  select(-yrstr, -yrstrse, -yrend, -yrendse, -pvl) %>% 
  mutate(
    slo = as.character(round(slo, 2)), 
    salithr = gsub('^.*_', '', salithr),
    tempthr = gsub('^.*_', '', tempthr), 
    thrtyp = factor(thrtyp, levels = c('tempcnt', 'salicnt', 'bothcnt'), 
                    labels = c('Temperature', 'Salinity', 'Both')),
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) %>% 
  mutate(
    slo = ifelse(is.na(slo), '-', slo)
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

# gam summaries -------------------------------------------------------------------------------

load(file = here::here('data/sgmods.RData'))
load(file = here::here('data/sgmodsum.RData'))

suppepccap1 <- paste0('Summary of smoother terms in the Generalized Additive Model used to evaluate seagrass response for the EPC data in relation to temperature (temp) and salinity (sal) stress metrics, with additional smoothers for year (yr) and light attenuation (la).  Separate smoothers were fit for each bay segment.  $s$ = invidual smoother, $ti$ = interaction term. OTB: Old Tampa Bay, HB: Hillsborough Bay, MTB: Middle Tampa Bay. ', modtxt_fun(sgmodsum)$epcmod1, '.')  
suppepcmod1tab <- gam_table(sgmods$epcmod1, cap = suppepccap1)
save(suppepcmod1tab, file = here('tabs/suppepcmod1tab.RData'))

suppepccap2 <- paste0('Summary of smoother terms in the Generalized Additive Model used to evaluate seagrass response for the EPC data in relation to "both" stressors (both temperature above and salinity below thresholds), with additional smoothers for year (yr) and light attenuation (la).  Separate smoothers were fit for each bay segment.  $s$ = invidual smoother, $ti$ = interaction term. OTB: Old Tampa Bay, HB: Hillsborough Bay, MTB: Middle Tampa Bay. ', modtxt_fun(sgmodsum)$epcmod2, '.')    
suppepcmod2tab <- gam_table(sgmods$epcmod2, cap = suppepccap2)
save(suppepcmod2tab, file = here('tabs/suppepcmod2tab.RData'))

suppfimcap <- paste0('Summary of smoother terms in the Generalized Additive Model used to evaluate seagrass response for the FIM data in relation to temperature (temp), salinity (sal), and year (yr).  Separate smoothers were fit for each bay segment.  $s$ = invidual smoother, $ti$ = interaction term. OTB: Old Tampa Bay, HB: Hillsborough Bay, MTB: Middle Tampa Bay. ' , modtxt_fun(sgmodsum)$fimmod, '.')   
suppfimmodtab <- gam_table(sgmods$fimmod, cap = suppfimcap)
save(suppfimmodtab, file = here('tabs/suppfimmodtab.RData'))

supppincocap <- paste0('Summary of smoother terms in the Generalized Additive Model used to evaluate seagrass response for the PDEM data in relation to temperature (temp), salinity (sal), and year (yr).  $s$ = invidual smoother, $ti$ = interaction term. ', modtxt_fun(sgmodsum)$pincomod, '.')    
supppincomodtab <- gam_table(sgmods$pincomod, cap = supppincocap)
save(supppincomodtab, file = here('tabs/supppincomodtab.RData'))
