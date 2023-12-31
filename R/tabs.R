library(tidyverse)
library(here)
library(mgcv)
library(broom)
library(tbeptools)

source(here('R/funcs.R'))

# linear trend summaries ----------------------------------------------------------------------

data(epcdata)
data(fimsgtempdat)
data(pincotemp)

epctmp <- epcdata %>% 
  select(bay_segment, epchc_station, SampleTime, yr, matches('Top|Bottom')) %>% 
  filter(yr < 2023) %>% 
  pivot_longer(names_to = 'var', values_to = 'val', matches('Top|Bottom')) %>% 
  mutate(
    var = factor(var, 
                 levels = c(c("Sal_Top_ppth", "Sal_Bottom_ppth", "Temp_Water_Top_degC", "Temp_Water_Bottom_degC"
                 )), 
                 labels = c("Sal_Top", "Sal_Bottom", "Temp_Top", "Temp_Bottom")
    ),
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) %>% 
  separate(var, c('var', 'loc')) %>% 
  mutate(
    var = factor(var, levels = c('Temp', 'Sal'), labels = c('temp', 'sal'))
  ) %>% 
  filter(loc %in% 'Bottom') %>% 
  filter(!is.na(val)) %>% 
  summarise(
    avev = mean(val, na.rm = T),
    .by = c('bay_segment', 'yr', 'var') 
  ) %>% 
  mutate(
    org = 'EPC'
  )

fimtmp <- fimsgtempdat %>% 
  select(date, temp, sal, bay_segment) %>% 
  mutate(
    yr = year(date), 
    mo = month(date)
  ) %>% 
  pivot_longer(temp:sal, names_to = 'var', values_to = 'val') %>% 
  summarise(
    avev = mean(val, na.rm = T),
    .by = c(bay_segment, yr, var)
  ) %>% 
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')),
    var = factor(var, levels = c('temp', 'sal')), 
    org = 'FIM'
  )

pincotmp <- pincotemp %>% 
  pivot_longer(temp:sal, names_to = 'var', values_to = 'val') %>% 
  summarise(
    avev = mean(val, na.rm = T), 
    .by = c(yr, var)
  ) %>% 
  mutate(
    var = factor(var, levels = c('temp', 'sal')), 
    bay_segment = 'OTB', 
    org = 'PDEM'
  )

trnds <- bind_rows(epctmp, fimtmp, pincotmp) %>% 
  group_nest(org, bay_segment, var) %>% 
  crossing(yrstr = c(1976, 1996, 2004)) %>% 
  mutate(
    i = 1:n(),
    data = purrr::pmap(list(i, yrstr, data), function(i, yrstr, data){
      
      cat(i, '\n')
      
      out <- tibble(
        slo = NA, 
        pval = NA, 
        strvest = NA,
        endvest = NA
      )
      
      if(!yrstr %in% unique(data$yr))
        return(out)
      
      tomod <- data %>% 
        filter(yr >= yrstr) %>% 
        arrange(yr)
      
      # model fit and results
      mod <- lm(avev ~ yr, data = tomod)
      coef <- summary(mod)$coefficients
      prds <- data.frame(predict(mod, se = T))
      prds$ci <- 1.96 * prds$se.fit
      strv <- prds[1, ]
      endv <- prds[nrow(prds), ]
      
      # output
      out$slo <- coef[2, 1]
      out$pval <- coef[2, 4]
      out$strvest <- strv$fit
      out$endvest <- endv$fit
      
      return(out)
      
    })
  ) %>% 
  unnest('data')

totab <- trnds %>% 
  mutate_at(vars(strvest, endvest), round, 1) %>% 
  mutate(
    chng = endvest - strvest,
    slo = round(slo, 2), 
    pval = p_ast(pval), 
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')), 
    org = factor(org, levels = c('EPC', 'FIM', 'PDEM'))
  ) %>% 
  unite('slo', slo, pval, sep = '') %>% 
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
  filter(`Change / year` != 'NANA') %>% 
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


# glm performance -----------------------------------------------------------------------------

load(file = here("data/sgmods.RData"))

# epcmod1
modtab <- sgmods %>% 
  lapply(tidy) %>% 
  enframe() %>% 
  unnest(value) %>%
  filter(term != '(Intercept)') %>% 
  mutate(p.value = p_ast2(p.value)) %>% 
  mutate_if(is.numeric, formatC, digits = 2) %>% 
  mutate(
    std.error = paste0('(', std.error, ')'), 
    name = case_when(
      name == 'epcmod1' ~ 'EPC 1', 
      name == 'epcmod2' ~ 'EPC 2',
      name == 'fimmod' ~ 'FIM',
      name == 'pincomod' ~ 'PDEM'
    ), 
    term = gsub('temp', 'Temp', term), 
    term = gsub('sal', 'Sal', term),
    term = gsub('both', 'Both', term),
    term = gsub('yrcatDecline \\(post - 2016\\)', 'Time_post', term),
    term = gsub('yrcatDecline \\(post - 2016\\)', 'Time_post', term), 
    term = gsub('bay_segment', 'Baysegment_', term),
    name = ifelse(duplicated(name), '', name)
  ) %>% 
  unite('estimate', estimate, std.error, sep = ' ') %>% 
  unite('statistic', statistic, p.value, sep = '') %>%
  select(Model = name, Term = term, Estimate = estimate, `t-statistic` = statistic) %>% 
  knitr::kable()

save(modtab, file = here('tabs/modtab.RData'))

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

