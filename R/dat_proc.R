library(wqtrends)
library(tbeptools)
library(tidyverse)
library(lubridate)

tomod <- epcdata %>% 
  select(bay_segment, station = epchc_station, date = SampleTime, value = Temp_Water_Top_degC) %>% 
  filter(year(date) < 2023) %>% 
  filter(year(date) > 1975) %>% 
  mutate(
    date = as.Date(date), 
    cont_year = decimal_date(date), 
    yr = year(date),
    param = 'Surface temperature (C)'
  ) %>% 
  filter(bay_segment == 'OTB') %>% 
  group_nest(bay_segment, station) %>% 
  mutate(
    mod = purrr::map(data, anlz_gam, trans = 'ident'), 
    fit = purrr::map(mod, anlz_fit), 
    prd = purrr::map(mod, anlz_prdday)
  )

toplo <- tomod %>% 
  select(-mod) %>% 
  mutate(
    cnts = purrr::map(prd, function(x){
      
      tibble(thr = c(10:35)) %>%
        group_nest(thr) %>%
        mutate(data = list(x)) %>%
        mutate(
          data = purrr::pmap(list(data, thr), function(data, thr){
            
            data %>%
              filter(yr < 2023) %>%
              summarise(
                cnt = sum(value > thr),
                .by = 'yr'
              )
            
          })
        ) %>% 
        unnest('data')
      
    })
  )


ggplot(toplo, aes(x = yr, y = cnt, group = thr)) +
  geom_point() + 
  geom_smooth(method = 'lm', se = F) +
  facet_grid(~thr) + 
  labs(
    y = 'Days per year with temp exceeding threshold'
  )


