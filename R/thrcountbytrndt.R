

##
# approximate sample month of transects by year, segment

transectocc <- anlz_transectocc(transect) 
transectavespp <- transectocc %>% 
  anlz_transectavespp(by_seg = TRUE) %>% 
  filter(bay_segment %in% c('HB', 'OTB', 'MTB')) %>% 
  filter(Savspecies %in% 'total') %>% 
  select(yr, bay_segment, foest)

# mean month of tranect sampling in HB
trnptssub <- trnpts %>%
  st_set_geometry(NULL) %>%
  select(
    Transect = TRAN_ID,
    bay_segment
  ) %>%
  filter(bay_segment %in% c('OTB', 'HB', 'MTB'))
trndts <- transectocc %>%
  ungroup() %>% 
  filter(Savspecies == 'total') %>%
  inner_join(trnptssub, by = 'Transect') %>%
  mutate(
    month = month(Date), 
    year = year(Date)
  ) %>% 
  summarise(
    date = ymd(mean(Date)), # use ymd because averaging dates leaves decimals 
    .by = c('bay_segment', 'year')
  ) %>% 
  mutate(trndt = paste0('trndt', 1:n()), .by = 'bay_segment') %>% 
  select(-year)

##
# 

# thresholds to count by year, temp is above, salinity is below
thrs <- list(
  temp = c(29, 30, 31),
  sali = c(15, 20, 25)
)

# GAM models for salinity, temp by station with daily predictions for POR
load(file = here('data/saliprd.RData'))
load(file = here('data/tempprd.RData'))

# lookup table for threshold type by parameter
cntsel <- tibble(
  name = c('temp', 'sali'),
  tosel = c('cntabv', 'cntbel')
)

# get T/F vectors for predictions above/below thresholds by day
sumdat <- list(
  sali = saliprd,
  temp = tempprd
) %>%
  enframe() %>% 
  unnest('value') %>% 
  filter(param %in% c('salibot', 'tempbot')) %>% 
  full_join(cntsel, by = 'name') %>% 
  mutate(
    cnts = purrr::pmap(list(name, tosel, prd), function(name, tosel, prd){
      
      out <- thrs[[name]] %>%
        tibble(thr = .) %>%
        group_nest(thr) %>%
        mutate(data = list(prd)) %>%
        unnest('data') %>% 
        filter(yr < 2023) %>% 
        mutate(
          cntbel = value <= thr,
          cntabv = value >= thr
        )
      
      names(out)[names(out) %in% tosel] <- 'cnt'
      
      out <- out %>% 
        select(thr, date, doy, cont_year,  yr, value, cnt)
      
      return(out)
      
    })
  )

# get counts above/below thresholds by year, station, also counts for both thresholds
cmbdat <- sumdat %>% 
  select(-prd) %>% 
  unnest('cnts') %>% 
  nest(.by = c('name', 'thr')) %>% 
  unite(thr, c('name', 'thr'))
salicnt <- cmbdat %>% 
  filter(grepl('sali', thr))
tempcnt <- cmbdat %>% 
  filter(grepl('temp', thr))

cmbdat <- expand_grid(salicnt, tempcnt, .name_repair = 'unique')
names(cmbdat) <- c('salithr', 'salicnt', 'tempthr', 'tempcnt')
thrdat <- cmbdat %>% 
  mutate(
    cmbcnt = purrr::pmap(list(salicnt, tempcnt), function(salicnt, tempcnt){
      
      salicnt <- salicnt %>%
        select(bay_segment, station, date, doy, yr, salicnt = cnt)
      tempcnt <- tempcnt %>%
        select(station, date, tempcnt = cnt)
      out <- inner_join(salicnt, tempcnt, by = c('station', 'date')) %>%
        mutate(
          bothcnt = salicnt & tempcnt
        )
      
      return(out)
      
    })
  ) %>% 
  select(-salicnt, -tempcnt) %>% 
  unnest('cmbcnt') %>% 
  pivot_longer(names_to = 'thrtyp', values_to = 'cnt', matches('cnt')) %>% 
  filter(yr >= (min(year(trndts$date)) - 1)) %>% 
  filter(!bay_segment == 'LTB') %>% 
  arrange(bay_segment, station, thrtyp, salithr, tempthr, date)

thrtrndat <- thrdat %>% 
  left_join(trndts, by = c('bay_segment', 'date'), relationship = 'many-to-one') %>% 
  group_by(bay_segment, station, salithr, tempthr, thrtyp) %>% 
  fill(trndt, .direction = 'up') %>% 
  ungroup() %>% 
  filter(!is.na(trndt)) %>% # counts occurring after max date in trnds
  mutate(
    dycnt = rev(1:n()), 
    trndt = max(date),
    .by = c(bay_segment, station, thrtyp, salithr, tempthr, trndt)
  ) 