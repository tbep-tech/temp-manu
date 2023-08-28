library(wqtrends)
library(tbeptools)
library(tidyverse)
library(lubridate)
library(haven)
library(here)
library(sf)
library(lmerTest)
library(modelbased)

source(here('R/funcs.R'))

# temp GAMs for each station, save file -------------------------------------------------------

# save model files for each station, separate for bottom/surface temp
epcdata %>% 
  select(bay_segment, station = epchc_station, date = SampleTime, temptop = Temp_Water_Top_degC, tempbot = Temp_Water_Bottom_degC) %>% 
  filter(year(date) < 2023) %>% 
  filter(year(date) > 1975) %>% 
  pivot_longer(matches('temp'), names_to = 'param', values_to = 'value') %>% 
  mutate(
    stationdup = station,
    paramdup = param,
    date = as.Date(date), 
    cont_year = decimal_date(date), 
    yr = year(date), 
  ) %>% 
  unite('ind', bay_segment, station, param) %>% 
  group_nest(ind) %>% 
  mutate(
    cnt = 1:nrow(.)
  ) %>% 
  unite('ind', cnt, ind, sep = '-') %>% 
  deframe() %>% 
  iwalk(function(data, ind){
    
    cat(ind, '\n')
    
    yrcnt <- data %>% 
      na.omit() %>% 
      summarize(
        cnt = n(),
        .by = 'yr'
      ) %>% 
      arrange(yr) %>% 
      mutate(flt = cnt >= 5)
    yrmin <- yrcnt$yr[which(yrcnt$flt)[1]]
    
    mod <- data %>% 
      filter(yr >= yrmin) %>% 
      filter(!is.na(value)) %>% 
      anlz_gam(trans = 'ident')
    
    ind <- gsub('^.*-', '', ind)
    assign(ind, mod)
    fl <- here(paste0('data/', ind, '.RData'))
    save(list = ind, file = fl, compress = 'xz')
    
  })

# get daily temp predictions from GAM files ---------------------------------------------------

fls <- list.files(here('data'), pattern = 'temp', full.names = T)
obs <- gsub('\\.RData$', '', basename(fls))

tempprd <- tibble(obs = obs) %>% 
  mutate(fit = NA, prd = NA)
for(i in seq_along(fls)){
  
  cat(i, 'of', length(fls), '\n')
  
  fl <- fls[i]
  ob <- obs[i]
  
  load(file = fl)
  toadd <- get(ob)
  tempprd[tempprd$obs == ob, 'fit'][[1]] <- list(anlz_fit(toadd))
  tempprd[tempprd$obs == ob, 'prd'][[1]] <- list(anlz_prdday(toadd))
  
  rm(toadd)
  rm(list = ob)
  
}

tempprd <- tempprd %>% 
  separate(obs, c('bay_segment', 'station', 'param')) %>% 
  unnest('fit')

save(tempprd, file = here('data/tempprd.RData'))

# sali GAMs for each station, save file -------------------------------------------------------

# save model files for each station, separate for bottom/surface temp
epcdata %>% 
  select(bay_segment, station = epchc_station, date = SampleTime, salitop = Sal_Top_ppth, salibot = Sal_Bottom_ppth) %>% 
  filter(year(date) < 2023) %>% 
  filter(year(date) > 1975) %>% 
  pivot_longer(matches('sal'), names_to = 'param', values_to = 'value') %>% 
  mutate(
    stationdup = station,
    paramdup = param,
    date = as.Date(date), 
    cont_year = decimal_date(date), 
    yr = year(date), 
  ) %>% 
  unite('ind', bay_segment, station, param) %>% 
  group_nest(ind) %>% 
  mutate(
    cnt = 1:nrow(.)
  ) %>% 
  unite('ind', cnt, ind, sep = '-') %>% 
  deframe() %>% 
  iwalk(function(data, ind){
    
    cat(ind, '\n')
    
    yrcnt <- data %>% 
      na.omit() %>% 
      summarize(
        cnt = n(),
        .by = 'yr'
      ) %>% 
      arrange(yr) %>% 
      mutate(flt = cnt >= 5)
    yrmin <- yrcnt$yr[which(yrcnt$flt)[1]]
    
    mod <- data %>% 
      filter(yr >= yrmin) %>% 
      filter(!is.na(value)) %>% 
      anlz_gam(trans = 'ident')
    
    ind <- gsub('^.*-', '', ind)
    assign(ind, mod)
    fl <- here(paste0('data/', ind, '.RData'))
    save(list = ind, file = fl, compress = 'xz')
    
  })

# get daily sali predictions from GAM files ---------------------------------------------------

fls <- list.files(here('data'), pattern = 'sali', full.names = T)
obs <- gsub('\\.RData$', '', basename(fls))

saliprd <- tibble(obs = obs) %>% 
  mutate(fit = NA, prd = NA)
for(i in seq_along(fls)){
  
  cat(i, 'of', length(fls), '\n')
  
  fl <- fls[i]
  ob <- obs[i]
  
  load(file = fl)
  toadd <- get(ob)
  saliprd[saliprd$obs == ob, 'fit'][[1]] <- list(anlz_fit(toadd))
  saliprd[saliprd$obs == ob, 'prd'][[1]] <- list(anlz_prdday(toadd))
  
  rm(toadd)
  rm(list = ob)
  
}

saliprd <- saliprd %>% 
  separate(obs, c('bay_segment', 'station', 'param')) %>% 
  unnest('fit')

save(saliprd, file = here('data/saliprd.RData'))

# chla GAMs for each station, save file -------------------------------------------------------

# save model files for each station, separate for bottom/surface temp
epcdata %>% 
  select(bay_segment, station = epchc_station, date = SampleTime, chlatop = chla) %>% 
  filter(year(date) < 2023) %>% 
  filter(year(date) > 1975) %>% 
  pivot_longer(matches('chla'), names_to = 'param', values_to = 'value') %>% 
  mutate(
    stationdup = station,
    paramdup = param,
    date = as.Date(date), 
    cont_year = decimal_date(date), 
    yr = year(date), 
  ) %>% 
  unite('ind', bay_segment, station, param) %>% 
  group_nest(ind) %>% 
  mutate(
    cnt = 1:nrow(.)
  ) %>% 
  unite('ind', cnt, ind, sep = '-') %>% 
  deframe() %>% 
  iwalk(function(data, ind){
    
    cat(ind, '\n')
    
    yrcnt <- data %>% 
      na.omit() %>% 
      summarize(
        cnt = n(),
        .by = 'yr'
      ) %>% 
      arrange(yr) %>% 
      mutate(flt = cnt >= 5)
    yrmin <- yrcnt$yr[which(yrcnt$flt)[1]]

    mod <- data %>% 
      filter(yr >= yrmin) %>% 
      filter(!is.na(value)) %>% 
      anlz_gam(trans = 'log10')
    
    ind <- gsub('^.*-', '', ind)
    assign(ind, mod)
    fl <- here(paste0('data/', ind, '.RData'))
    save(list = ind, file = fl, compress = 'xz')
    
  })

# get daily chla predictions from GAM files ---------------------------------------------------

fls <- list.files(here('data'), pattern = 'chla', full.names = T)
obs <- gsub('\\.RData$', '', basename(fls))

chlaprd <- tibble(obs = obs) %>% 
  mutate(fit = NA, prd = NA)
for(i in seq_along(fls)){
  
  cat(i, 'of', length(fls), '\n')
  
  fl <- fls[i]
  ob <- obs[i]
  
  load(file = fl)
  toadd <- get(ob)
  chlaprd[chlaprd$obs == ob, 'fit'][[1]] <- list(anlz_fit(toadd))
  chlaprd[chlaprd$obs == ob, 'prd'][[1]] <- list(anlz_prdday(toadd))
  
  rm(toadd)
  rm(list = ob)
  
}

chlaprd <- chlaprd %>% 
  separate(obs, c('bay_segment', 'station', 'param')) %>% 
  unnest('fit')

save(chlaprd, file = here('data/chlaprd.RData'))

# count of days per year per station above/below thresholds -----------------------------------

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
  summarise(
    cnt = sum(cnt), 
    .by = c(bay_segment, station, salithr, tempthr, thrtyp, yr)
  )

save(thrdat, file = here('data/thrdat.RData'), compress = 'xz')

# mixef mods of threshold counts over time ----------------------------------------------------

load(file = here("data/thrdat.RData'"))

# create mixed effects models of cnt ~ yr with station as random intercept
mixmods <- thrdat %>% 
  group_by(bay_segment, salithr, tempthr, thrtyp) %>% 
  nest() %>% 
  summarise(
    mod = purrr::map(data, function(x){
      
      mod <- try(lmer(cnt ~ yr + (1|station), data = x, REML = F), silent = T)
      if(inherits(mod, 'try-error'))
        return(NA)
      
      effs <- estimate_means(mod, 'yr')
      
      summod <- summary(mod)$coefficients
      
      out <- tibble(
        slo = summod['yr', 'Estimate'],
        pvl = summod['yr', 'Pr(>|t|)'], 
        yrstr = effs$Mean[1],
        yrstrse = effs$SE[1], 
        yrend = effs$Mean[nrow(effs)],
        yrendse = effs$SE[nrow(effs)]
      )
      
      return(out)
      
    })
  ) %>% 
  unnest('mod') %>% 
  ungroup()

save(mixmods, file = here('data/mixmods.RData'), compress = 'xz')

# fim data ------------------------------------------------------------------------------------

prj <- 4326

# import original SAS data
phyraw <- read_sas('https://github.com/tbep-tech/fim-macroalgae/raw/main/data/raw/tbm_physical.sas7bdat')
hydraw <- read_sas('https://github.com/tbep-tech/fim-macroalgae/raw/main/data/raw/tbm_hydrolab.sas7bdat')

# physical data
# select AM projects (monhtly FIM sampling)
# zones A-E for TB proper, then clip by TB segments
# year >= 1998, but make sure to use only 2005 to present for trawls (300/301)
phydat <- phyraw %>% 
  mutate(date = ymd(date)) %>% 
  # filter(Project_1 == 'AM'| Project_2 == 'AM' | Project_3 == 'AM') %>% 
  filter(Zone %in% c('A', 'B', 'C', 'D', 'E')) %>% 
  filter(!is.na(Longitude) | !is.na(Latitude)) %>% 
  st_as_sf(coords = c('Longitude', 'Latitude'), crs = prj, remove = F) %>% 
  .[tbseg, ] %>% 
  select(date, Reference)

hyddat <- hydraw %>% 
  filter(Beg_end %in% 'B') %>% # each location has a beginning and end log, take one 
  filter(Depth == max(Depth), .by = 'Reference') %>% # some have depth profile, take bottom
  select(Reference, depth = Depth, temp = Temperature, sal = Salinity)

fimtempdat <- phydat %>% 
  inner_join(hyddat, by = 'Reference') %>% 
  filter(year(date) > 1995) %>% # only spring/fall sampling prior to 1996
  st_intersection(tbseg[, 'bay_segment']) 

save(fimtempdat, file = here('data/fimtempdat.RData'), compress = 'xz')
