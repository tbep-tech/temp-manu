library(wqtrends)
library(tbeptools)
library(tidyverse)
library(lubridate)
library(haven)
library(here)
library(sf)

# GAMs for each station, save file ------------------------------------------------------------

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
      anlz_gam(trans = 'ident')
    
    ind <- gsub('^.*-', '', ind)
    assign(ind, mod)
    fl <- here(paste0('data/', ind, '.RData'))
    save(list = ind, file = fl, compress = 'xz')
    
  })

# get daily predictions from GAM files --------------------------------------------------------

fls <- list.files(here('data'), pattern = 'temp', full.names = T)
obs <- gsub('\\.RData$', '', basename(fls))

moddat <- tibble(obs = obs) %>% 
  mutate(fit = NA, prd = NA)
for(i in seq_along(fls)){
  
  cat(i, 'of', length(fls), '\n')
  
  fl <- fls[i]
  ob <- obs[i]
  
  load(file = fl)
  toadd <- get(ob)
  moddat[moddat$obs == ob, 'fit'][[1]] <- list(anlz_fit(toadd))
  moddat[moddat$obs == ob, 'prd'][[1]] <- list(anlz_prdday(toadd))
  
  rm(toadd)
  rm(list = ob)
  
}

moddat <- moddat %>% 
  separate(obs, c('bay_segment', 'station', 'param')) %>% 
  unnest('fit')

save(moddat, file = here('data/moddat.RData'))

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
  filter(Beg_end %in% 'B') %>% 
  filter(Depth == max(Depth), .by = 'Reference') %>% 
  select(Reference, depth = Depth, temp = Temperature, sal = Salinity)

fimdat <- phydat %>% 
  inner_join(hyddat, by = 'Reference') %>% 
  filter(year(date) > 1995) %>% # only spring/fall sampling prior to 1996
  st_intersection(tbseg[, 'bay_segment']) 

save(fimdat, file = here('data/fimdat.RData'), compress = 'xz')
