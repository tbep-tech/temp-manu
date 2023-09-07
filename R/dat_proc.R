# setup ---------------------------------------------------------------------------------------

library(wqtrends)
library(tbeptools)
library(tidyverse)
library(lubridate)
library(haven)
library(here)
library(sf)
library(lmerTest)
library(modelbased)
library(mgcv)
library(SPEI)
library(rnoaa)
library(dataRetrieval)

source(here('R/funcs.R'))

# SPEI, monthly rain, and monthly air temp ----------------------------------------------------

noaa_key <- Sys.getenv('NOAA_KEY')

yrs <- seq(1975, 2022)

stameta <- ghcnd_stations()
sta <- 'USW00012842'
stalat <- stameta %>% 
  filter(id %in% sta) %>% 
  pull(latitude) %>% 
  unique()

res <- yrs %>%
  tibble::enframe('name', 'year') %>%
  dplyr::group_by(name) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    ests = purrr::map(data, function(x){
      
      yr <- x$year
      cat(yr, '\n')
      
      start <- paste0(yr, "-01-01")
      end <- paste0(yr, "-12-31")
      
      # download NOAA UWS rainfall and temperature station data
      tia_dat <- meteo_pull_monitors(monitors = sta, 
                                     date_min = start, date_max = end)
      
      out <- tia_dat
      
      return(out)
      
    })
  )

# tavg does not start until 1998, so reproduce from tmin/tmax
resdat <- res %>% 
  unnest('data') %>% 
  unnest('ests') %>% 
  ungroup() %>% 
  select(date, prcp, tmax, tmin) %>% 
  mutate(
    yr = year(date),
    mo = month(date),
    tavg_c = (tmax + tmin) / 2,
    tavg_c = tavg_c / 10,
    tavg_c = na.approx(tavg_c),
    precip_mm = prcp / 10, 
    precip_in = precip_mm * 0.0393701
  ) %>% 
  select(date, yr, mo, tavg_c, precip_mm, precip_in)

# monthly avg, sum
speidat <- resdat %>% 
  summarise(
    precip_mm = sum(precip_mm),
    precip_in = sum(precip_in), 
    tavg_c = mean(tavg_c),
    date = min(date),
    .by = c('yr', 'mo')
  ) %>% 
  mutate(
    PET = thornthwaite(tavg_c, stalat),
    BAL = precip_mm - PET,
    spei = as.numeric(spei(BAL, 12)$fitted), 
    spi = spi(precip_mm, 12)$fitted,
    speisign = ifelse(spei >= 0, 'wet', 'dry'), 
    spisign = ifelse(spi >= 0, 'wet', 'dry')
  )

save(speidat, file = here('data/speidat.RData'))

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
    cnt = runfunc(cnt), 
    .by = c(bay_segment, station, salithr, tempthr, thrtyp, yr)
  )

save(thrdat, file = here('data/thrdat.RData'), compress = 'xz')

# count of days per year per station above chl threshold --------------------------------------

# thresholds to count by year, segment specific
thrs <- targets %>% 
  filter(bay_segment %in% c('OTB', 'HB', 'MTB', 'LTB')) %>% 
  select(bay_segment, thr = chla_thresh)

# GAM models for chla by station with daily predictions for POR
load(file = here('data/chlaprd.RData'))

# get T/F vectors for predictions above thresholds by day
chlthrdat <- chlaprd %>% 
  left_join(thrs, by = 'bay_segment') %>% 
  mutate(
    cnts = purrr::pmap(list(thr, prd), function(thr, prd){

      out <- prd %>% 
        filter(yr < 2023) %>% 
        mutate(
          cnt = value >= thr
        ) %>% 
        summarise(
          cnt = runfunc(cnt), 
          .by = c('yr')
        )
      
      return(out)
      
    })
  ) %>% 
  select(bay_segment, station, chla_thr = thr, cnts) %>% 
  mutate(
    chla_thr = paste0('chla_', chla_thr), 
    thrtyp = 'chlacnt'
  ) %>% 
  unnest('cnts')

save(chlthrdat, file = here('data/chlthrdat.RData'), compress = 'xz')

# mixef mods of threshold counts over time ----------------------------------------------------

load(file = here("data/thrdat.RData"))

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

# combined transect, chl, sal, temp data ------------------------------------------------------

load(file = here('data/thrdat.RData'))
load(file = here('data/chlthrdat.RData'))

seglng <- c('Old Tampa Bay', 'Hillsborough Bay', 'Middle Tampa Bay', 'Lower Tampa Bay')
segshr <- c('OTB', 'HB', 'MTB', 'LTB')

##
# seagrass fo
transectocc <- anlz_transectocc(transect)
transectave <- anlz_transectave(transectocc)

fodat <- transectave %>% 
  select(bay_segment, yr, total = foest) %>% 
  filter(bay_segment %in% segshr) %>% 
  mutate(
    total = total / 100
  ) %>% 
  ungroup()

trnptsshed <- trnpts %>% 
  sf::st_set_geometry(NULL) %>% 
  select(Transect = TRAN_ID, bay_segment) %>% 
  unique()

bbdat <- transectocc %>% 
  ungroup() %>% 
  filter(Savspecies %in% c('Halodule', 'Syringodium', 'Thalassia')) %>% 
  filter(
    bbest > 0, 
    .by = c('Date', 'Transect')
  ) %>% 
  inner_join(trnptsshed, by = 'Transect') %>% 
  filter(!bay_segment %in% c('BCB', 'TCB', 'MR')) %>% 
  mutate(
    yr = year(Date)
  ) %>% 
  summarise(
    bbave = mean(bbest),
    .by = c('bay_segment', 'yr')
  )

# ##
# # coverage data
# 
# load(file = url('https://github.com/tbep-tech/tbep-os-presentations/raw/master/data/sgsegest.RData'))
# 
# sgsegest <- sgsegest %>% 
#   mutate(
#     segment = factor(segment, 
#                      levels = c("Old Tampa Bay", "Hillsborough Bay", "Middle Tampa Bay", "Lower Tampa Bay", 
#                                 "Boca Ciega Bay", "Terra Ceia Bay", "Manatee River"),
#                      labels = c('OTB', 'HB', 'MTB', 'LTB', 'BCB', 'TCB', 'MR'))
#   ) %>% 
#   filter(segment %in% segshr) %>% 
#   select(
#     bay_segment = segment,
#     yr = year, 
#     acres
#   )

##
# threshold count data

# these had strong models
salthr <- '25'
tmpthr <- '30'

cntsaltmp <- thrdat %>% 
  filter(salithr == paste0('sali_', salthr)) %>% 
  filter(tempthr == paste0('temp_', tmpthr)) %>% 
  summarise(
    cnt = mean(cnt),
    .by = c('thrtyp', 'bay_segment', 'yr')
  )

cntchla <- chlthrdat %>% 
  summarise(
    cnt = mean(cnt), 
    .by = c('thrtyp', 'bay_segment', 'yr')
  )

cntdat <- bind_rows(cntsaltmp, cntchla) %>% 
  pivot_wider(names_from = thrtyp, values_from = cnt)

##
# combine

cmbdat <- fodat %>% 
  inner_join(cntdat, by = c('yr', 'bay_segment')) %>% 
  inner_join(bbdat, by = c('yr', 'bay_segment')) %>% 
  mutate(
    bay_segment = factor(bay_segment, 
                         levels = segshr)
  ) %>% 
  rename(
    Year = yr,
    Sal = salicnt, 
    Temp = tempcnt, 
    Both = bothcnt, 
    Chla = chlacnt
  )

save(cmbdat, file = here('data/cmbdat.RData'))

# gam for cmbdat ------------------------------------------------------------------------------

load(file = here('data/cmbdat.RData'))

tomod <- cmbdat %>% 
  filter(!bay_segment %in% 'LTB')

cmbmod <- gam(total ~ ti(Year) + ti(Temp) + ti(Sal) + ti(Both) + ti(Temp, Year) + ti(Sal, Year) + ti(Both, Year), data = tomod)

# cmbmod <- gam(total ~ te(Temp, Year) + te(Sal, Year), data = tomod)

save(cmbmod, file = 'data/cmbmod.RData')


# USGS gage data ------------------------------------------------------------------------------

yrs <- seq(1975, 2022)

sites <- list(
  Hillsborough = "02303000",
  Alafia = "02301500",
  LittleManatee = "02300500",
  Tarpon = "02307498",
  Ward = "02300042",
  Manatee = "02299950"
)

res <- yrs %>%
  tibble::enframe('name', 'year') %>%
  dplyr::group_by(name) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    ests = purrr::map(data, function(x){
      
      yr <- x$year
      cat(yr, '\n')
      
      start <- paste0(yr, "-01-01")
      end <- paste0(yr, "-12-31")
      
      ##
      # flow 
      
      # download USGS streamflow data
      outflow <- lapply(sites, 
                        function(x) dataRetrieval::readNWISdv(x, "00060", start, end)
      )
      
      # combine flow
      # remove empty elements
      nodat <- lapply(outflow, function(x) nrow(x) == 0) %>% unlist
      if(!any(!nodat)) {
        outflow <- NULL
      } else {
        outflow <- outflow[!nodat]
        outflow <- outflow %>%
          enframe() %>% 
          unnest('value') %>% 
          rename(
            site = site_no,
            date = Date
          ) %>% 
          mutate(
            val = X_00060_00003 * 86400 / 35.3147 # ft3/s to m3/d
          ) %>% 
          summarise(val = sum(val, na.rm = T), .by = c(name, site, date)) %>% 
          mutate(
            var = 'flow_m3d'
          )
      }
      
      ##
      # water temp
      
      # download USGS temp data
      outtemp <- lapply(sites, 
                        function(x) readNWISdata(sites = x, service = "dv", parameterCd = "00010", startDate = start, endDate = end)
      )
      
      # remove empty list elements
      nodat <- lapply(outtemp, function(x) nrow(x) == 0) %>% unlist
      if(!any(!nodat)) {
        outtemp <- NULL
      } else {
        outtemp <- outtemp[!nodat]
        
        # create avg temp from min max 
        nms <- lapply(outtemp, 
                      function(x) 'X_00010_00001' %in% names(x) & 'X_00010_00002' %in% names(x)
        ) %>% unlist()
        if(any(nms)){
          outtemp[nms] <- lapply(outtemp[nms], function(x){
            
            x <- x %>% 
              mutate(
                X_00010_00011 = (X_00010_00001 + X_00010_00002) / 2
              ) %>% 
              select(-X_00010_00001, -X_00010_00002)
            
            return(x)
          })
        }
        
        outtemp <- outtemp %>% 
          enframe() %>% 
          unnest('value') %>% 
          rename(
            site = site_no, 
            date = dateTime,
            val = X_00010_00011
          ) %>%  
          summarise(val = mean(val, na.rm = T), .by = c(name, site, date)) %>% 
          mutate(
            var = 'temp_c'
          )
      }
      
      # combine all
      out <- bind_rows(outflow, outtemp)
      
      return(out)
      
    })
  )

gagedat <- res %>%
  unnest('data') %>%
  ungroup() %>% 
  select(-name) %>% 
  unnest('ests') %>%
  select(name, site, year, date, val, var)

save(gagedat, file = here('data/gagedat.RData'), compress = 'xz')

# ports water temp data -----------------------------------------------------------------------

noaa_key <- Sys.getenv('NOAA_KEY')

ports <- list(
  stpete = '8726520', 
  oldporttampa = '8726607', 
  manatee = '8726384', 
  mckay = '8726667'#, 
  # eastbay = '8726674' # only a few years of water temp data
)

yrs <- tibble(
  dts = seq.Date(from = as.Date('1975-01-01'), to = as.Date('2022-12-31'), 'days')
) %>% 
  mutate(
    mo = month(dts), 
    yr = year(dts)
  ) %>% 
  crossing(ports %>% enframe %>% unnest('value')) %>% 
  group_nest(yr, mo, name) %>% 
  mutate(
    res = purrr::pmap(list(name, yr, mo, data), function(name, yr, mo, data){
      
      cat(name, yr, mo)
      
      nm <- unique(data$value) %>% as.numeric()
      rng <- range(data$dts) %>% as.character() %>% gsub('\\-', '', .) %>% as.numeric()
      
      res <- try(coops_search(station_name = nm, begin_date = rng[1],
                              end_date = rng[2], product = "water_temperature"), silent = T)
      
      if(inherits(res, 'try-error')){
        cat('\tno data\n')
        return(NULL)
      }
      
      cat('\n')
      
      return(res$data)
      
    })
  )

portsdat <- yrs %>% 
  mutate(
    flt = purrr::map(res, is.data.frame)
  ) %>% 
  unnest(flt) %>% 
  filter(flt) %>% 
  select(yr, mo, name, res) %>% 
  unnest(res) %>% 
  filter(!is.na(v)) %>% 
  mutate(
    dy = day(t), 
  ) %>% 
  mutate(
    cnt = length(unique(hour(t))), 
    .by = c(name, yr, mo, dy)
  ) %>% 
  filter(cnt >= 20) %>% 
  summarise(
    temp_c = median(v), 
    date = as.Date(min(t)), 
    .by = c(name, yr, mo, dy)
  ) %>% 
  filter(!(yr < 2003 & name == 'stpete')) # wonky data before 2004

save(portsdat, file = here('data/portsdat.RData'))

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
