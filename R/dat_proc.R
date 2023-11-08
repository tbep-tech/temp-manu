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
library(readxl)

source(here('R/funcs.R'))
yrsel <- c(1998, 2022)

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


# gam performance, temp/sal only --------------------------------------------------------------

fls <- list.files(here('data'), 'salibot|tempbot')

gamfitnest <- fls %>% 
  tibble(fl = .) %>% 
  group_nest(fl) %>% 
  mutate(
    data = purrr::pmap(list(fl, data), function(fl, data){
      
      cat(fl, '\n')
      
      obj <- gsub('\\.RData', '', fl)
      load(file = paste0('data/', fl))
      mod <- get(obj)
      anlz_fit(mod)
      
    })
  )

gamfit <- gamfitnest %>% 
  unnest('data') %>% 
  separate(fl, into = c('bay_segment', 'epchc_station', 'parameter'), sep = '_', convert = T) %>% 
  left_join(stations, by = c('bay_segment', 'epchc_station')) %>% 
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')), 
    parameter = factor(parameter, levels = c('tempbot.RData', 'salibot.RData'), 
                       labels = c('temp', 'sal')
    ), 
    AIC = round(AIC, 0)
  ) %>% 
  mutate_at(vars(GCV, R2), round, 2) %>% 
  mutate_at(vars(Latitude, Longitude), round, 3) %>% 
  arrange(parameter, bay_segment, epchc_station) %>% 
  select(parameter, bay_segment, epchc_station, Longitude, Latitude, everything()) 

save(gamfit, file = here('data/gamfit.RData'))

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

# counts by day above/below thresholds referenced to transect dates ---------------------------

# note that these include all days and threhsolds
# data are not summarized by year as in thrdat

##
# approximate sample month of transects by year, segment

transectocc <- anlz_transectocc(transect) 
transectavespp <- transectocc %>% 
  anlz_transectavespp(by_seg = TRUE) %>% 
  filter(bay_segment %in% c('HB', 'OTB', 'MTB', 'LTB')) %>% 
  filter(Savspecies %in% 'total') %>% 
  select(yr, bay_segment, foest)

# mean month of tranect sampling in HB
trnptssub <- trnpts %>%
  st_set_geometry(NULL) %>%
  select(
    Transect = TRAN_ID,
    bay_segment
  ) %>%
  filter(bay_segment %in% c('OTB', 'HB', 'MTB', 'LTB'))
trndts <- transectocc %>%
  ungroup() %>% 
  filter(Savspecies == 'total') %>%
  inner_join(trnptssub, by = 'Transect') %>%
  mutate(
    year = year(Date)
  ) %>% 
  summarise(
    date = ymd(median(Date)), # use ymd because averaging dates leaves decimals 
    # avedate = ymd(mean(Date)),
    # mindate = ymd(min(Date)),
    # maxdate = ymd(max(Date)),
    .by = c('bay_segment', 'year')
  ) %>% 
  mutate(trndt = paste0('trndt', 1:n()), .by = 'bay_segment') %>% 
  select(-year)

##
# threshold counts by day

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
  arrange(bay_segment, station, thrtyp, salithr, tempthr, date)

##
# reference counts to transect sample dates by bay segment

thralltrndat <- thrdat %>% 
  left_join(trndts, by = c('bay_segment', 'date'), relationship = 'many-to-one') %>% 
  group_by(bay_segment, station, salithr, tempthr, thrtyp) %>% 
  fill(trndt, .direction = 'up') %>% 
  ungroup() %>% 
  filter(!is.na(trndt)) %>% # counts occurring after max date in trnds
  mutate(
    dycnt = rev(1:n()),
    trndtcnt = trndt,
    trndt = max(date),
    trnyr = max(year(date)),
    .by = c(bay_segment, station, thrtyp, salithr, tempthr, trndt)
  ) %>% 
  filter(!(dycnt > 365 & trndtcnt == 'trndt1')) %>%  # remove day cnts beyond one year for starting year
  select(-doy, -trndtcnt)

save(thralltrndat, file = here('data/thralltrndat.RData'))

# mixef mods of threshold counts over time ----------------------------------------------------

load(file = here("data/thrdat.RData"))

# create mixed effects models of cnt ~ yr with station as random intercept
mixmods <- thrdat %>% 
  group_by(bay_segment, tempthr, salithr, thrtyp) %>% 
  nest() %>% 
  summarise(
    mod = purrr::map(data, function(x){
      
      mod <- try(lmer(cnt ~ yr + (1|station), data = x, REML = F), silent = T)
      if(inherits(mod, 'try-error'))
        return(NA)
      
      effs <- estimate_means(mod, at = 'yr = c(1976, 2022)')
      
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

# supp mixef mods of threshold counts over time ----------------------------------------------------

load(file = here("data/thrdat.RData"))

# create mixed effects models of cnt ~ yr with station as random intercept
suppmixmods <- thrdat %>% 
  group_by(bay_segment, tempthr, salithr, thrtyp) %>% 
  nest() %>% 
  summarise(
    mod = purrr::map(data, function(x){

      tomod <- x %>% 
        filter(yr >= yrsel[1] & yr <= yrsel[2])
      
      mod <- try(lmer(cnt ~ yr + (1|station), data = tomod, REML = F), silent = T)
      if(inherits(mod, 'try-error'))
        return(NA)
      
      effs <- estimate_means(mod, at = 'yr = c(1998, 2022)')
    
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

save(suppmixmods, file = here('data/suppmixmods.RData'), compress = 'xz')

# mixef mod predictions for select thresolds --------------------------------------------------

load(file = here('data/thrdat.RData'))

salthr <- '25'
tmpthr <- '30'

mixmodprds <- thrdat %>% 
  filter(salithr == paste0('sali_', salthr)) %>% 
  filter(tempthr == paste0('temp_', tmpthr)) %>% 
  group_by(thrtyp, bay_segment) %>% 
  nest() %>% 
  mutate(
    mod = purrr::map(data, function(data){
      
      mod <- try(lmer(cnt ~ yr + (1|station), data = data, REML = F), silent = T)
      if(inherits(mod, 'try-error'))
        return(NULL)
      
      return(mod)
      
    }),
    data = purrr::pmap(list(data, mod), function(data, mod){
      if(!is.null(mod))
        out <- bind_cols(data, prd = predict(mod))
      else
        out <- data %>% 
          mutate(prd = NA)
      return(out)
    }),
    fix = purrr::map(mod, function(mod){
      
      if(!is.null(mod)){
      
        yrpd <- sort(unique(mod@frame$yr))
        atv <- paste0('yr=c(', paste(yrpd, collapse = ', '), ')')
        fixef <- estimate_means(mod, at = atv)
        
        out <- tibble(
          yr = fixef$yr,
          prd = fixef$Mean
        )
        
        return(out)
        
      }
      
    }),
    slo = purrr::map(mod, function(mod){
      
      if(!is.null(mod)){
        
        summod <- summary(mod)$coefficients
        
        pvl <- summod['yr', 'Pr(>|t|)']
        if(pvl >= 0.05)
          return('')
        
        out <- summod['yr', 'Estimate'] %>%
          round(2) %>%
          as.character()
        
        return(out)
        
      }
      
    })
  ) %>% 
  ungroup() %>% 
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')), 
    thrtyp = factor(thrtyp, 
                    levels = c('tempcnt', 'salicnt', 'bothcnt'), 
                    labels = c(paste('Temperature >', tmpthr), paste('Salinity <', salthr), 'Both')
    )
  )

save(mixmodprds, file = here('data/mixmodprds.RData'))

# supp mixef mod predictions for select thresolds ---------------------------------------------

load(file = here('data/thrdat.RData'))

salthr <- '25'
tmpthr <- '30'

suppmixmodprds <- thrdat %>% 
  filter(yr >= yrsel[1] & yr <= yrsel[2]) %>% 
  filter(salithr == paste0('sali_', salthr)) %>% 
  filter(tempthr == paste0('temp_', tmpthr)) %>% 
  group_by(thrtyp, bay_segment) %>% 
  nest() %>% 
  mutate(
    mod = purrr::map(data, function(data){
      
      mod <- try(lmer(cnt ~ yr + (1|station), data = data, REML = F), silent = T)
      if(inherits(mod, 'try-error'))
        return(NULL)
      
      return(mod)
      
    }),
    data = purrr::pmap(list(data, mod), function(data, mod){
      if(!is.null(mod))
        out <- bind_cols(data, prd = predict(mod))
      else
        out <- data %>% 
          mutate(prd = NA)
      return(out)
    }),
    fix = purrr::map(mod, function(mod){
      
      if(!is.null(mod)){
        
        yrpd <- sort(unique(mod@frame$yr))
        atv <- paste0('yr=c(', paste(yrpd, collapse = ', '), ')')
        fixef <- estimate_means(mod, at = atv)
        
        out <- tibble(
          yr = fixef$yr,
          prd = fixef$Mean
        )
        
        return(out)
        
      }
      
    }),
    slo = purrr::map(mod, function(mod){
      
      if(!is.null(mod)){
        
        summod <- summary(mod)$coefficients
        
        pvl <- summod['yr', 'Pr(>|t|)']
        if(pvl >= 0.05)
          return('')
        
        out <- summod['yr', 'Estimate'] %>%
          round(2) %>%
          as.character()
        
        return(out)
        
      }
      
    })
  ) %>% 
  ungroup() %>% 
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')), 
    thrtyp = factor(thrtyp, 
                    levels = c('tempcnt', 'salicnt', 'bothcnt'), 
                    labels = c(paste('Temperature >', tmpthr), paste('Salinity <', salthr), 'Both')
    )
  )

save(suppmixmodprds, file = here('data/suppmixmodprds.RData'))

# combined epc metric, chl, sal, temp data ----------------------------------------------------

load(file = here('data/thralltrndat.RData'))

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

##
# threshold count data

# these had strong models
salthr <- '25'
tmpthr <- '30'

# dys <- seq(30, 365, by = 30)
# for(i in dys){
#   
#   thrtrndat <- thralltrndat %>% 
#     filter(salithr == paste0('sali_', salthr)) %>% 
#     filter(tempthr == paste0('temp_', tmpthr)) %>% 
#     filter(dycnt <= i) %>% 
#     summarise(
#       cnt = runfunc(cnt),
#       .by = c('bay_segment', 'station', 'thrtyp', 'trnyr')
#     ) %>% 
#     summarise(
#       cnt = mean(cnt), 
#       .by = c('bay_segment', 'thrtyp', 'trnyr')
#     ) %>% 
#     rename(yr = trnyr)
#   
#   ##
#   # combine
#   
#   cmbdat <- fodat %>% 
#     inner_join(thrtrndat, by = c('yr', 'bay_segment')) %>% 
#     inner_join(bbdat, by = c('yr', 'bay_segment')) %>% 
#     mutate(
#       bay_segment = factor(bay_segment, 
#                            levels = segshr)
#     ) %>% 
#     pivot_wider(names_from = 'thrtyp', values_from = 'cnt') %>% 
#     rename(
#       Year = yr,
#       Sal = salicnt, 
#       Temp = tempcnt, 
#       Both = bothcnt
#     ) %>% 
#     arrange(bay_segment, Year) %>% 
#     mutate(
#       chng = c(NA, sign(diff(total))), 
#       chng = ifelse(chng == -1, 1, 0), 
#       .by = 'bay_segment'
#     ) %>% 
#     filter(!is.na(chng))
#   
#   tomod <- cmbdat %>% 
#     mutate(
#       pchg = (total - lag(total)) / lag(total), 
#       chg = c(NA, sign(diff(total))), 
#       chg = ifelse(chg == -1, 1, 0),
#       .by = bay_segment, 
#       yrcat = cut(Year, breaks = c(-Inf, 2016, Inf), labels = c('Recovery (pre - 2016)', 'Decline (2016 - present)'), right = F)
#     ) %>% 
#     filter(!is.na(pchg)) %>% 
#     filter(bay_segment != 'LTB')
#   
#   combmod <- glm(chg ~ bay_segment*Sal*yrcat + bay_segment*Temp*yrcat, data = tomod, family = binomial('logit')) %>% 
#     step(trace = 0)
#   
#   cat(i, '\n')
#   # print(with(summary(combmod), 1 - deviance/null.deviance))
#   print(coefficients(summary(combmod)))
#   
# }

# use i = 360

##
# % change

i <- 360

thrtrndat <- thralltrndat %>% 
  filter(salithr == paste0('sali_', salthr)) %>% 
  filter(tempthr == paste0('temp_', tmpthr)) %>% 
  filter(dycnt <= i) %>% 
  summarise(
    cnt = runfunc(cnt),
    .by = c('bay_segment', 'station', 'thrtyp', 'trnyr')
  ) %>% 
  summarise(
    cnt = mean(cnt), 
    .by = c('bay_segment', 'thrtyp', 'trnyr')
  ) %>% 
  rename(yr = trnyr)

epccmbdat <- fodat %>% 
  inner_join(thrtrndat, by = c('yr', 'bay_segment')) %>% 
  inner_join(bbdat, by = c('yr', 'bay_segment')) %>% 
  mutate(
    bay_segment = factor(bay_segment, 
                         levels = segshr)
  ) %>% 
  pivot_wider(names_from = 'thrtyp', values_from = 'cnt') %>% 
  rename(
    sal = salicnt, 
    temp = tempcnt, 
    both = bothcnt
  ) %>% 
  arrange(bay_segment, yr)

save(epccmbdat, file = here('data/epccmbdat.RData'))

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
habraw <- read_sas('https://github.com/tbep-tech/fim-macroalgae/raw/main/data/raw/tbm_habitat.sas7bdat')
fimraw <- read_sas('https://github.com/tbep-tech/fim-macroalgae/raw/main/data/raw/fim_codes.sas7bdat')

# codes from fimraw, these are inclusive of bottom veg, not bycatch, latter would include more floating macro 
sgcode <- c('GM', 'GU', 'HA', 'HC', 'HE', 'HI', 'HJ', 'RU', 'SY', 'TH')
mccode <- c("AB", "AC", "AG", "AM", "AR", "AT", "AU", "BA", "CA", 
            "GR", "HM", "LA")

# physical data
# filter zones A-E for TB proper, then clip by TB segments
# filter gear type 20 (nearshore seine)
# filter referenc with loc info

phydat <- phyraw %>% 
  mutate(
    date = ymd(date), 
    starttime = as.character(starttime),
    hr = gsub('^([[:digit:]]+):.*$', '\\1', starttime),
    hr = as.numeric(hr)
  ) %>% 
  filter(Zone %in% c('A', 'B', 'C', 'D', 'E')) %>% 
  filter(Gear == 20) %>% 
  filter(!is.na(Longitude) | !is.na(Latitude)) %>% 
  filter(!is.na(BottomVegCover)) %>% 
  st_as_sf(coords = c('Longitude', 'Latitude'), crs = prj, remove = F) %>% 
  .[tbseg, ] %>% 
  select(Reference, date, hr, BottomVegCover)

# sal, temp data (no location)
hyddat <- hydraw %>% 
  filter(Beg_end %in% 'B') %>% # each location has a beginning and end log, take one 
  filter(Depth == max(Depth), .by = 'Reference') %>% # some have depth profile, take bottom
  select(Reference, depth = Depth, temp = Temperature, sal = Salinity)

# combine sal, temp data with location data, specific to tb
fimtempdat <- phydat %>% 
  inner_join(hyddat, by = 'Reference') %>% 
  filter(year(date) > 1995) %>% # only spring/fall sampling prior to 1996
  st_intersection(tbseg[, 'bay_segment']) 

# habitat data, summarized as sg present/abs
# must multiple total bottom veg cover from hydat by bottom veg ratio to get true cover
# sites with > 50% sg are T, otherwise F, none is truly no seagrass
# sites with only macro are not included
habdat <- habraw %>% 
  filter(Reference %in% fimtempdat$Reference) %>% 
  mutate(
    bvegcat = case_when(
      BottomVeg %in% sgcode ~ 'sg', # seagrass codes
      BottomVeg %in% mccode ~ 'mc', # mc codes
      BottomVeg %in% c('NO', 'N') ~ 'none'
    ) 
  ) %>%
  select(Reference, BottomVeg, BottomVegRatio, bvegcat) %>% 
  filter(!is.na(bvegcat) | bvegcat == 'mc') %>% 
  left_join(phydat, by = 'Reference') %>% 
  select(-geometry) %>% 
  mutate(
    vegcov = BottomVegRatio * BottomVegCover / 10
  ) %>% 
  summarise(
    sgpres = case_when(
      all(is.na(vegcov)) ~ 0, 
      sum(vegcov, na.rm = T) > 50 ~ 1,
      sum(vegcov, na.rm = T) <= 50 ~ 0
    ),
    sgcov = sum(vegcov, na.rm = T),
    .by = Reference
  )

# combine seagrass p/a with fimtempdat
fimsgtempdat <- fimtempdat %>% 
  left_join(habdat, by = 'Reference') %>% 
  filter(!is.na(sgpres)) %>% # will still include those with algae only
  mutate(
    lon = st_coordinates(.)[, 1], 
    lat = st_coordinates(.)[, 2], 
    yr = year(date)
  ) %>% 
  st_set_geometry(NULL)

save(fimsgtempdat, file = here('data/fimsgtempdat.RData'), compress = 'xz')

# pinellas data -------------------------------------------------------------------------------

# all otb samples
# four sites (polygons) E1:E4, 32 samples per year (eight per site)
tmp <- read_excel(here('data-raw/pinellas_from_sd.xlsx'), 
                  col_types = 'text')

pincotemp <- tmp %>% 
  select(
    site = Site, 
    sample = Sample,
    date = Date,
    hr = Time, 
    lat = Latitude, 
    lon = Longitude, 
    level = Level, 
    temp  = Temp_Water, 
    sal = Salinity, 
    SAV, 
    SAV_Type
  ) %>%
  mutate(
    date = mdy(date),
    hr = gsub('\\:.*$', '', hr), 
    hr = case_when(
      as.numeric(hr) %% 1 != 0 ~ round(as.numeric(hr) * 24),
      T ~ as.numeric(hr)
    ),
    SAV = case_when(
      SAV %in% c('n', 'N') ~ 0, 
      SAV %in% c('T', 'Y') ~ 1, 
      T ~ NaN
    ), 
    level = factor(level, levels = c('Surface', 'Middle', 'Bottom')), 
    levelint = as.numeric(level), 
    lon = ifelse(sign(as.numeric(lon)) == 1, -1 * as.numeric(lon), as.numeric(lon))
  ) %>% 
  filter(SAV %in% c(0, 1)) %>% 
  mutate_at(vars(lat, lon, temp, sal), as.numeric) %>% 
  filter(lat < 35) # one outlier at 72

# get salinity, temperature by lowest level (not always bottom)
saltemp <- pincotemp %>% 
  select(-SAV, -SAV_Type) %>% 
  mutate(
    levelmax = max(levelint), 
    .by = c(site, sample, date)
  ) %>% 
  filter(levelint == levelmax) %>% 
  select(-levelint, -levelmax)

# identify sites w/ and w/o seagrass (sav 1 could also include macro)
# allsg sites are those with only sg species found (will be zero if any macro found)
sav <- pincotemp %>% 
  select(-temp, -sal, -level, -levelint) %>% 
  unique() %>% 
  mutate(SAV_Type = strsplit(SAV_Type, ',|,\\s')) %>% 
  unnest('SAV_Type') %>% 
  mutate(
    allsg = case_when(
      any(!SAV_Type %in% c('H', 'Halophila', 'Hd', 'Hs', 'Hw', 'HW', 'R', 'Rm', 'RM', 'S', 'Sf', 'SF', 'T', 'Tt')) ~ 0, 
      all(is.na(SAV_Type)) ~ 0,
      T ~ 1
      ), 
    .by = c('site', 'sample', 'date')
  ) %>% 
  select(-SAV, -SAV_Type) %>% 
  unique()

# rejoin temp, sal with sav ifno
pincotemp <-  saltemp %>%
  inner_join(sav, by = c('site', 'sample', 'date', 'hr', 'lat', 'lon')) %>% 
  mutate(
    yr = year(date), 
    mo = month(date), 
    bay_segment = 'OTB'
  ) %>% 
  filter(yr > 2003) # only a few
  
save(pincotemp, file = here('data/pincotemp.RData'))

# seagrass declines models, epc, fim, and pinco -----------------------------------------------

load(file  = here('data/epccmbdat.RData'))
load(file = here('data/fimsgtempdat.RData'))
load(file = here('data/pincotemp.RData'))

modprep <- function(dat){

  dat %>% 
    filter(!bay_segment %in% c('LTB')) %>% 
    mutate(
      yrcat = case_when(
        yr <= 2016 ~ 'pre', 
        yr > 2016 ~ 'post'
      ), 
      yrcat = factor(yrcat, levels = c('pre', 'post'))
    )

}

epctomod <- modprep(epccmbdat)
fimtomod <- modprep(fimsgtempdat)
pincotomod <- modprep(pincotemp)

epcmod <- lm(total ~ sal * temp * yrcat + bay_segment, data = epctomod) %>% 
  step()
fimmod <- lm(sgcov ~ sal * temp * yrcat + bay_segment, data = fimtomod) %>% 
  step()
pincomod <- glm(allsg ~ sal * temp * yrcat, data = pincotomod, family = binomial(link = 'logit')) %>%
  step()

sgmods <- list(epcmod = epcmod, fimmod = fimmod, pincomo = pincomod)

save(sgmods, file = here('data/sgmods.RData'))

