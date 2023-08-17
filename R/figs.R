# setup ---------------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(EnvStats)
library(tbeptools)
library(ggmap)
library(sf)
library(patchwork)
library(here)
library(wqtrends)

yrsel <- c(1974, 2022)
tempthr <- c(29, 30, 31)
tempcol <- c('coral', 'red2', 'darkred')
salithr <- c(20, 25, 30)
salicol <- c('navyblue', 'dodgerblue2', 'slategray3')

# seagrass loss -------------------------------------------------------------------------------

load(file = url('https://github.com/tbep-tech/tbep-os-presentations/raw/master/data/sgsegest.RData'))

sgsegest <- sgsegest %>% 
  mutate(
    segment = factor(segment, 
                     levels = c("Old Tampa Bay", "Hillsborough Bay", "Middle Tampa Bay", "Lower Tampa Bay", 
                                "Boca Ciega Bay", "Terra Ceia Bay", "Manatee River"),
                     labels = c('OTB', 'HB', 'MTB', 'LTB', 'BCB', 'TCB', 'MR'))
  )

# segment coverage targets in 1k acres
segtrgs <- tibble(
  segment = factor(c(levels(sgsegest$segment), 'Total')), 
  trgs = c(11.1, 1.751, 9.4, 7.4, 8.8, 1.1, 0.449, 40)
)  

savval <- c('Halodule', 'Syringodium', 'Thalassia')
savleg <- expression(italic('H. wrightii'), italic('S. filiforme'), italic('T. testudinum'))
savcol <- c('#ED90A4', '#CCA65A', '#7EBA68')
names(savcol) <- savval

toplo1 <- sgsegest %>%
  filter(!segment %in% c('BCB', 'TCB', 'MR')) %>%
  mutate(acres = acres / 1000) %>%
  mutate(segment = forcats::fct_drop(segment))

subsegtrgs <- segtrgs %>%
  filter(segment %in% levels(toplo1$segment))

p1 <- ggplot(toplo1, aes(x = factor(year), y = acres)) +
  geom_bar(fill = '#00806E', stat = 'identity', colour = 'black', width = 0.6) +
  geom_hline(data = subsegtrgs, aes(yintercept = trgs, color = 'Target')) +
  scale_color_manual(values = 'red') +
  facet_wrap(~segment, ncol = 4) +
  theme_bw() + 
  theme(panel.grid.minor =element_blank(),
        panel.grid.major.x =element_blank(),
        # plot.background = element_rect(fill = NA, color = NA),
        axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(colour = 'black', angle = 45, size = 6, hjust = 1),
        strip.background = element_blank(),
        # strip.text = element_text(hjust = 0, size = 13),
        legend.position = 'none'
  ) +
  labs(
    y = 'Seagrass Coverage (x1,000 acres)',
    x = NULL,
    color = NULL,
    subtitle = '(a) Coverage changes for all species by bay segment',
    caption = expression(italic('Source: Southwest Florida Water Management District'))
  )

transectocc <- anlz_transectocc(transect)
transectavespp <- anlz_transectavespp(transectocc, by_seg = TRUE)

toplo2 <- transectave %>% 
  filter(bay_segment %in% c('OTB', 'HB', 'MTB', 'LTB')) %>% 
  filter(Savspecies %in% c('Halodule', 'Syringodium', 'Thalassia'))

p2 <- ggplot(toplo2, aes(x = yr, y = foest, color = Savspecies)) + 
  geom_line() + 
  geom_point() + 
  scale_x_continuous(breaks = seq(min(toplo2$yr), max(toplo2$yr), by = 2)) +
  scale_color_manual(
    values = savcol, 
    labels = savleg
  ) +
  facet_wrap(~bay_segment, ncol = 4) + 
  theme_bw() + 
  theme(
    legend.position = 'top',
    strip.background = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(colour = 'black', angle = 45, size = 6, hjust = 1),
  ) + 
  labs(
    y = 'Frequency Occurrence', 
    x = NULL, 
    color = NULL,
    subtitle = '(b) Frequency occurrence changes by species by bay segment',
    caption = expression(italic('Source: Interagency Seagrass Monitoring Program'))
  )

p <- p1 + p2 + plot_layout(ncol = 1)

png(here('figs/seagrasschg.png'), height = 6, width = 8, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# kendall all months --------------------------------------------------------------------------

sktres <- epcdata %>% 
  select(bay_segment, station = epchc_station, SampleTime, lon = Longitude, lat = Latitude, yr, 
         mo, matches('Top|Bottom')) %>% 
  filter(yr >= yrsel[1] & yr <= yrsel[2]) %>% 
  pivot_longer(matches('Top|Bottom'), names_to = 'var', values_to = 'val') %>% 
  nest(.by = c('bay_segment', 'var', 'station', 'lon', 'lat')) %>% 
  mutate(
    skt = purrr::pmap(list(station, var, data), function(station, var, data){
      
      cat(station, var, '\n')
      
      # yr selection
      tomod <- data %>% 
        arrange(yr, mo) %>% 
        na.omit()
      
      # run tests
      ests <- kendallSeasonalTrendTest(val ~ mo + yr, data = tomod)
      
      out <- tibble(
        pval = ests$p.value[2][[1]], 
        slos = ests$estimate[2][[1]],
        n = nrow(tomod)
      )
      
      return(out)
      
    })
  ) %>% 
  select(-data) %>% 
  unnest(skt)

# basemap
dat_ext <- unname(st_bbox(st_buffer(tbseg, dist = 2000)))
bsmap1 <- get_stamenmap(bbox = dat_ext, maptype = 'terrain-background', zoom = 11)

# change opacity of basemap
mapatt <- attributes(bsmap1)
bsmap1_transparent <- matrix(adjustcolor(bsmap1,
                                         alpha.f = 0.4),
                             nrow = nrow(bsmap1))
attributes(bsmap1_transparent) <- mapatt

leglab <- expression(paste(yr^{-1}))
colrng <- range(sktres$slos) %>% 
  abs %>% 
  max
colrng <- c(-1 * colrng, colrng)

pthm <- theme_bw(base_family = 'serif', base_size = 14) +
  theme(
    legend.position = 'right',
    # legend.box = 'vertical',
    strip.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 14)
  )

toplo <- sktres %>%
  mutate(
    pvalcol = ifelse(pval < 0.05, T, F),
    coefsgn = sign(slos),
    coefsgn = factor(coefsgn, levels = c('1', '-1'), labels = c('inc', 'dec')), 
    loc =  gsub("^.*_(.*)_.*$", "\\1", var),
    var = gsub('^(.*?)_.*$', '\\1', var, perl = T),
    loc = factor(loc, levels = c('Top', 'Bottom'))
  )

p <- ggmap(bsmap1_transparent) +
  geom_point(data = toplo, aes(x = lon, y = lat, size = abs(slos), fill = slos, shape = coefsgn, color = pvalcol), stroke = 1) +
  facet_grid(loc ~ var) +
  # scale_fill_gradient2(leglab, low = 'blue', mid = 'grey',  high = 'tomato1', midpoint = 0) +
  scale_fill_gradientn(leglab, limits = colrng, colors = c('blue', 'grey', 'tomato1')) +
  scale_color_manual(values = c(scales::alpha('black', 0), 'black'), guide = F, drop = FALSE) +
  coord_map() +
  scale_shape_manual(values = c(24, 25), guide = 'none', drop = FALSE) +
  pthm +
  scale_size(range = c(2, 6), guide = F) +
  labs(
    caption = 'Outlines indicate p < 0.05'
  )

png(here('figs/sktall.png'), height = 7, width = 6, family = 'serif', units = 'in', res = 300)
p
dev.off()

# kendall by month ----------------------------------------------------------------------------

ktres <- epcdata %>% 
  select(bay_segment, station = epchc_station, SampleTime, lon = Longitude, lat = Latitude, yr, 
         mo, matches('Top|Bottom')) %>% 
  filter(yr >= yrsel[1] & yr <= yrsel[2]) %>% 
  pivot_longer(matches('Top|Bottom'), names_to = 'var', values_to = 'val') %>% 
  nest(.by = c('bay_segment', 'var', 'station', 'lon', 'lat', 'mo')) %>% 
  mutate(
    kt = purrr::pmap(list(station, var, mo, data), function(station, var, mo, data){
      
      cat(station, var, mo, '\n')
      
      # yr selection
      tomod <- data %>% 
        arrange(yr) %>% 
        na.omit()
      
      out <- tibble(
        pval = NA,
        slos = NA,
        n = NA
      )
      
      # run tests
      ests <- try(kendallTrendTest(val ~ yr, data = tomod), silent = T)
      if(inherits(ests, 'try-error'))
        return(out)
      
      out$pval <- ests$p.value[1]
      out$slos <- ests$estimate[2]
      out$n <- nrow(tomod)
      
      return(out)
      
    })
  ) %>% 
  select(-data) %>% 
  unnest(kt)

# basemap
dat_ext <- unname(st_bbox(st_buffer(tbseg, dist = 2000)))
bsmap1 <- get_stamenmap(bbox = dat_ext, maptype = 'terrain-background', zoom = 11)

# change opacity of basemap
mapatt <- attributes(bsmap1)
bsmap1_transparent <- matrix(adjustcolor(bsmap1,
                                         alpha.f = 0.4),
                             nrow = nrow(bsmap1))
attributes(bsmap1_transparent) <- mapatt

leglab <- expression(paste(yr^{-1}))
colrng <- range(ktres$slos) %>% 
  abs %>% 
  max
colrng <- c(-1 * colrng, colrng)

pthm <- theme_bw(base_family = 'serif', base_size = 14) +
  theme(
    legend.position = 'right',
    # legend.box = 'vertical',
    strip.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 14)
  )

ktres %>%
  mutate(
    pvalcol = ifelse(pval < 0.05, T, F),
    coefsgn = sign(slos),
    coefsgn = factor(coefsgn, levels = c('1', '-1'), labels = c('inc', 'dec')), 
    mo = sprintf('%02d', mo), 
    sz = scales::rescale(abs(slos), c(2, 6)),
    loc =  gsub("^.*_(.*)_.*$", "\\1", var),
    var = gsub('^(.*?)_.*$', '\\1', var, perl = T),
    loc = factor(loc, levels = c('Top', 'Bottom'))
  ) %>% 
  group_nest(mo) %>% 
  deframe() %>% 
  iwalk(function(x, mo){

    p <- ggmap(bsmap1_transparent) +
      geom_point(data = x, aes(x = lon, y = lat, fill = slos, shape = coefsgn, color = pvalcol), stroke = 1, 
                 size = x$sz) +
      facet_grid(loc ~ var) +
      scale_fill_gradientn(leglab, limits = colrng, colors = c('blue', 'grey', 'tomato1')) +
      # scale_fill_gradient2(leglab, low = 'blue', mid = 'grey',  high = 'tomato1', midpoint = 0) +
      scale_color_manual(values = c(scales::alpha('black', 0), 'black'), guide = F, drop = FALSE) +
      coord_map() +
      scale_shape_manual(values = c(24, 25), guide = 'none', drop = FALSE) +
      pthm +
      # scale_size(range = c(1, 6), guide = F) +
      labs(
        caption = 'Outlines indicate p < 0.05'
      )
    
    fl <- paste0('figs/skt_',sprintf('%02d', as.numeric(mo)), '.png')
    png(here(fl), height = 7, width = 6, family = 'serif', units = 'in', res = 300)
    print(p)
    dev.off()
    
  })

# segment summary of monthly kendall by temp loc ----------------------------------------------

ktres %>% 
  mutate(
    mo = month(mo, label = T), 
    bay_segment = factor(bay_segment, levels = c('LTB', 'MTB', 'HB', 'OTB')), 
    loc =  gsub("^.*_(.*)_.*$", "\\1", var),
    varsimp = gsub('^(.*?)_.*$', '\\1', var, perl = T),
    loc = factor(loc, levels = c('Top', 'Bottom'))
  ) %>%
  nest(.by = c('bay_segment', 'mo', 'varsimp', 'var', 'loc')) %>% 
  mutate(
    sum = purrr::pmap(list(varsimp, data), function(varsimp, data){

      data %>% 
        summarise(
          nsig = case_when(
            varsimp == 'Sal' ~ -1 * sum(pval < 0.05 & slos < 0),
            varsimp == 'Temp' ~ sum(pval < 0.05 & slos > 0)
          ),
          cnt = n(),
          nsigper = round(100 * nsig / cnt, 0),
          perlab = ifelse(nsigper == 0, '', as.character(abs(nsigper)))
        )
      
    })
  ) %>% 
  select(-data) %>% 
  unnest('sum') %>% 
  group_nest(loc) %>% 
  deframe %>% 
  iwalk(function(x, loc){

    p <- ggplot(x, aes(x = mo, y = bay_segment, fill = nsigper)) + 
      geom_tile(color = 'darkgrey') +
      geom_text(aes(label = perlab)) +
      scale_x_discrete(expand = c(0, 0)) + 
      scale_y_discrete(expand = c(0, 0)) + 
      scale_fill_gradientn(leglab, limits = c(-100, 100), colors = c('blue', 'grey', 'tomato1')) +
      facet_wrap(~ varsimp, ncol = 1, scales = 'free_x') + 
      theme(
        strip.background = element_blank(), 
        strip.text = element_text(hjust = 0, size = 12), 
        legend.position = 'none'
      ) + 
      labs(
        x = NULL, 
        y = 'Bay segment'
      )
    
    fl <- paste0('figs/segsum_', loc, '.png')
    png(here(fl), height = 4, width = 6, family = 'serif', units = 'in', res = 300)
    print(p)
    dev.off()
    
  })

# temp station gam --------------------------------------------------------------------------

fl <- 'OTB_66_tempbot'
load(file = here(paste0('data/', fl, '.RData')))

toplo <- get(fl)

p <- show_prdseries(toplo, ylab = 'Bottom temperature (\u00B0 C)') + 
  labs(subtitle ='Station 66, OTB')

png(here('figs/tempgamex.png'), height = 3, width = 6, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# temp thresholds doy -------------------------------------------------------------------------

# model file with daily temp predictions for each station, bottom/surface temp
load(file = here('data/tempprd.RData'))

tempsum <- tempprd %>%
  mutate(
    cnts = purrr::map(prd, function(x){
      
      tibble(thr = tempthr) %>%
        group_nest(thr) %>%
        mutate(data = list(x)) %>%
        mutate(
          data = purrr::pmap(list(data, thr), function(data, thr){
            
            data %>%
              filter(yr < 2023) %>%
              summarise(
                cnt = sum(value > thr),
                frt = date[which(value >= thr)[1]],
                lst = date[rev(which(value >= thr))[1]],
                .by = 'yr'
              )
            
          })
        ) %>% 
        unnest('data')
      
    })
  )

toplo1 <- tempsum %>% 
  filter(param %in% 'tempbot') %>% 
  select(-AIC, -GCV, -R2, -prd) %>% 
  unnest('cnts') %>% 
  summarise(
    frt = median(frt, na.rm = T), 
    lst = median(lst, na.rm = T),
    cnt = mean(cnt),
    .by = c('bay_segment', 'thr', 'yr', 'param')
  ) %>% 
  mutate(
    frt = yday(frt), 
    frt = as.Date(date_decimal(2020 + frt / 365)),
    lst = yday(lst), 
    lst = as.Date(date_decimal(2020 + lst / 365)),
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')), 
    thr = paste0(thr, '\u00B0 C')
  ) %>% 
  pivot_longer(frt:lst, names_to = 'dts', values_to = 'doy') 

toplo2 <- toplo1 %>% 
  group_nest(bay_segment, param, thr, dts) %>% 
  mutate(
    reg = purrr::map(data, function(x){
      
      dat <- x[x$doy > 0, ]
      mod <- lm(doy ~ yr, dat)
      prddat <- tibble(yr = range(dat$yr, na.rm = T)) %>% 
        mutate(
          doy = predict(mod, newdata = .), 
          doy = as.Date(doy, origin = '1969-12-31')
        )
      
      return(prddat)
      
    })
  ) %>% 
  select(-data) %>% 
  unnest('reg') %>% 
  mutate(
    ind = case_when(
      dts == 'frt' & yr == max(yr) ~ 1, 
      dts == 'frt' & yr == min(yr) ~ 2,
      dts == 'lst' & yr == min(yr) ~ 3,
      dts == 'lst' & yr == max(yr) ~ 4
    ),
    .by = c(bay_segment, param, thr)
  ) %>% 
  arrange(bay_segment, param, thr, ind)

p <- ggplot(toplo1, aes(x = yr, y = doy)) +
  geom_polygon(data = toplo2, aes(x = yr, y = doy, fill = thr), alpha = 0.5) + 
  geom_point(size = 0.5, aes(color = thr)) +
  geom_smooth(aes(group = dts, color = thr), method = 'lm', se = F, lwd = 1) +
  facet_grid(bay_segment ~ thr) +
  scale_colour_manual(values = tempcol) +
  scale_fill_manual(values = tempcol) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_date(date_breaks = '2 months', date_labels = paste('%b', '1st')) +
  coord_cartesian(
    ylim = as.Date(c('2020-05-01', '2020-11-01'))
  ) +
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(), 
    strip.background = element_blank(), 
    strip.text = element_text(size = 12) ,
    axis.text = element_text(size = 9), 
    legend.position = 'none'
  ) +
  labs(
    x = NULL, 
    y = "Day of year", 
    title = "Seasonal duration of time above temperature thresholds", 
    subtitle = "First and last date exceeded by year"
  )

png(here('figs/tempexceed.png'), height = 7, width = 6, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# temp thresholds length of time --------------------------------------------------------------

# model file with daily temp predictions for each station, bottom/surface temp
load(file = here('data/tempprd.RData'))

tempsum <- tempprd %>%
  mutate(
    cnts = purrr::map(prd, function(x){
      
      tibble(thr = tempthr) %>%
        group_nest(thr) %>%
        mutate(data = list(x)) %>%
        mutate(
          data = purrr::pmap(list(data, thr), function(data, thr){
            
            data %>%
              filter(yr < 2023) %>%
              summarise(
                cnt = sum(value > thr),
                frt = date[which(value >= thr)[1]],
                lst = date[rev(which(value >= thr))[1]],
                .by = 'yr'
              )
            
          })
        ) %>% 
        unnest('data')
      
    })
  )

toplo <- tempsum %>% 
  select(-AIC, -GCV, -R2, -prd) %>% 
  unnest('cnts') %>% 
  summarise(
    frt = median(frt, na.rm = T), 
    lst = median(lst, na.rm = T),
    cnt = median(cnt),
    .by = c('bay_segment', 'thr', 'yr', 'param')
  ) %>% 
  mutate(
    frt = yday(frt), 
    lst = yday(lst), 
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')), 
    thr = paste0(thr, '\u00B0'),
    param = factor(param, levels = c('temptop', 'tempbot'), labels = c('Top', 'Bottom'))
    # cnt = lst - frt
  ) 

p <- ggplot(toplo, aes(x = yr, y = cnt, group = thr, color = thr)) + 
  geom_point(size = 0.5) + 
  geom_smooth(formula = y ~ x, method = 'lm', se = F) + 
  facet_grid(param ~ bay_segment) + 
  scale_colour_manual(values = tempcol) +
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 9), 
    axis.text.x = element_text(size = 8),
    legend.position = 'top'
  ) +
  labs(
    x = NULL,
    y = "Number of days", 
    title = "Number of days above threshold over time", 
    subtitle = "Results for top and bottom water temperatures by bay segment", 
    color = "Threshold (C)"
  )

png(here('figs/tempcount.png'), height = 7, width = 7, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# sali station gam --------------------------------------------------------------------------

fl <- 'OTB_66_salibot'
load(file = here(paste0('data/', fl, '.RData')))

toplo <- get(fl)

p <- show_prdseries(toplo, ylab = 'Bottom salinity (psu)') + 
  labs(subtitle ='Station 66, OTB')

png(here('figs/saligamex.png'), height = 3, width = 6, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# sali thresholds doy -------------------------------------------------------------------------

# model file with daily sali predictions for each station, bottom/surface sali
load(file = here('data/saliprd.RData'))

salisum <- saliprd %>%
  mutate(
    cnts = purrr::map(prd, function(x){
      
      tibble(thr = salithr) %>%
        group_nest(thr) %>%
        mutate(data = list(x)) %>%
        mutate(
          data = purrr::pmap(list(data, thr), function(data, thr){
            
            data %>%
              filter(yr < 2023) %>%
              summarise(
                cnt = sum(value < thr),
                frt = date[which(value <= thr)[1]],
                lst = date[rev(which(value <= thr))[1]],
                .by = 'yr'
              )
            
          })
        ) %>% 
        unnest('data')
      
    })
  )

toplo1 <- salisum %>% 
  filter(param %in% 'salibot') %>% 
  select(-AIC, -GCV, -R2, -prd) %>% 
  unnest('cnts') %>% 
  summarise(
    frt = median(frt, na.rm = T), 
    lst = median(lst, na.rm = T),
    cnt = mean(cnt),
    .by = c('bay_segment', 'thr', 'yr', 'param')
  ) %>% 
  mutate(
    frt = yday(frt), 
    frt = as.Date(date_decimal(2020 + frt / 365)),
    lst = yday(lst), 
    lst = as.Date(date_decimal(2020 + lst / 365)),
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')), 
    thr = paste0(thr, 'psu')
  ) %>% 
  pivot_longer(frt:lst, names_to = 'dts', values_to = 'doy') 

toplo2 <- toplo1 %>% 
  group_nest(bay_segment, param, thr, dts) %>% 
  mutate(
    reg = purrr::map(data, function(x){
      
      dat <- x[x$doy > 0, ]
      mod <- lm(doy ~ yr, dat)
      prddat <- tibble(yr = range(dat$yr, na.rm = T)) %>% 
        mutate(
          doy = predict(mod, newdata = .), 
          doy = as.Date(doy, origin = '1969-12-31')
        )
      
      return(prddat)
      
    })
  ) %>% 
  select(-data) %>% 
  unnest('reg') %>% 
  mutate(
    ind = case_when(
      dts == 'frt' & yr == max(yr) ~ 1, 
      dts == 'frt' & yr == min(yr) ~ 2,
      dts == 'lst' & yr == min(yr) ~ 3,
      dts == 'lst' & yr == max(yr) ~ 4
    ),
    .by = c(bay_segment, param, thr)
  ) %>% 
  arrange(bay_segment, param, thr, ind)

p <- ggplot(toplo1, aes(x = yr, y = doy)) +
  geom_polygon(data = toplo2, aes(x = yr, y = doy, fill = thr), alpha = 0.5) +
  geom_point(size = 0.5, aes(color = thr)) +
  geom_smooth(aes(group = dts, color = thr), method = 'lm', se = F, lwd = 1) +
  facet_grid(bay_segment ~ thr) +
  scale_colour_manual(values = salicol) +
  scale_fill_manual(values = salicol) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_date(date_breaks = '2 months', date_labels = paste('%b', '1st')) +
  # coord_cartesian(
  #   ylim = as.Date(c('2020-05-01', '2020-11-01'))
  # ) +
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(), 
    strip.background = element_blank(), 
    strip.text = element_text(size = 12) ,
    axis.text = element_text(size = 9), 
    legend.position = 'none'
  ) +
  labs(
    x = NULL, 
    y = "Day of year", 
    title = "Seasonal duration of time below salinity thresholds", 
    subtitle = "First and last date exceeded by year"
  )

png(here('figs/saliexceed.png'), height = 7, width = 6, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# sali thresholds length of time --------------------------------------------------------------

# model file with daily sali predictions for each station, bottom/surface temp
load(file = here('data/saliprd.RData'))

salisum <- saliprd %>%
  mutate(
    cnts = purrr::map(prd, function(x){
      
      tibble(thr = salithr) %>%
        group_nest(thr) %>%
        mutate(data = list(x)) %>%
        mutate(
          data = purrr::pmap(list(data, thr), function(data, thr){
            
            data %>%
              filter(yr < 2023) %>%
              summarise(
                cnt = sum(value <= thr),
                frt = date[which(value <= thr)[1]],
                lst = date[rev(which(value <= thr))[1]],
                .by = 'yr'
              )
            
          })
        ) %>% 
        unnest('data')
      
    })
  )

toplo <- salisum %>% 
  select(-AIC, -GCV, -R2, -prd) %>% 
  unnest('cnts') %>% 
  summarise(
    frt = median(frt, na.rm = T), 
    lst = median(lst, na.rm = T),
    cnt = median(cnt),
    .by = c('bay_segment', 'thr', 'yr', 'param')
  ) %>% 
  mutate(
    frt = yday(frt), 
    lst = yday(lst), 
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')), 
    thr = paste0(thr, 'psu'),
    param = factor(param, levels = c('salitop', 'salibot'), labels = c('Top', 'Bottom'))
    # cnt = lst - frt
  ) 

p <- ggplot(toplo, aes(x = yr, y = cnt, group = thr, color = thr)) + 
  geom_point(size = 0.5) + 
  geom_smooth(formula = y ~ x, method = 'lm', se = F) + 
  facet_grid(param ~ bay_segment) + 
  scale_colour_manual(values = salicol) +
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 9), 
    axis.text.x = element_text(size = 8),
    legend.position = 'top'
  ) +
  labs(
    x = NULL,
    y = "Number of days", 
    title = "Number of days below threshold over time", 
    subtitle = "Results for top and bottom water salinity by bay segment", 
    color = "Threshold (psu)"
  )

png(here('figs/salicount.png'), height = 7, width = 7, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()