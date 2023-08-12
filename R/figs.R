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
        pval = ests$p.value[2], 
        slos = ests$estimate[2],
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


# station gam -------------------------------------------------------------------------------

fl <- 'OTB_66_temptop'
load(file = here(paste0('data/', fl, '.RData')))

toplo <- get(fl)

p <- show_prdseries(toplo, ylab = 'Surface tempeature (\u00B0 C') + 
  labs(subtitle ='Station 66, OTB')

png(here('figs/gamex.png'), height = 3, width = 6, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# temp thresholds doy -------------------------------------------------------------------------

# model file with daily temp predictions for each station, bottom/surface temp
load(file = here('data/moddat.RData'))

modsum <- moddat %>%
  mutate(
    cnts = purrr::map(prd, function(x){
      
      tibble(thr = c(20, 25, 30)) %>%
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

toplo1 <- modsum %>% 
  filter(param %in% 'temptop') %>% 
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
  arrange(param, thr, doy, yr)

p <- ggplot(toplo1, aes(x = yr, y = doy)) +
  geom_polygon(data = toplo2, aes(x = yr, y = doy, fill = thr), alpha = 0.5) + 
  geom_point(size = 0.5, aes(color = thr)) +
  geom_smooth(aes(group = dts, color = thr), method = 'lm', se = F, lwd = 1) +
  facet_grid(bay_segment ~ thr) +
  scale_colour_manual(values = c('coral', 'red2', 'darkred')) +
  scale_fill_manual(values = c('coral', 'red2', 'darkred')) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_date(date_breaks = '3 months', date_labels = paste('%b', '1st')) +
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

png(here('figs/threxceed.png'), height = 7, width = 6, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# temp thresholds length of time --------------------------------------------------------------

# model file with daily temp predictions for each station, bottom/surface temp
load(file = here('data/moddat.RData'))

modsum <- moddat %>%
  mutate(
    cnts = purrr::map(prd, function(x){
      
      tibble(thr = c(20, 25, 30)) %>%
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

toplo <- modsum %>% 
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
  scale_colour_manual(values = c('coral', 'red2', 'darkred')) +
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


png(here('figs/thrcount.png'), height = 7, width = 7, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()
