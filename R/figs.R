library(tidyverse)
library(lubridate)
library(EnvStats)
library(tbeptools)
library(ggmap)
library(sf)
library(patchwork)
library(here)

# kendall all months --------------------------------------------------------------------------

sktres <- epcdata %>% 
  select(bay_segment, station = epchc_station, SampleTime, lon = Longitude, lat = Latitude, yr, 
         mo, matches('Temp')) %>% 
  pivot_longer(matches('Temp'), names_to = 'var', values_to = 'val') %>% 
  nest(.by = c('bay_segment', 'var', 'station', 'lon', 'lat')) %>% 
  crossing(
    ., 
    tibble(
      per = c('all', 'rec', 'cur'),
      yrs = c('1974-2023', '2000-2016', '2016-2023')
    )
  ) %>% 
  mutate(
    skt = purrr::pmap(list(station, var, yrs, data), function(station, var, yrs, data){
      
      cat(station, var, yrs, '\n')
      
      yrsel <- strsplit(yrs, '-') %>% 
        .[[1]] %>% 
        as.numeric()
      
      # yr selection
      tomod <- data %>% 
        arrange(yr, mo) %>% 
        na.omit() %>% 
        filter(yr >= yrsel[1] & yr < yrsel[2])
      
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

leglab <- "\u00B0C / year"

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
    var = factor(var, 
      levels = c("Temp_Water_Top_degC", "Temp_Water_Mid_degC", "Temp_Water_Bottom_degC"), 
      labels = c('Top', 'Mid', 'Bottom')
    ),
    yrs = factor(yrs, 
      levels = c('1974-2023', '2000-2016', '2016-2023'), 
      labels = c('1974-2022', '2000-2015', '2016-2022')
    )
  ) 

p <- ggmap(bsmap1_transparent) +
  geom_point(data = toplo, aes(x = lon, y = lat, size = abs(slos), fill = slos, shape = coefsgn, color = pvalcol), stroke = 1) +
  facet_grid(var ~ yrs) +
  # scale_fill_gradient2(leglab, low = 'green', mid = 'grey',  high = 'tomato1', midpoint = 0) +
  scale_fill_gradientn(leglab, limits = c(-0.25, 0.25), colors = c('green', 'grey', 'tomato1')) +
  scale_color_manual(values = c(scales::alpha('black', 0), 'black'), guide = 'none', drop = FALSE) +
  coord_map() +
  scale_shape_manual(values = c(24, 25), guide = 'none', drop = FALSE) +
  pthm +
  scale_size(range = c(3, 8), guide = F) +
  labs(
    caption = 'Outlines indicate p < 0.05'
  )

png(here('figs/sktall.png'), height = 10, width = 9, family = 'serif', units = 'in', res = 300)
p
dev.off()

# kendall by month ----------------------------------------------------------------------------

sktres <- epcdata %>% 
  select(bay_segment, station = epchc_station, SampleTime, lon = Longitude, lat = Latitude, yr, 
         mo, matches('Temp')) %>% 
  pivot_longer(matches('Temp'), names_to = 'var', values_to = 'val') %>% 
  nest(.by = c('bay_segment', 'var', 'station', 'lon', 'lat', 'mo')) %>% 
  crossing(
    ., 
    tibble(
      per = c('all', 'rec', 'cur'),
      yrs = c('1974-2023', '2000-2016', '2016-2023')
    )
  ) %>% 
  mutate(
    skt = purrr::pmap(list(station, var, yrs, mo, data), function(station, var, yrs, mo, data){
      
      cat(station, var, yrs, mo, '\n')
      
      yrsel <- strsplit(yrs, '-') %>% 
        .[[1]] %>% 
        as.numeric()
      
      # yr selection
      tomod <- data %>% 
        arrange(yr) %>% 
        na.omit() %>% 
        filter(yr >= yrsel[1] & yr < yrsel[2])
      
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

leglab <- "\u00B0C / year"

pthm <- theme_bw(base_family = 'serif', base_size = 14) +
  theme(
    legend.position = 'right',
    # legend.box = 'vertical',
    strip.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 14)
  )

sktres %>%
  filter(slos < 1.5) %>% 
  mutate(
    pvalcol = ifelse(pval < 0.05, T, F),
    coefsgn = sign(slos),
    coefsgn = factor(coefsgn, levels = c('1', '-1'), labels = c('inc', 'dec')), 
    var = factor(var, 
                 levels = c("Temp_Water_Top_degC", "Temp_Water_Mid_degC", "Temp_Water_Bottom_degC"), 
                 labels = c('Top', 'Mid', 'Bottom')
    ),
    mo = sprintf('%02d', mo), 
    sz = scales::rescale(abs(slos), c(3, 8))
  ) %>% 
  group_nest(mo) %>% 
  deframe() %>% 
  iwalk(function(x, mo){

    p <- ggmap(bsmap1_transparent) +
      geom_point(data = x, aes(x = lon, y = lat, fill = slos, shape = coefsgn, color = pvalcol), stroke = 1, 
                 size = x$sz) +
      facet_grid(var ~ yrs) +
      scale_fill_gradientn(leglab, limits = c(-1.5, 1.5), colors = c('green', 'grey', 'tomato1')) +
      # scale_fill_gradient2(leglab, low = 'green', mid = 'grey',  high = 'tomato1', midpoint = 0) +
      scale_color_manual(values = c(scales::alpha('black', 0), 'black'), guide = 'none', drop = FALSE) +
      coord_map() +
      scale_shape_manual(values = c(24, 25), guide = 'none', drop = FALSE) +
      pthm +
      # scale_size(range = c(1, 6), guide = F) +
      labs(
        caption = 'Outlines indicate p < 0.05'
      )
    
    fl <- paste0('figs/skt_',sprintf('%02d', as.numeric(mo)), '.png')
    png(here(fl), height = 11, width = 10, family = 'serif', units = 'in', res = 300)
    print(p)
    dev.off()
    
  })
