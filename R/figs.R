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
library(ggspatial)
library(modelbased)
library(lmerTest)
library(vegan)
library(ggord)
library(scales)
library(visreg)
library(car)
library(mgcv)

source(here('R/funcs.R'))

seglng <- c('Old Tampa Bay', 'Hillsborough Bay', 'Middle Tampa Bay', 'Lower Tampa Bay')
segshr <- c('OTB', 'HB', 'MTB', 'LTB')

yrsel1 <- c(1975, 2022)
yrsel2 <- c(1998, 2022)
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

savval <- c('total', 'Halodule', 'Syringodium', 'Thalassia')
savleg <- expression('total', italic('H. wrightii'), italic('S. filiforme'), italic('T. testudinum'))
savcol <- c('darkgrey', '#ED90A4', '#CCA65A', '#7EBA68')
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

toplo2 <- transectavespp %>% 
  filter(bay_segment %in% c('OTB', 'HB', 'MTB', 'LTB')) %>% 
  filter(Savspecies %in% c('total', 'Halodule', 'Syringodium', 'Thalassia'))

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
  coord_cartesian(ylim = c(0, 1)) + 
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

# meteorological, salinity, water temp obs data -----------------------------------------------

load(file = here('data/speidat.RData'))
load(file = url("https://github.com/tbep-tech/load-estimates/raw/main/data/totanndat.RData"))

# hydro load, otb, mtb, ltb, hb
hydrodat <- totanndat %>% 
  select(year, bay_segment, hy_load) %>% 
  filter(bay_segment %in% c('Hillsborough Bay', 'Old Tampa Bay', 'Lower Tampa Bay', 'Middle Tampa Bay')) %>% 
  summarise(
    hy_load = sum(hy_load) / 1e3, 
    .by = 'year'
  ) 

toplo1 <- speidat %>% 
  select(yr, precip_mm, tavg_c) %>% 
  summarise(
    precip_mm = sum(precip_mm), 
    tavg_c = mean(tavg_c), 
    .by = 'yr'
  )

thm <- theme_minimal() + 
  theme(
    strip.placement = 'outside', 
    panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 12), 
    legend.text = element_text(size= 12)
  )

p1 <- ggplot(toplo1, aes(x = yr, y = tavg_c)) + 
  geom_line() + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F, formula = y~x, color = 'red2') +
  thm +
  theme(
    axis.text.x = element_blank()
  ) +
  labs(
    x = NULL, 
    y = 'Air temp. (\u00B0 C)'
  )

p2 <- ggplot(toplo1, aes(x = yr, y = precip_mm / 1e3)) + 
  geom_line() + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F, formula = y~x, color = 'blue') +
  thm +
  theme(
    axis.text.x = element_blank()
  ) +
  labs(
    x = NULL, 
    y = 'Precip. (m)'
  )

p3 <- ggplot(hydrodat, aes(x = year, y = hy_load)) +
  geom_line() + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F, formula = y~x, color = 'blue') +
  coord_cartesian(xlim = range(speidat$yr)) + 
  thm +
  theme(
    axis.text.x = element_blank()
  ) +
  labs(
    x = NULL, 
    y = expression(paste('Hyd. load (', 10^3, ' t/yr)'))
  )

p4 <- ggplot(speidat, aes(x = date, y = spi, fill = spisign)) + 
  geom_col(width = 35) +
  scale_fill_manual(values = c('grey70', 'grey30'), guide = 'none') + 
  thm  +
  theme(
    axis.text.x = element_text(size = 12)
  ) +
  coord_cartesian(xlim = range(speidat$date)) +
  labs(
    x = NULL, 
    y = 'SPI (z-values)'
  )

toplo <- epcdata %>% 
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
    var = factor(var, levels = c('Temp', 'Sal'), labels = c('Water temp. (\u00B0C)', 'Salinity (ppt)')),
    loc = factor(loc, levels = c('Top', 'Bottom'))
  ) %>% 
  filter(!is.na(val)) %>% 
  summarise(
    avev = mean(val, na.rm = T),
    lov = t.test(val, na.rm = T)$conf.int[1],
    hiv = t.test(val, na.rm = T)$conf.int[2],
    .by = c('bay_segment', 'yr', 'var', 'loc') 
  )

wd <- 0.5

toplo5 <- toplo %>% 
  filter(var == 'Water temp. (\u00B0C)') %>% 
  filter(!(yr %in% c(1982, 1985) & bay_segment == 'OTB')) %>%  # missing months create outliers
  filter(!(yr %in% 1975 & bay_segment == 'HB')) # missing months create outliers
p5 <- ggplot(toplo5, aes(x = yr, y = avev, group = loc, color = loc)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), position = position_dodge2(width = wd), show.legend = F, alpha = 0.7) + 
  geom_point(position = position_dodge2(width = wd), size = 0.5) +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  geom_smooth(method = 'lm', formula = y ~ x, se = F) +
  facet_grid(~ bay_segment) +
  scale_color_manual(values = c( 'steelblue4', 'steelblue1')) +
  thm +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_blank()
  ) +
  labs(
    x = NULL, 
    y = 'Water temp. (\u00B0C)', 
    color = NULL, 
    shape = NULL
  )
toplo6 <- toplo %>% 
  filter(var == 'Salinity (ppt)')
p6 <- ggplot(toplo6, aes(x = yr, y = avev, group = loc, color = loc)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), position = position_dodge2(width = wd), show.legend = F, alpha = 0.7) + 
  geom_point(position = position_dodge2(width = wd), size = 0.5) +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  geom_smooth(method = 'lm', formula = y ~ x, se = F) +
  facet_grid(~ bay_segment) +
  scale_color_manual(values = c('steelblue4', 'steelblue1')) +
  thm +
  theme(
    strip.text = element_blank()
  ) +
  labs(
    x = NULL, 
    y = 'Salinity (ppt)', 
    color = NULL, 
    shape = NULL
  )

p <- p1 + p2 + p3 + p4 + (p5 + p6 + plot_layout(ncol = 1, guides = 'collect')) + plot_layout(ncol = 1, heights = c(1, 1, 1, 1, 2.5)) & theme(legend.position = 'top')

png(here('figs/meteowqraw.png'), height = 8.5, width = 7, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()

# supp meteorological, salinity, water temp obs data ------------------------------------------

load(file = here('data/speidat.RData'))
load(file = url("https://github.com/tbep-tech/load-estimates/raw/main/data/totanndat.RData"))

# hydro load, otb, mtb, ltb, hb
hydrodat <- totanndat %>% 
  select(year, bay_segment, hy_load) %>% 
  filter(bay_segment %in% c('Hillsborough Bay', 'Old Tampa Bay', 'Lower Tampa Bay', 'Middle Tampa Bay')) %>% 
  summarise(
    hy_load = sum(hy_load) / 1e3, 
    .by = 'year'
  ) %>% 
  filter(year >= yrsel2[1] & year <= yrsel2[2])

speidat <- speidat%>% 
  filter(yr >= yrsel2[1] & yr <= yrsel2[2])

toplo1 <- speidat %>% 
  select(yr, precip_mm, tavg_c) %>% 
  summarise(
    precip_mm = sum(precip_mm), 
    tavg_c = mean(tavg_c), 
    .by = 'yr'
  )

thm <- theme_minimal() + 
  theme(
    strip.placement = 'outside', 
    panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 12), 
    legend.text = element_text(size= 12)
  )

p1 <- ggplot(toplo1, aes(x = yr, y = tavg_c)) + 
  geom_line() + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F, formula = y~x, color = 'red2') +
  thm +
  theme(
    axis.text.x = element_blank()
  ) +
  labs(
    x = NULL, 
    y = 'Air temp. (\u00B0 C)'
  )

p2 <- ggplot(toplo1, aes(x = yr, y = precip_mm / 1e3)) + 
  geom_line() + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F, formula = y~x, color = 'blue') +
  thm +
  theme(
    axis.text.x = element_blank()
  ) +
  labs(
    x = NULL, 
    y = 'Precip. (m)'
  )

p3 <- ggplot(hydrodat, aes(x = year, y = hy_load)) +
  geom_line() + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F, formula = y~x, color = 'blue') +
  coord_cartesian(xlim = range(toplo1$yr)) + 
  thm +
  theme(
    axis.text.x = element_blank()
  ) +
  labs(
    x = NULL, 
    y = expression(paste('Hyd. load (', 10^3, ' t/yr)'))
  )

p4 <- ggplot(speidat, aes(x = date, y = spi, fill = spisign)) + 
  geom_col(width = 35) +
  scale_fill_manual(values = c('grey70', 'grey30'), guide = 'none') + 
  thm +
  theme(
    axis.text.x = element_text(size = 12)
  ) +
  coord_cartesian(xlim = range(speidat$date)) +
  labs(
    x = NULL, 
    y = 'SPI (z-values)'
  )

toplo <- epcdata %>% 
  select(bay_segment, epchc_station, SampleTime, yr, matches('Top|Bottom')) %>% 
  filter(yr >= yrsel2[1] & yr <= yrsel2[2]) %>% 
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
    var = factor(var, levels = c('Temp', 'Sal'), labels = c('Water temp. (\u00B0C)', 'Salinity (ppt)')),
    loc = factor(loc, levels = c('Top', 'Bottom'))
  ) %>% 
  filter(!is.na(val)) %>% 
  summarise(
    avev = mean(val, na.rm = T),
    lov = t.test(val, na.rm = T)$conf.int[1],
    hiv = t.test(val, na.rm = T)$conf.int[2],
    .by = c('bay_segment', 'yr', 'var', 'loc') 
  )

wd <- 0.5

toplo5 <- toplo %>% 
  filter(var == 'Water temp. (\u00B0C)') %>% 
  filter(!(yr %in% c(1982, 1985) & bay_segment == 'OTB')) %>%  # missing months create outliers
  filter(!(yr %in% 1975 & bay_segment == 'HB')) # missing months create outliers
p5 <- ggplot(toplo5, aes(x = yr, y = avev, group = loc, color = loc)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), position = position_dodge2(width = wd), show.legend = F, alpha = 0.7) + 
  geom_point(position = position_dodge2(width = wd), size = 0.5) +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  geom_smooth(method = 'lm', formula = y ~ x, se = F) +
  facet_grid(~ bay_segment) +
  scale_color_manual(values = c( 'steelblue4', 'steelblue1')) +
  thm +
  theme(
    axis.text.x = element_blank(),
    strip.text = element_text(size = 12)
    ) +
  labs(
    x = NULL, 
    y = 'Water temp. (\u00B0C)', 
    color = NULL, 
    shape = NULL
  )
toplo6 <- toplo %>% 
  filter(var == 'Salinity (ppt)')
p6 <- ggplot(toplo6, aes(x = yr, y = avev, group = loc, color = loc)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), position = position_dodge2(width = wd), show.legend = F, alpha = 0.7) + 
  geom_point(position = position_dodge2(width = wd), size = 0.5) +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  geom_smooth(method = 'lm', formula = y ~ x, se = F) +
  facet_grid(~ bay_segment) +
  scale_color_manual(values = c('steelblue4', 'steelblue1')) +
  thm +
  theme(
    strip.text = element_blank()
  ) +
  labs(
    x = NULL, 
    y = 'Salinity (ppt)', 
    color = NULL, 
    shape = NULL
  )

p <- p1 + p2 + p3 + p4 + (p5 + p6 + plot_layout(ncol = 1, guides = 'collect')) + plot_layout(ncol = 1, heights = c(1, 1, 1, 1, 2.5)) & theme(legend.position = 'top')

png(here('figs/suppmeteowqraw.png'), height = 8.5, width = 7, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()

# kendall -------------------------------------------------------------------------------------

leglab <- expression(paste(yr^{-1}))

# kendall all years
sktres <- epcdata %>% 
  select(bay_segment, station = epchc_station, SampleTime, lon = Longitude, lat = Latitude, yr, 
         mo, matches('Top|Bottom')) %>% 
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

#kendall by month
ktres <- epcdata %>% 
  select(bay_segment, station = epchc_station, SampleTime, lon = Longitude, lat = Latitude, yr, 
         mo, matches('Top|Bottom')) %>% 
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

# plots of segment summaries by month
ktresplo <- ktres %>% 
  mutate(
    mo = month(mo, label = T), 
    bay_segment = factor(bay_segment, levels = c('LTB', 'MTB', 'HB', 'OTB')), 
    loc =  gsub("^.*_(.*)_.*$", "\\1", var),
    varsimp = gsub('^(.*?)_.*$', '\\1', var, perl = T),
    varsimp = factor(varsimp, levels = c('Temp', 'Sal'), labels = c('Temperature', 'Salinity')),
    loc = factor(loc, levels = c('Top', 'Bottom'))
  ) %>%
  nest(.by = c('bay_segment', 'mo', 'varsimp', 'var', 'loc')) %>% 
  mutate(
    sum = purrr::pmap(list(varsimp, data), function(varsimp, data){
      
      data %>% 
        summarise(
          nsig = case_when(
            varsimp == 'Salinity' ~ -1 * sum(pval < 0.05 & slos < 0),
            varsimp == 'Temperature' ~ sum(pval < 0.05 & slos > 0)
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
  mutate(
    plo = purrr::pmap(list(loc, data), function(loc, data){

      p <- ggplot(data, aes(x = mo, y = bay_segment, fill = nsigper)) + 
        geom_tile(color = 'darkgrey') +
        geom_text(aes(label = perlab), color = 'white', fontface = 'bold') +
        scale_x_discrete(expand = c(0, 0)) + 
        scale_y_discrete(expand = c(0, 0)) + 
        scale_fill_gradientn(leglab, limits = c(-100, 100), colors = c('blue', 'grey', 'tomato1')) +
        facet_wrap(~ varsimp, ncol = 1, scales = 'free_x') + 
        theme(
          strip.background = element_blank(), 
          strip.text = element_text(hjust = 0, size = 12), 
          legend.position = 'none',
          axis.text = element_text(size = 10)
        ) + 
        labs(
          x = NULL, 
          y = 'Bay segment'
        )
      
      return(p)

    })
  )

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

toplo <- sktres %>%
  mutate(
    pvalcol = ifelse(pval < 0.05, T, F),
    coefsgn = sign(slos),
    coefsgn = factor(coefsgn, levels = c('1', '-1'), labels = c('inc', 'dec')), 
    loc =  gsub("^.*_(.*)_.*$", "\\1", var),
    var = gsub('^(.*?)_.*$', '\\1', var, perl = T),
    var = factor(var, levels = c('Temp', 'Sal'), labels = c('Temperature', 'Salinity')),
    loc = factor(loc, levels = c('Top', 'Bottom'))
  )

pthm <- theme_bw(base_family = 'serif') +
  theme(
    legend.position = c(0.9, 0.1), 
    # legend.box = 'horizontal',
    strip.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    legend.text = element_text(size = 8),
    legend.background = element_rect(fill = NA)
  )

p1 <- ggmap(bsmap1_transparent) +
  geom_sf(data = tbseg, inherit.aes = F) +
  geom_point(data = toplo, aes(x = lon, y = lat, size = abs(slos), fill = slos, shape = coefsgn, color = pvalcol), stroke = 1) +
  facet_grid(loc ~ var) +
  # scale_fill_gradient2(leglab, low = 'blue', mid = 'grey',  high = 'tomato1', midpoint = 0) +
  scale_fill_gradientn(leglab, limits = colrng, colors = c('blue', 'grey', 'tomato1')) +
  scale_color_manual(values = c(scales::alpha('black', 0), 'black'), guide = 'none', drop = FALSE) +
  # coord_map() +
  scale_shape_manual(values = c(24, 25), drop = FALSE, guide = 'none') +
  pthm +
  scale_size(range = c(0.75, 4), guide = 'none') +
  guides(fill = guide_colourbar(barwidth = 0.4, barheight = 2.5)) + 
  labs(
    title = paste0('(a) Change per year, ', yrsel1[1], '-', yrsel1[2]),
    caption = 'Outlines indicate p < 0.05'
  )

p2 <- ktresplo %>% 
  filter(loc == 'Top') %>% 
  pull(plo) %>% 
  .[[1]] + labs(title = '(b) Top, % stations with significant trends by month')

p3 <- ktresplo %>% 
  filter(loc == 'Bottom') %>% 
  pull(plo) %>% 
  .[[1]] + labs(title = '(c) Bottom, % stations with significant trends by month')

p <- p1 + (p2 + p3 + plot_layout(ncol = 1)) + plot_layout(ncol = 2, width = c(0.9, 1))

png(here('figs/kendall.png'), height = 7, width = 9.5, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# supp kendall --------------------------------------------------------------------------------

leglab <- expression(paste(yr^{-1}))

# kendall all years
sktres <- epcdata %>% 
  select(bay_segment, station = epchc_station, SampleTime, lon = Longitude, lat = Latitude, yr, 
         mo, matches('Top|Bottom')) %>% 
  filter(yr >= yrsel2[1] & yr <= yrsel2[2]) %>% 
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

#kendall by month
ktres <- epcdata %>% 
  select(bay_segment, station = epchc_station, SampleTime, lon = Longitude, lat = Latitude, yr, 
         mo, matches('Top|Bottom')) %>% 
  filter(yr >= yrsel2[1] & yr <= yrsel2[2]) %>% 
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

# plots of segment summaries by month
ktresplo <- ktres %>% 
  mutate(
    mo = month(mo, label = T), 
    bay_segment = factor(bay_segment, levels = c('LTB', 'MTB', 'HB', 'OTB')), 
    loc =  gsub("^.*_(.*)_.*$", "\\1", var),
    varsimp = gsub('^(.*?)_.*$', '\\1', var, perl = T),
    varsimp = factor(varsimp, levels = c('Temp', 'Sal'), labels = c('Temperature', 'Salinity')),
    loc = factor(loc, levels = c('Top', 'Bottom'))
  ) %>%
  nest(.by = c('bay_segment', 'mo', 'varsimp', 'var', 'loc')) %>% 
  mutate(
    sum = purrr::pmap(list(varsimp, data), function(varsimp, data){
      
      data %>% 
        summarise(
          nsig = case_when(
            varsimp == 'Salinity' ~ -1 * sum(pval < 0.05 & slos < 0),
            varsimp == 'Temperature' ~ sum(pval < 0.05 & slos > 0)
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
  mutate(
    plo = purrr::pmap(list(loc, data), function(loc, data){
      
      p <- ggplot(data, aes(x = mo, y = bay_segment, fill = nsigper)) + 
        geom_tile(color = 'darkgrey') +
        geom_text(aes(label = perlab), color = 'white', fontface = 'bold') +
        scale_x_discrete(expand = c(0, 0)) + 
        scale_y_discrete(expand = c(0, 0)) + 
        scale_fill_gradientn(leglab, limits = c(-100, 100), colors = c('blue', 'grey', 'tomato1')) +
        facet_wrap(~ varsimp, ncol = 1, scales = 'free_x') + 
        theme(
          strip.background = element_blank(), 
          strip.text = element_text(hjust = 0, size = 12), 
          legend.position = 'none',
          axis.text = element_text(size = 10)
        ) + 
        labs(
          x = NULL, 
          y = 'Bay segment'
        )
      
      return(p)
      
    })
  )

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

toplo <- sktres %>%
  mutate(
    pvalcol = ifelse(pval < 0.05, T, F),
    coefsgn = sign(slos),
    coefsgn = factor(coefsgn, levels = c('1', '-1'), labels = c('inc', 'dec')), 
    loc =  gsub("^.*_(.*)_.*$", "\\1", var),
    var = gsub('^(.*?)_.*$', '\\1', var, perl = T),
    var = factor(var, levels = c('Temp', 'Sal'), labels = c('Temperature', 'Salinity')),
    loc = factor(loc, levels = c('Top', 'Bottom'))
  )

pthm <- theme_bw(base_family = 'serif') +
  theme(
    legend.position = c(0.9, 0.1), 
    # legend.box = 'horizontal',
    strip.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    legend.text = element_text(size = 8),
    legend.background = element_rect(fill = NA)
  )

p1 <- ggmap(bsmap1_transparent) +
  geom_sf(data = tbseg, inherit.aes = F) +
  geom_point(data = toplo, aes(x = lon, y = lat, size = abs(slos), fill = slos, shape = coefsgn, color = pvalcol), stroke = 1) +
  facet_grid(loc ~ var) +
  # scale_fill_gradient2(leglab, low = 'blue', mid = 'grey',  high = 'tomato1', midpoint = 0) +
  scale_fill_gradientn(leglab, limits = colrng, colors = c('blue', 'grey', 'tomato1')) +
  scale_color_manual(values = c(scales::alpha('black', 0), 'black'), guide = 'none', drop = FALSE) +
  # coord_map() +
  scale_shape_manual(values = c(24, 25), drop = FALSE, guide = 'none') +
  pthm +
  scale_size(range = c(0.75, 4), guide = 'none') +
  guides(fill = guide_colourbar(barwidth = 0.4, barheight = 2.5)) + 
  labs(
    title = paste0('(a) Change per year, ', yrsel2[1], '-', yrsel2[2]),
    caption = 'Outlines indicate p < 0.05'
  )

p2 <- ktresplo %>% 
  filter(loc == 'Top') %>% 
  pull(plo) %>% 
  .[[1]] + labs(title = '(b) Top, % stations with significant trends by month')

p3 <- ktresplo %>% 
  filter(loc == 'Bottom') %>% 
  pull(plo) %>% 
  .[[1]] + labs(title = '(c) Bottom, % stations with significant trends by month')

p <- p1 + (p2 + p3 + plot_layout(ncol = 1)) + plot_layout(ncol = 2, width = c(0.9, 1))

png(here('figs/suppkendall.png'), height = 7, width = 9.5, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# supp fim temp, sal trends -------------------------------------------------------------------

load(file = here('data/fimsgtempdat.RData'))

toplo <- fimsgtempdat %>% 
  select(date, temp, sal, bay_segment) %>% 
  mutate(
    yr = year(date), 
    mo = month(date)
  ) %>% 
  pivot_longer(temp:sal) %>% 
  summarise(
    avev = mean(value, na.rm = T), 
    lov = t.test(value)$conf.int[1], 
    hiv = t.test(value)$conf.int[2], 
    cnt = n(),
    .by = c(bay_segment, yr, name)
  ) %>% 
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')),
    name = factor(name, levels = c('temp', 'sal'), labels = c('Water temp. (\u00B0C)', 'Salinity (ppt)'))
  )

thm <- theme_minimal() + 
  theme(
    strip.placement = 'outside', 
    panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 11), 
    legend.text = element_text(size= 12), 
    axis.text.x = element_text(size = 7)
  )

p <- ggplot(toplo, aes(x = yr, y = avev)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), show.legend = F, alpha = 0.7, color = 'steelblue4') + 
  geom_point(size = 0.75, color = 'steelblue4') +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  geom_smooth(method = 'lm', formula = y ~ x, se = F, color = 'steelblue4') +
  facet_grid(name ~ bay_segment, switch = 'y', scales = 'free_y') +
  thm +
  theme(
    strip.text = element_text(size = 12)
  ) +
  labs(
    x = NULL, 
    y = NULL,
    color = NULL, 
    shape = NULL
  )

png(here('figs/suppfimtrnds.png'), height = 3.5, width = 7, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()

# mixeff example plot -------------------------------------------------------------------------

load(file = here('data/mixmodprds.RData'))

toplo1 <- mixmodprds %>% 
  select(-mod, -fix, -slo) %>% 
  unnest('data') %>% 
  filter(
    !(cnt > 100 & thrtyp == paste('Temperature >', tmpthr)) # outliers
  ) %>% 
  filter(
    !(cnt > 100 & thrtyp == 'Both') # outliers
  )
toplo2 <- mixmodprds %>% 
  select(-mod, -data, -slo) %>% 
  unnest('fix')
toplo3 <- mixmodprds %>% 
  select(-mod, -data, -fix) %>% 
  unnest('slo') %>% 
  filter(slo != '')

p <- ggplot(toplo1, aes(x = yr, y = cnt)) + 
  geom_point(size = 0.7, alpha = 0.6, color = 'darkgrey') + 
  geom_line(aes(y = prd, group = station), color = 'darkgrey', linewidth = 0.5) + 
  geom_line(data = toplo2, aes(y = prd, col = thrtyp), show.legend = F, linewidth = 1.5) + 
  geom_label(data = toplo3, aes(x = 2022, y = 0, label = slo, col = thrtyp), show.legend = F, 
             vjust = 0, hjust = 1, fontface = 'italic', size = 3, label.r = unit(0, "lines"),
             label.padding = unit(0.1, "lines")) +
  scale_color_manual(values = c('red2', 'dodgerblue2', 'black')) +
  # coord_cartesian(ylim = c(-10, NA)) +
  facet_grid(thrtyp ~ bay_segment, scales = 'free_y') + 
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(), 
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  ) + 
  labs(
    x = NULL, 
    y = expression(paste("Days ", year^-1))
  )

png(here('figs/mixeff.png'), height = 6, width = 7, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# supp mixeff example plot --------------------------------------------------------------------

load(file = here('data/suppmixmodprds.RData'))

toplo1 <- suppmixmodprds %>% 
  select(-mod, -fix, -slo) %>% 
  unnest('data') %>% 
  filter(
    !(cnt > 100 & thrtyp == paste('Temperature >', tmpthr)) # outliers
  ) %>% 
  filter(
    !(cnt > 100 & thrtyp == 'Both') # outliers
  ) %>% 
  filter(yr >= yrsel2[1] & yr <= yrsel2[2])
toplo2 <- suppmixmodprds %>% 
  select(-mod, -data, -slo) %>% 
  unnest('fix')
toplo3 <- suppmixmodprds %>% 
  select(-mod, -data, -fix) %>% 
  unnest('slo') %>% 
  filter(slo != '')

p <- ggplot(toplo1, aes(x = yr, y = cnt)) + 
  geom_point(size = 0.7, alpha = 0.6, color = 'darkgrey') + 
  geom_line(aes(y = prd, group = station), color = 'darkgrey', linewidth = 0.5) + 
  geom_line(data = toplo2, aes(y = prd, col = thrtyp), show.legend = F, linewidth = 1.5) + 
  geom_label(data = toplo3, aes(x = 2022, y = 0, label = slo, col = thrtyp), show.legend = F, 
             vjust = 0, hjust = 1, fontface = 'italic', size = 3, label.r = unit(0, "lines"),
             label.padding = unit(0.1, "lines")) +
  scale_color_manual(values = c('red2', 'dodgerblue2', 'black')) +
  coord_cartesian(xlim = c(1998, 2022)) +
  facet_grid(thrtyp ~ bay_segment, scales = 'free_y') + 
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(), 
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  ) + 
  labs(
    x = NULL, 
    y = expression(paste("Days ", year^-1))
  )

png(here('figs/suppmixeff.png'), height = 6, width = 7, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# nms -----------------------------------------------------------------------------------------

load(file = here('data/cmbdat.RData'))

cmbdat <- cmbdat %>% 
  filter(!bay_segment %in% 'LTB')

toord1 <- cmbdat %>% 
  dplyr::select(-Year, -Chla, -bbave, -Both, -bay_segment) %>% 
  mutate_all(rescale, to = c(0, 1))

orddat1 <- metaMDS(toord1, k = 2, trymax = 100)

toord2 <- cmbdat %>% 
  dplyr::select(-Year, -Chla, -bbave, -Sal, -Temp, -bay_segment) %>% 
  mutate_all(rescale, to = c(0, 1))

orddat2 <- metaMDS(toord2, k = 2, trymax = 100)

grps <- factor(cmbdat$bay_segment)

vec_ext <- 1.75
coord_fix <- F
size <- toord1$total
repel <- F
force <- T
arrow <- 0.2
txt <- 3
alpha <- 0.8
ext <- 1.25
exp <- 0.1
parse <- F
ellipse <- T
cols <- c("#3B9AB2", "#EBCC2A", "#F21A00")# c("#1F78B4", "#33A02C", "#E31A1C")
vec_lab <- list(
  'total'= 'Freq Occ',
  'Sal' = 'Sal', 
  'Temp' = 'Temp', 
  'Both' = 'Both'
)
grp_title <- 'Bay segment'
sizelab <- 'Freq Occ'

p1 <- ggord(orddat1, axes = c('1', '2'), cols = cols, ellipse = ellipse, grp_in = grps,
            parse = parse, vec_ext = vec_ext, coord_fix = coord_fix, size = size, 
            repel = repel, arrow = arrow, txt = txt, alpha = alpha, ext = ext, 
            exp = exp, grp_title = grp_title, sizelab = sizelab, vec_lab = vec_lab, 
            force = force)
p2 <- ggord(orddat2, axes = c('1', '2'), cols = cols, ellipse = ellipse, grp_in = grps,
            parse = parse, vec_ext = vec_ext, coord_fix = coord_fix, size = size, 
            repel = repel, arrow = arrow, txt = txt, alpha = alpha, ext = ext, 
            exp = exp, grp_title = grp_title, sizelab = sizelab, vec_lab = vec_lab, 
            force = force)

p <- p1 + p2 + plot_layout(ncol = 2, guides = 'collect')

png(here('figs/nms.png'), height = 4, width = 9, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# seagrass decline models ---------------------------------------------------------------------

load(file = here('data/pchgmod.RData'))
load(file = here('data/binomod.RData'))

toplo1 <- getprd_fun(pchgmod)
toplo2 <- getprd_fun(binomod)

p1 <- ggplot(toplo1, aes(x = Sal)) + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = 'lightgrey', alpha = 0.5) + 
  geom_line(aes(y = visregFit)) + 
  geom_point(data = tomod, aes(x = Sal, y = pchg), color = 'dodgerblue2') +
  facet_grid(bay_segment ~ yrcat) + 
  labs(
    y = '% change in seagrass', 
    x = '# days salinity < threshold', 
    title = '(a) Percent change in seagrass between years'
  ) + 
  theme_bw() + 
  theme(
    strip.background = element_blank(), 
    panel.grid.minor = element_blank()
  )

p2 <- ggplot(toplo2, aes(x = Sal)) + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = 'lightgrey', alpha = 0.5) + 
  geom_line(aes(y = visregFit)) + 
  geom_rug(data = tomod[tomod$chg == 0,], aes(x = Sal, y = chg), sides = 'b', linewidth = 1, color = 'blue') +
  geom_rug(data = tomod[tomod$chg == 1,], aes(x = Sal, y = chg), sides = 't', linewidth = 1, color = 'red') +
  facet_grid(bay_segment ~ yrcat) + 
  scale_y_continuous(labels = scales::percent) +
  labs(
    y = 'Probability of seagrass decline', 
    x = '# days salinity < threshold', 
    title = '(b) Probability of seagrass decline between years'
  ) + 
  theme_bw() + 
  theme(
    strip.background = element_blank(), 
    panel.grid.minor = element_blank()
  )

p <- p1 + p2 + plot_layout(ncol = 2) & 
  theme(
    strip.text = element_text(size = 12), 
    axis.title = element_text(size = 13), 
    axis.text = element_text(size = 11)
    )

png(here('figs/sgmod.png'), height = 7, width = 9, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# fim gam plots -------------------------------------------------------------------------------

load(file = here("data/fimsalmod.RData"))
load(file = here("data/fimtempmod.RData"))

toplo1 <- getprd_fun2(fimsalmod, 'sal') %>% 
  filter(sal > 10 & sal < 30)
toplo2 <- getprd_fun2(fimtempmod, 'temp')

thm <- theme_bw() + 
  theme(
    panel.grid.minor = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(size = 12), 
    axis.text = element_text(size = 11), 
    axis.title = element_text(size = 12)
  )

# separate plots by bay segment
p1 <- ggplot(toplo1, aes(x = sal, y = visregFit)) + 
  geom_point(data = fimsalmod$model, aes(y = sgcov), size = 0.5, alpha = 0.25) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = 'dodgerblue2', alpha = 0.2) +
  geom_line(color = 'dodgerblue2', size = 1) +
  coord_cartesian(
    ylim = c(0, max(toplo1$visregUpr)),
    xlim = c(10, 30)
  ) +
  facet_grid(bay_segment ~ yrcat, scales = 'free_y') +
  thm + 
  labs(
    x = 'Salinity (psu)', 
    y = 'Seagrass % cover', 
    title = '(a) Seagrass cover vs salinity'
  )

p2 <- ggplot(toplo2, aes(x = temp, y = visregFit)) + 
  geom_point(data = fimtempmod$model, aes(y = sgcov), size = 0.5, alpha = 0.25) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = 'red2', alpha = 0.2) +
  geom_line(color = 'red2', size = 1) +
  coord_cartesian(ylim = c(0, max(toplo2$visregUpr))) +
  facet_grid(bay_segment ~ yrcat, scales = 'free_y') +
  thm + 
  labs(
    x = 'Water temp. (\u00B0 C)', 
    y = 'Seagrass % cover',
    title = '(b) Seagrass cover vs temperature'
  )

p <- p1 + p2 + plot_layout(ncol = 2)

png(here('figs/sggammod.png'), height = 7, width = 9, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# hillsborough flow v temp over time ----------------------------------------------------------

load(file = here('data/gagedat.RData'))

toplo <- gagedat %>%
  filter(name == 'Hillsborough') %>%
  pivot_wider(names_from = 'var', values_from = 'val') %>%
  na.omit() %>%
  mutate(
    yrgroup = ifelse(year < 2000, '1975 - 1999', '2000 - 2022'),
    mo = month(date, label = T, abbr = F)
  ) %>%
  # filter(flow_m3d < quantile(flow_m3d, 0.9)) %>%
  filter(mo %in% c('February', 'August')) #%>%
# filter(temp_c > 1)

p <- ggplot(toplo, aes(x = flow_m3d / 1e6, y = temp_c, group = yrgroup, fill = yrgroup, color = yrgroup)) +
  geom_point(size = 0.6, show.legend = F) +
  scale_x_log10() +
  facet_grid(~mo) +
  theme_minimal() + 
  scale_color_manual(values = c('tomato1', 'darkred')) +
  scale_fill_manual(values = c('tomato1', 'darkred')) +
  theme(
    legend.position = 'top', 
    panel.grid.minor = element_blank(), 
    strip.text = element_text(size = 11)
  ) +
  geom_smooth(se = T, method = 'lm', formula = y ~ x) +
  labs(
    x = expression(paste('Hillsborough River discharge (', 10^6~m^3, ' / yr)')), 
    y = 'Water temp. (\u00B0 C)',
    color = 'Time period',
    fill = 'Time period'
  )

png(here('figs/flowtemp.png'), height = 4, width = 6, family = 'serif', units = 'in', res = 600)
print(p)
dev.off()

# # ports data threshold counts -----------------------------------------------------------------
# 
# load(file = here('data/portsdat.RData'))
# 
# # thresholds to count by year, temp is above
# thrs <- tibble(
#   temp = c(29, 30, 31)
# )
# 
# # get T/F vectors for predictions above/below thresholds by day
# portscnt <- portsdat %>%
#   group_nest(name) %>%
#   crossing(thrs) %>%
#   mutate(
#     data = purrr::pmap(list(data, temp), function(data, temp){
# 
#       out <- data %>%
#         mutate(
#           cnt = temp_c >= temp
#         )
# 
#       return(out)
# 
#     })
#   ) %>%
#   unnest('data') %>%
#   summarise(
#     runcnt = runfunc(cnt),
#     sumcnt = sum(cnt),
#     .by = c(name, yr, temp)
#   )
# 
# # ggplot(portscnt, aes(x = yr, y = runcnt, col = factor(temp))) +
# #   geom_point() +
# #   geom_smooth(se = F, method = 'lm', formula = y ~ x) +
# #   facet_wrap(~name)
# 
# ggplot(portsdat, aes(x = date, y = temp_c)) +
#   geom_line() +
#   facet_wrap(~name, ncol = 1)
# 
# ##
# # epc stations closes to gages (verified from maps)
# 
# epccnts <- tibble(
#     fl = c('MTB_28_tempbot', 'LTB_90_temptop', 'HB_52_temptop', 'OTB_36_temptop'),
#     name = c('stpete', 'manatee', 'mckay', 'oldporttampa')
#   ) %>% 
#   group_nest(name) %>% 
#   mutate(
#     mod = purrr::map(data, function(x){
#   
#       load(file = here(paste0('data/', x, '.RData')))
#       
#       out <- get(x[[1]])
#       
#       return(out)
#       
#     }),
#     prd = purrr::map(mod, anlz_prdday)
#   ) %>% 
#   crossing(
#     temp = c(29, 30, 31)
#   ) %>%
#   mutate(
#     cnt = purrr::pmap(list(prd, temp), function(prd, temp){
# 
#       out <- prd %>% 
#         mutate(
#           cnt = value > temp
#         ) %>% 
#         summarise(
#           modsumcnt = sum(cnt),
#           modruncnt = runfunc(cnt),
#           .by = yr
#         )
#       
#       return(out)
#       
#     })
#   ) %>% 
#   select(-prd, -mod) %>% 
#   unnest(data) %>% 
#   unnest(cnt)
# 
# cmbcnts <- inner_join(portscnt, epccnts, by = c('name', 'temp', 'yr')) 
# 
# toplo <- cmbcnts %>%
#   select(-runcnt, -modruncnt) %>% 
#   pivot_longer(names_to = 'var', values_to = 'val', cols = c(sumcnt, modsumcnt))
# 
# # 1:1 scatterplots  
# ggplot(cmbcnts, aes(x = sumcnt, y = modsumcnt, group = temp, color = factor(temp))) + 
#   geom_point() + 
#   geom_smooth(method = 'lm', se = F) +
#   facet_wrap(~name)
# 
# # line plot by year
# ggplot(toplo, aes(x = yr, y = val, group = var, color = var)) + 
#   geom_point() +
#   # geom_line() + 
#   geom_smooth(method = 'lm', se = F) + 
#   facet_wrap(temp ~ name)
# 
# # comparing ports data to bay segment averages
# data(thrdat)
# 
# modcnt2 <- thrdat %>% 
#   filter(thrtyp == 'tempcnt') %>% 
#   summarise(
#     modcnt = mean(cnt), 
#     .by = c(bay_segment, yr, tempthr)
#   ) %>% 
#   mutate(
#     temp = as.numeric(sub('^temp_', '', tempthr)),
#     name = case_when(
#       bay_segment == 'HB' ~ 'mckay', 
#       bay_segment == 'LTB' ~ 'manatee', 
#       bay_segment == 'OTB' ~ 'oldporttampa', 
#       bay_segment == 'MTB' ~ 'stpete'
#     )
#   ) %>% 
#   select(name, temp, yr, modcnt)
# 
# cmbcnt2 <- inner_join(portscnt, modcnt2, by = c('name', 'temp', 'yr'))
# 
# toplo <- cmbcnt2 %>% 
#   select(-runcnt) %>% 
#   pivot_longer(names_to = 'var', values_to = 'val', sumcnt:modcnt)
# 
# ggplot(toplo, aes(x = yr, y = val, color = var)) + 
#   geom_point() + 
#   geom_line() +
#   facet_wrap(temp~name)
