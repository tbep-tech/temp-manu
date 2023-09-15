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

yrsel <- c(1998, 2022)
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

p1 <- ggplot(toplo1, aes(x = yr, y = precip_mm / 1e3)) + 
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

p2 <- ggplot(hydrodat, aes(x = year, y = hy_load)) +
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

p3 <- ggplot(speidat, aes(x = date, y = spi, fill = spisign)) + 
  geom_col(width = 35) +
  scale_fill_manual(values = c('grey70', 'grey30'), guide = 'none') + 
  thm + 
  theme(
    axis.text.x = element_blank()
  ) +
  coord_cartesian(xlim = range(speidat$date)) +
  labs(
    x = NULL, 
    y = 'SPI (z-values)'
  )

p4 <- ggplot(toplo1, aes(x = yr, y = tavg_c)) + 
  geom_line() + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F, formula = y~x, color = 'red2') +
  thm +
  theme(
    axis.text.x = element_text(size = 12)
  ) +
  labs(
    x = NULL, 
    y = 'Air temp. (\u00B0 C)'
  )

toplo <- epcdata %>% 
  select(bay_segment, epchc_station, SampleTime, yr, matches('Top|Bottom')) %>% 
  filter(yr >= yrsel[1] & yr <= yrsel[2]) %>% 
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
    var = factor(var, levels = c('Sal', 'Temp'), labels = c('Salinity (ppt)', 'Water temp. (\u00B0C)')),
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
  filter(var == 'Salinity (ppt)')
p5 <- ggplot(toplo5, aes(x = yr, y = avev, group = loc, color = loc)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), position = position_dodge2(width = wd), show.legend = F, alpha = 0.7) + 
  geom_point(position = position_dodge2(width = wd), size = 0.5) +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  geom_smooth(method = 'lm', formula = y ~ x, se = F) +
  facet_grid(~ bay_segment) +
  scale_color_manual(values = c('steelblue1', 'steelblue4')) +
  thm +
  theme(
    axis.text.x = element_blank(), 
    strip.text = element_text(size = 12)
  ) +
  labs(
    x = NULL, 
    y = 'Salinity (ppt)', 
    color = NULL, 
    shape = NULL
  )
toplo6 <- toplo %>% 
  filter(var == 'Water temp. (\u00B0C)') %>% 
  filter(!(yr %in% c(1982, 1985) & bay_segment == 'OTB')) %>%  # missing months create outliers
  filter(!(yr %in% 1975 & bay_segment == 'HB')) # missing months create outliers
p6 <- ggplot(toplo6, aes(x = yr, y = avev, group = loc, color = loc)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), position = position_dodge2(width = wd), show.legend = F, alpha = 0.7) + 
  geom_point(position = position_dodge2(width = wd), size = 0.5) +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  geom_smooth(method = 'lm', formula = y ~ x, se = F) +
  facet_grid(~ bay_segment) +
  scale_color_manual(values = c('steelblue1', 'steelblue4')) +
  thm +
  theme(strip.text = element_blank()) +
  labs(
    x = NULL, 
    y = 'Water temp. (\u00B0C)', 
    color = NULL, 
    shape = NULL
  )

p <- p1 + p2 + p3 + p4 + (p5 + p6 + plot_layout(ncol = 1, guides = 'collect')) + plot_layout(ncol = 1, heights = c(1, 1, 1, 1, 2.5)) & theme(legend.position = 'top')

png(here('figs/meteowqraw.png'), height = 8.5, width = 7, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()

# kendall -------------------------------------------------------------------------------------

leglab <- expression(paste(yr^{-1}))

# kendall all years
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

#kendall by month
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

# plots of segment summaries by month
ktresplo <- ktres %>% 
  mutate(
    mo = month(mo, label = T), 
    bay_segment = factor(bay_segment, levels = c('LTB', 'MTB', 'HB', 'OTB')), 
    loc =  gsub("^.*_(.*)_.*$", "\\1", var),
    varsimp = gsub('^(.*?)_.*$', '\\1', var, perl = T),
    varsimp = factor(varsimp, levels = c('Sal', 'Temp'), labels = c('Salinity', 'Temperature')),
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
    var = factor(var, levels = c('Sal', 'Temp'), labels = c('Salinity', 'Temperature')),
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
    title = paste0('(a) Change per year, ', yrsel[1], '-', yrsel[2]),
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
toplo2 <- modprd %>% 
  select(-mod, -data, -slo) %>% 
  unnest('fix')
toplo3 <- modprd %>% 
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
  scale_color_manual(values = c('dodgerblue2', 'red2', 'black')) +
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

# gam results ---------------------------------------------------------------------------------

load(file = here('data/cmbmod.RData'))

# load(file = here('data/cmbdat.RData'))
# 
# tomod <- cmbdat %>%
#   filter(!bay_segment %in% 'LTB')
# 
# vifmod <- glm(total ~ Sal + Temp + Both + Year, data = tomod)
# vif(vifmod)

yrbrks <- c(2000, 2010, 2016, 2020, 2022)
# visreg(cmbmod, 'Both', by = 'Year', breaks = yrbrks, scale = 'response', rug = 1)
yrvis <- visreg(cmbmod, 'Year', scale = 'response', rug = 1, plot = F)
yrfit <- yrvis$fit
yrres <- yrvis$res
tempvis <- visreg(cmbmod, 'Temp', by = 'Year', breaks = yrbrks, scale = 'response', rug = 1, plot = F)
tempfit <- tempvis$fit
tempres <- tempvis$res
salivis <- visreg(cmbmod, 'Sal', by = 'Year', breaks = yrbrks, scale = 'response', rug = 1, plot = F)
salifit <- salivis$fit
salires <- salivis$res

alph <- 0.5
ylb <- 'Freq. Occ.'
lwd <- 0.8
ptsz <- 1
thm <- theme_minimal() + 
  theme(
    panel.grid.minor = element_blank(), 
    strip.text = element_text(size = 11),
    panel.spacing.x = unit(0.5, 'lines')
  )
salicol <- 'dodgerblue2'
tempcol <- 'red2'
ylim <- c(0, 1)

p1 <- ggplot(yrfit, aes(x = Year)) + 
  geom_point(data = yrres, aes(y = visregRes), size = ptsz) + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = alph) + 
  geom_line(aes(y = visregLwr), linetype = 'dashed', linewidth = lwd) +
  geom_line(aes(y = visregUpr), linetype = 'dashed', linewidth = lwd) +
  geom_line(aes(y = visregFit), linewidth = lwd) + 
  coord_cartesian(ylim = ylim) +
  thm +
  labs(
    y = ylb, 
    x = 'Year'
  )

p2 <- ggplot(tempfit, aes(x = Temp)) + 
  geom_point(data = tempres, aes(y = visregRes), size = ptsz) + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = alph, color = tempcol, fill = tempcol) + 
  geom_line(aes(y = visregLwr), linetype = 'dashed', linewidth = lwd) +
  geom_line(aes(y = visregUpr), linetype = 'dashed', linewidth = lwd) +
  geom_line(aes(y = visregFit), linewidth = lwd) + 
  coord_cartesian(ylim = ylim) +
  facet_wrap(~Year, ncol = length(unique(tempres$Year))) +
  scale_x_continuous(n.breaks = 3) +
  thm +
  theme(axis.text.x = element_text(size = 9)) +
  labs(
    y = ylb, 
    x = '# days temperature > threshold'
  )

p3 <- ggplot(salifit, aes(x = Sal)) + 
  geom_point(data = salires, aes(y = visregRes), size = ptsz) + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = alph, color = salicol, fill = salicol) + 
  geom_line(aes(y = visregLwr), linetype = 'dashed', linewidth = lwd) +
  geom_line(aes(y = visregUpr), linetype = 'dashed', linewidth = lwd) +
  geom_line(aes(y = visregFit), linewidth = lwd) + 
  coord_cartesian(ylim = ylim) +
  facet_wrap(~Year, ncol = length(unique(salires$Year))) +
  scale_x_continuous(n.breaks = 3) +
  thm +
  theme(axis.text.x = element_text(size = 9)) +
  labs(
    y = ylb, 
    x = '# days salinity < threshold'
  )

p <- p1 + (p2 + p3 + plot_layout(ncol = 1)) + plot_layout(ncol = 2, widths = c(0.5, 1))

png(here('figs/gamres.png'), height = 5, width = 8.5, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# # hillsborough flow v temp over time ----------------------------------------------------------
# 
# load(file = here('data/gagedat.RData'))
# 
# toplo <- gagedat %>% 
#   filter(name == 'Hillsborough') %>% 
#   pivot_wider(names_from = 'var', values_from = 'val') %>% 
#   na.omit() %>% 
#   mutate(
#     yrgroup = ifelse(year < 2000, 'past', 'present'),
#     mo = month(date)
#   ) %>% 
#   # filter(flow_m3d < quantile(flow_m3d, 0.9)) %>%
#   filter(mo %in% c(2, 8)) #%>%
# # filter(temp_c > 1)
# 
# ggplot(toplo, aes(x = flow_m3d, y = temp_c, group = yrgroup, color = yrgroup)) + 
#   geom_point(size = 0.6) +  
#   scale_x_log10() +
#   facet_grid(~mo) +
#   geom_smooth(se = T, method = 'lm', formula = y~ x)

# # ports data threshold counts -----------------------------------------------------------------
# 
# # thresholds to count by year, temp is above
# thrs <- tibble(
#   temp = c(29, 30, 31)
# )
# 
# # get T/F vectors for predictions above/below thresholds by day
# sumdat <- portsdat %>% 
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
# ggplot(sumdat, aes(x = yr, y = runcnt, col = factor(temp))) + 
#   geom_point() + 
#   geom_smooth(se = F, method = 'lm', formula = y ~ x) + 
#   facet_wrap(~name)
# 
# ggplot(portsdat, aes(x = date, y = temp_c)) + 
#   geom_line() + 
#   facet_wrap(~name, ncol = 1)



# # alafia hydro load v s2t3 fo chng ------------------------------------------------------------
#
# data(gagedat)
# 
# toplo <- gagedat %>%
#   filter(var == 'flow_m3d') %>%
#   summarise(
#     # loqnt = quantile(val, 0.99),
#     meqnt = quantile(val, 0.90),
#     # hiqnt = quantile(val, 0.85),
#     .by = c(name, year)
#   ) %>%
#   pivot_longer(names_to = 'met', values_to = 'val', matches('qnt')) %>%
#   filter(name %in% c('Alafia', 'Hillsborough'))
# 
# ggplot(toplo, aes(x = year, y = val, color = met)) +
#   geom_point() +
#   geom_smooth(method = 'lm', se = F) +
#   facet_grid(~name, scales = 'free_y')
# # scale_y_log10()
# 
# 
# trnptssub <- trnpts %>%
#   st_set_geometry(NULL) %>%
#   select(
#     Transect = TRAN_ID,
#     bay_segment
#   ) %>%
#   filter(bay_segment == 'HB')
# transectocc <- anlz_transectocc(transect)
# tots <- transectocc %>%
#   filter(Savspecies == 'total') %>%
#   inner_join(trnptssub, by = 'Transect') %>%
#   filter(Transect %in% c('S2T5', 'S2T4', 'S2T3', 'S2T2'))
# 
# ggplot(tots, aes(x = Date, y = foest, group = Transect)) +
#   geom_line() +
#   geom_point()
# facet_grid(~Transect)
# 
# s2t3 <- tots %>%
#   ungroup() %>%
#   filter(Transect == 'S2T3') %>%
#   select(Date, foest, bbest)
# 
# alaf <- gagedat %>%
#   filter(name == 'Alafia') %>%
#   filter(year >= min(year(s2t3$Date))) %>%
#   filter(var == 'flow_m3d') %>%
#   mutate(
#     trndt = cut(as.numeric(as.Date(date)), c(-Inf, unique(s2t3$Date), Inf),
#                 labels = c(as.character(unique(s2t3$Date)), 'beyond')),
#   ) %>%
#   filter(!trndt %in% 'beyond') %>%
#   mutate(trndt = ymd(trndt)) %>%
#   mutate(
#     mocnt = trndt - floor_date(as.Date(date), 'month'),
#     mocnt = factor(mocnt, levels = unique(mocnt), labels = rev(c(1:length(unique(mocnt))))),
#     mocnt = as.numeric(as.character(mocnt)),
#     .by = trndt
#   )
# 
# for(i in 1:12){
# 
#   alaf2 <- alaf %>%
#     filter(mocnt <= i) %>%
#     summarise(
#       ave = exp(mean(log(val))),
#       ld = sum(val),
#       qnt90 = quantile(val, 0.9),
#       .by = 'trndt'
#     )
# 
#   tomod <- s2t3 %>%
#     ungroup() %>%
#     full_join(alaf2, by = c('Date' = 'trndt')) %>%
#     mutate(
#       fodif = c(NA, diff(foest)),
#       fochg = ifelse(sign(fodif) == -1, 1, 0),
#       bbdif = c(NA, diff(bbest)),
#       bbchg = ifelse(sign(bbdif) == -1, 1, 0)
#     ) %>%
#     filter(!is.na(fochg))
# 
#   mod <- glm(fochg ~ ld, data = tomod, family = binomial('logit'))
#   # visreg(mod, 'ld', scale = 'response')
#   print(coefficients(summary(mod))[2, 4])
# 
# }
# 
# alaf2 <- alaf %>% 
#   filter(mocnt <= 8) %>%
#   summarise(
#     ave = exp(mean(log(val))),
#     ld = sum(val),
#     qnt90 = quantile(val, 0.9), 
#     .by = 'trndt'
#   )
# 
# tomod <- s2t3 %>% 
#   ungroup() %>% 
#   full_join(alaf2, by = c('Date' = 'trndt')) %>% 
#   mutate(
#     fodif = c(NA, diff(foest)),
#     fochg = ifelse(sign(fodif) == -1, 1, 0), 
#     bbdif = c(NA, diff(bbest)), 
#     bbchg = ifelse(sign(bbdif) == -1, 1, 0)
#   ) %>% 
#   filter(!is.na(fochg))
# 
# mod <- glm(fochg ~ ld, data = tomod, family = binomial('logit'))
# visreg(mod, 'ld', scale = 'response')
# 
# ggplot(tomod, aes(x = factor(Date), y = ld, fill = fodif)) + 
#   geom_col(color = 'black') +
#   scale_fill_gradient2(low = 'red', mid = 'white', high = 'green', midpoint = 0)
# 
# ggplot(tomod, aes(x = factor(Date), y = ld, fill = factor(fochg))) + 
#   geom_col(color = 'black')

# bay segment hydrologic load v foest ---------------------------------------------------------

load(file = url('https://github.com/tbep-tech/load-estimates/raw/main/data/mohydat.RData'))

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
    month = floor(mean(month)), 
    .by = c('bay_segment', 'year')
  ) %>% 
  mutate(trndt = paste0('trndt', 1:n()), .by = 'bay_segment')

lddat <- mohydat %>% 
  arrange(bay_segment, year, month) %>% 
  mutate(
    bay_segment = case_when(
      bay_segment == 'Old Tampa Bay' ~ 'OTB',
      bay_segment == 'Hillsborough Bay' ~ 'HB',
      bay_segment == 'Middle Tampa Bay' ~ 'MTB',
      T ~ 'remove'
    )
  ) %>% 
  filter(!bay_segment == 'remove') %>% 
  filter(year >= (min(trndts$year) - 1)) %>% 
  left_join(trndts, by = c('year', 'month', 'bay_segment')) %>% 
  group_by(bay_segment) %>% 
  fill(trndt, .direction = 'up') %>% 
  ungroup() %>% 
  mutate(
    mocnt = rev(1:n()), 
    trnyear = max(year),
    .by = c('trndt', 'bay_segment')
  ) %>% 
  filter(!is.na(trndt))

for(i in 1:12){
  
  lddat2 <- lddat %>%
    filter(mocnt <= i) %>%
    summarise(
      ave = mean(hy_load_106_m3_mo),
      ld = sum(hy_load_106_m3_mo ),
      qnt90 = quantile(hy_load_106_m3_mo, 0.9),
      .by = c('trnyear', 'bay_segment')
    )
  
  tomod <- transectavespp %>%
    ungroup() %>%
    inner_join(lddat2, by = c('yr' = 'trnyear', 'bay_segment')) %>%
    arrange(bay_segment, yr) %>% 
    mutate(
      fodif = c(NA, diff(foest)),
      fochg = ifelse(sign(fodif) == -1, 1, 0), 
      .by = bay_segment
    ) %>%
    filter(!is.na(fochg))
  
  # mod <- glm(fochg ~ ld, data = tomod, family = binomial('logit'))
  mod <- glm(foest ~ ld, data = tomod)
  visreg(mod, 'ld', scale = 'response')
  print(coefficients(summary(mod))[2, 4])
  
}

lddat2 <- lddat %>%
  filter(mocnt <= 12) %>%
  summarise(
    ave = mean(hy_load_106_m3_mo),
    ld = sum(hy_load_106_m3_mo ),
    qnt90 = quantile(hy_load_106_m3_mo, 0.9),
    .by = c('trnyear', 'bay_segment')
  )

tomod <- transectavespp %>%
  ungroup() %>%
  inner_join(lddat2, by = c('yr' = 'trnyear', 'bay_segment')) %>%
  arrange(bay_segment, yr) %>% 
  mutate(
    fodif = c(NA, diff(foest)),
    fochg = ifelse(sign(fodif) == -1, 1, 0), 
    .by = bay_segment
  ) %>%
  filter(!is.na(fochg))

ggplot(tomod, aes(x = ld, y = foest)) + 
  geom_point() +
  geom_smooth(method = 'lm', formula = y~x) 


# ---------------------------------------------------------------------------------------------

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

dys <- seq(30, 365, by = 30)
for(i in dys){
  
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
  
  ##
  # combine
  
  cmbdat <- fodat %>% 
    inner_join(thrtrndat, by = c('yr', 'bay_segment')) %>% 
    inner_join(bbdat, by = c('yr', 'bay_segment')) %>% 
    mutate(
      bay_segment = factor(bay_segment, 
                           levels = segshr)
    ) %>% 
    pivot_wider(names_from = 'thrtyp', values_from = 'cnt') %>% 
    rename(
      Year = yr,
      Sal = salicnt, 
      Temp = tempcnt, 
      Both = bothcnt
    ) %>% 
    arrange(bay_segment, Year) %>% 
    mutate(
      chng = c(NA, sign(diff(total))), 
      chng = ifelse(chng == -1, 1, 0), 
      .by = 'bay_segment'
    ) %>% 
    filter(!is.na(chng))
  
  tomod <- cmbdat %>% 
    mutate(
      pchg = (total - lag(total)) / lag(total), 
      .by = bay_segment
    ) %>% 
    filter(!is.na(pchg)) %>% 
    filter(bay_segment != 'LTB')
  
  combmod <- lm(pchg ~ bay_segment*Sal + Year*Sal + bay_segment*Temp + Year*Temp, data = tomod) %>% 
    step(trace = 0)
  
  cat(i, '\n')
  # print(summary(combmod)$s.table)
  print(coefficients(summary(combmod)))
  
}

# use i = 360

yrbrks <- c(2000, 2010, 2016, 2020, 2022)

visreg(combmod, 'Sal', by = 'Year', breaks = yrbrks, cond = list(bay_segment = 'HB'))
visreg(combmod, 'Temp', by = 'Year', breaks = yrbrks)

cmbmod <- gam(pchg ~ ti(Temp, by = bay_segment) + ti(Temp, Year), data = tomod)

visreg(cmbmod, 'Temp', by = 'Year', breaks = yrbrks, cond = list(bay_segment = 'HB'))

cmbmod <- gam(pchg ~ ti(Sal, by = bay_segment) + ti(Sal, Year), data = tomod)

visreg(cmbmod, 'Sal', by = 'Year', breaks = yrbrks, cond = list(bay_segment = 'OTB'))
