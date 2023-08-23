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

# raw temp and salinity changes ---------------------------------------------------------------

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
    var = factor(var, levels = c('Sal', 'Temp'), labels = c('Salinity (psu)', 'Temperature (\u00B0 C)')),
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
p <- ggplot(toplo, aes(x = yr, y = avev, group = loc, color = loc)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), position = position_dodge2(width = wd), show.legend = F, alpha = 0.7) + 
  geom_point(position = position_dodge2(width = wd)) +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  geom_smooth(method = 'lm', formula = y ~ x, se = F) +
  facet_grid(var ~ bay_segment, scales = 'free_y', switch = 'y') +
  scale_color_manual(values = c('steelblue1', 'steelblue4')) +
  theme_bw(base_size = 14) + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(), 
    strip.placement = 'outside', 
    legend.position = 'top', 
    strip.background = element_blank(),
    axis.text.x = element_text(colour = 'black', angle = 45, size = 9, hjust = 1)
  ) +
  labs(
    x = NULL, 
    y = NULL, 
    color = NULL, 
    shape = NULL
  )

png(here('figs/saltempraw.png'), height = 5, width = 9, family = 'serif', units = 'in', res = 300)
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
    title = '(a) Change per year, 1974-2022',
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

load(file = here('data/thrdat.RData'))

salthr <- '25'
tmpthr <- '30'

modprd <- thrdat %>% 
  filter(salithr == paste0('sali_', salthr)) %>% 
  filter(tempthr == paste0('temp_', tmpthr)) %>% 
  group_by(thrtyp, bay_segment) %>% 
  nest() %>% 
  mutate(
    mod = purrr::map(data, function(x){
      
      mod <- try(lmer(cnt ~ yr + (1|station), data = x, REML = F), silent = T)
      if(inherits(mod, 'try-error'))
        return(NA)
      
      return(mod)
      
    }), 
    data = purrr::pmap(list(data, mod), function(data, mod){
      bind_cols(data, prd = predict(mod))
    }),
    fix = purrr::map(mod, function(x){
      
      fixef <- estimate_means(x, 'yr')
      
      out <- tibble(
        yr = fixef$yr, 
        prd = fixef$Mean
      )
      
      return(out)
      
    }),
    slo = purrr::map(mod, function(x){
      
      summod <- summary(x)$coefficients
      
      pvl <- summod['yr', 'Pr(>|t|)']
      if(pvl >= 0.05)
        return('')
      
      out <- summod['yr', 'Estimate'] %>% 
        round(2) %>% 
        as.character()
      
      return(out)
      
    })
  ) %>% 
  ungroup() %>% 
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')), 
    thrtyp = factor(thrtyp, 
                    levels = c('salicnt', 'tempcnt', 'bothcnt'), 
                    labels = c(paste('Salinity <', salthr), paste('Temperature >', tmpthr), 'Both')
    )
  )

toplo1 <- modprd %>% 
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
