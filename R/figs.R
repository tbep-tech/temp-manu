# setup ---------------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(EnvStats)
library(tbeptools)
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
library(mixmeta)
library(maps)

source(here('R/funcs.R'))

seglng <- c('Old Tampa Bay', 'Hillsborough Bay', 'Middle Tampa Bay', 'Lower Tampa Bay')
segshr <- c('OTB', 'HB', 'MTB', 'LTB')

yrsel1 <- c(1975, 2022)
yrsel2 <- c(1998, 2022)
tempthr <- c(29, 30, 31)
tempcol <- c('coral', 'red2', 'darkred')
salithr <- c(20, 25, 30)
salicol <- c('navyblue', 'dodgerblue2', 'slategray3')

# use mid salinity for 1975, bottom missing
data(epcdata)
epcdat <- epcdata %>% 
  mutate(
    Sal_Bottom_ppth = case_when(
      yr == 1975 & is.na(Sal_Bottom_ppth) ~ Sal_Mid_ppth, 
      T ~ Sal_Bottom_ppth
    )
  ) %>% 
  filter(yr < 2023)

# map -----------------------------------------------------------------------------------------

fl <- paste0(tempdir(), '/sgdat2022.RData')
download.file('https://github.com/tbep-tech/hmpu-workflow/raw/master/data/sgdat2022.RData', destfile = fl)
load(file = fl)

data(pincotemp)
data(fimsgtempdat)

sgdat <- sgdat2022 %>% 
  filter(FLUCCSCODE %in% c(9113, 9116)) %>% 
  st_simplify(5, preserveTopology = F)

flpoly <- map_data('state', 'florida') %>% 
  st_as_sf(coords = c('long', 'lat'), crs = 4326) %>% 
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

bbox <- st_bbox(tbseg)

minset <- ggplot() + 
  geom_sf(data = flpoly, fill = 'grey', color = NA) + 
  geom_sf(data = st_as_sfc(bbox), fill = NA, color = 'black', linewidth = 0.5) + 
  theme_void() +
  theme( 
    panel.background = element_rect(fill = '#FFFFFF', colour = 'white'), 
    panel.border = element_rect(colour = 'black', fill = 'transparent')
  ) 

segcent <- tbseg %>% 
  st_centroid() 

thm <- theme(
  panel.grid = element_blank(), 
  axis.title = element_blank(), 
  axis.text.y = element_text(size = 6),
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6), 
  axis.ticks = element_blank()
)

m1 <- ggplot() + 
  ggspatial::annotation_map_tile(zoom = 11, type = 'cartolight', cachedir = system.file("rosm.cache", package = "ggspatial")) +
  geom_sf(data = sgdat, fill = 'darkgreen', color = NA, inherit.aes = F) +
  geom_sf(data = trnpts, color = 'black', inherit.aes = F) +
  annotation_north_arrow(location = 'tl', style = north_arrow_orienteering(fill = c('black', 'black'), text_col = NA), 
                         height = unit(0.5, "cm"), width = unit(0.5, "cm")) +
  annotation_scale(location = 'br', text_cex = 1) +
  geom_sf(data = tbseg, fill = NA, color = NA, inherit.aes = F) +
  geom_sf(data = tbseglines, color = 'black', inherit.aes = F) +
  geom_sf_text(data = segcent, aes(label = bay_segment), size = 4, color = 'black', inherit.aes = F) +
  # annotation_custom(ggplotGrob(minset), xmin = -9.185e6, xmax = -9.17e6, ymin = 3.22e6, ymax = 3.28e6) + 
  annotation_custom(ggplotGrob(minset), xmin = bbox[3] - 0.1, xmax = bbox[3] + 0.015, ymin = bbox[4] - 0.1, ymax = bbox[4] + 0.06) + 
  coord_sf(xlim = bbox[c('xmin', 'xmax')], ylim = bbox[c('ymin', 'ymax')], crs = 4326) +
  labs(
    subtitle = '(a) Bay segments, seagrass'
  ) +
  thm + 
  theme(
    axis.text.x = element_blank()
  )

# xnrg <- ggplot_build(m1)$layout$panel_scales_x[[1]]$range$range
# yrng <- ggplot_build(m1)$layout$panel_scales_y[[1]]$range$range

epcpts <- epcdat %>% 
  select(long = Longitude, lat = Latitude) %>% 
  unique() %>% 
  st_as_sf(coords = c('long', 'lat'), crs = 4326) 

m2 <- ggplot() + 
  ggspatial::annotation_map_tile(zoom = 11, type = 'cartolight', cachedir = system.file("rosm.cache", package = "ggspatial")) +
  # geom_sf(data = sgdat, fill = 'darkgreen', color = NA, inherit.aes = F) +
  geom_sf(data = tbseg, fill = NA, color = NA, inherit.aes = F) +
  geom_sf(data = epcpts, color = 'black', inherit.aes = F) + 
  coord_sf(xlim = bbox[c('xmin', 'xmax')], ylim = bbox[c('ymin', 'ymax')], crs = 4326) +
  labs(
    subtitle = '(b) EPC, 1975-2022'
  ) +
  thm +  
  theme(
    axis.text.y = element_blank(), 
    axis.text.x = element_blank()
  )

fimpts <- fimsgtempdat %>% 
  select(long = lon, lat) %>% 
  unique() %>% 
  st_as_sf(coords = c('long', 'lat'), crs = 4326) 

m3 <- ggplot() + 
  ggspatial::annotation_map_tile(zoom = 11, type = 'cartolight', cachedir = system.file("rosm.cache", package = "ggspatial")) +
  # geom_sf(data = sgdat, fill = 'darkgreen', color = NA, inherit.aes = F) +
  geom_sf(data = tbseg, fill = NA, color = NA, inherit.aes = F) +
  geom_sf(data = fimpts, color = 'black', inherit.aes = F, size = 0.5, alpha = 0.5) + 
  coord_sf(xlim = bbox[c('xmin', 'xmax')], ylim = bbox[c('ymin', 'ymax')], crs = 4326) +
  labs(
    subtitle = '(c) FIM, 1996-2022'
  ) +
  thm

pdempts <- pincotemp %>% 
  select(long = lon, lat) %>% 
  unique() %>% 
  st_as_sf(coords = c('long', 'lat'), crs = 4326)

m4 <- ggplot() + 
  ggspatial::annotation_map_tile(zoom = 11, type = 'cartolight', cachedir = system.file("rosm.cache", package = "ggspatial")) +
  # geom_sf(data = sgdat, fill = 'darkgreen', color = NA, inherit.aes = F) +
  geom_sf(data = tbseg, fill = NA, color = NA, inherit.aes = F) +
  geom_sf(data = pdempts, color = 'black', inherit.aes = F, size = 0.5, alpha = 0.5) + 
  coord_sf(xlim = bbox[c('xmin', 'xmax')], ylim = bbox[c('ymin', 'ymax')], crs = 4326) +
  labs(
    subtitle = '(d) PDEM, 2004-2022'
  ) +
  thm + 
  theme(
    axis.text.y = element_blank()
  )

m <- m1 + m2 + m3 + m4 + plot_layout(ncol = 2)

png(here('figs/map.png'), height = 6.5, width = 4, family = 'serif', units = 'in', res = 300)
print(m)
dev.off()

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
  select(yr, tavg_c) %>% 
  summarise(
    tavg_c = mean(tavg_c), 
    .by = 'yr'
  )

toplo2 <- speidat %>% 
  filter(mo %in% c(6:8)) %>%
  select(yr, precip_mm) %>% 
  summarise(
    precip_mm = sum(precip_mm), 
    .by = 'yr'
  )

thm <- theme_minimal() + 
  theme(
    strip.placement = 'outside', 
    panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 11), 
    axis.title.y = element_text(size = 10),
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
    y = 'Air temp. (\u00B0C)'
  )

p2 <- ggplot(toplo2, aes(x = yr, y = precip_mm / 1e3)) + 
  geom_line() + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F, formula = y~x, color = 'blue') +
  thm +
  theme(
    axis.text.x = element_blank()
  ) +
  labs(
    x = NULL, 
    y = 'Jun-Aug prec. (m)'
  )

# p3 <- ggplot(hydrodat, aes(x = year, y = hy_load)) +
#   geom_line() + 
#   geom_point() + 
#   geom_smooth(method = 'lm', se = F, formula = y~x, color = 'blue') +
#   coord_cartesian(xlim = range(speidat$yr)) + 
#   thm +
#   theme(
#     axis.text.x = element_blank()
#   ) +
#   labs(
#     x = NULL, 
#     y = expression(paste('Hy. load (', 10^3, ' t/y)'))
#   )

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

toplo <- epcdat %>% 
  select(bay_segment, epchc_station, SampleTime, yr, matches('Top|Bottom')) %>% 
  pivot_longer(names_to = 'var', values_to = 'val', matches('Top|Bottom')) %>% 
  mutate(
    var = factor(var, 
                 levels = c(c("Sal_Top_ppth", "Sal_Bottom_ppth", "Temp_Water_Top_degC", "Temp_Water_Bottom_degC"
                 )), 
                 labels = c("Sal_Top", "Sal_Bottom", "Temp_Top", "Temp_Bottom")
                 ),
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) %>% 
  separate(var, c('var', 'loc'), remove = T) %>% 
  mutate(
    var = factor(var, levels = c('Temp', 'Sal'), labels = c('Water temp. (\u00B0C)', 'Salinity (ppt)')),
    loc = factor(loc, levels = c('Top', 'Bottom'))
  ) %>% 
  filter(!is.na(val)) %>% 
  summarise(
    avev = mean(val, na.rm = T),
    se = sd(val, na.rm = T)^2,
    lov = t.test(val, na.rm = T)$conf.int[1],
    hiv = t.test(val, na.rm = T)$conf.int[2],
    .by = c('bay_segment', 'yr', 'var', 'loc') 
  )

# # get mixed-effect meta-analysis linear preds
# mixmets <- toplo %>% 
#   group_nest(bay_segment, var, loc) %>% 
#   mutate(
#     data = purrr::map(data, function(x){
#       
#       mod <- mixmeta(avev ~ yr, S = se, random = ~1|yr, data = x, method = 'reml')
#       
#       out <- x %>% 
#         mutate(
#           prd = predict(mod)
#         ) %>% 
#         select(yr, prd)
#       
#       return(out)
#       
#     }) 
#   ) %>% 
#   unnest(data)

wd <- 0.5

toplo5 <- toplo %>% 
  filter(var == 'Water temp. (\u00B0C)') %>% 
  filter(!(yr %in% c(1982, 1985) & bay_segment == 'OTB')) %>%  # missing months create outliers
  filter(!(yr %in% 1975 & bay_segment == 'HB')) # missing months create outliers
p5 <- ggplot(toplo5, aes(x = yr, y = avev, group = loc, color = loc)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), position = position_dodge2(width = wd), show.legend = F, alpha = 0.7) + 
  geom_point(position = position_dodge2(width = wd), size = 0.5) +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  geom_smooth(method = 'lm', se = F, formula = y ~ x) +
  # geom_line(data = mixmets %>% filter(var == 'Water temp. (\u00B0C)'), aes(y = prd)) +
  facet_grid(~ bay_segment) +
  scale_color_manual(values = c( 'red2', 'red4')) +
  thm +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_blank()
  ) +
  labs(
    x = NULL, 
    y = '\u00B0C', 
    color = 'Water temp.', 
    shape = NULL
  )
toplo6 <- toplo %>% 
  filter(var == 'Salinity (ppt)')
p6 <- ggplot(toplo6, aes(x = yr, y = avev, group = loc, color = loc)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), position = position_dodge2(width = wd), show.legend = F, alpha = 0.7) + 
  geom_point(position = position_dodge2(width = wd), size = 0.5) +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  # geom_line(data = mixmets %>% filter(var == 'Salinity (ppt)'), aes(y = prd)) +
  geom_smooth(method = 'lm', se = F, formula = y ~ x) +
  facet_grid(~ bay_segment) +
  scale_color_manual(values = c('dodgerblue1', 'dodgerblue4')) +
  thm +
  theme(
    strip.text = element_blank()
  ) +
  labs(
    x = NULL, 
    y = 'ppt', 
    color = 'Salinity', 
    shape = NULL
  )

p <- p1 + p2 + p4 + (p5 + p6 + plot_layout(ncol = 1, guides = 'collect')) + plot_layout(ncol = 1, heights = c(1, 1, 1, 2.75)) & theme(legend.position = 'top')

png(here('figs/meteowqraw.png'), height = 7.5, width = 7, family = 'serif', units = 'in', res = 500)
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
  select(yr, tavg_c) %>% 
  summarise(
    tavg_c = mean(tavg_c), 
    .by = 'yr'
  )

toplo2 <- speidat %>% 
  filter(mo %in% c(6:8)) %>% 
  select(yr, precip_mm) %>% 
  summarise(
    precip_mm = sum(precip_mm), 
    .by = 'yr'
  )

thm <- theme_minimal() + 
  theme(
    strip.placement = 'outside', 
    panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 11), 
    axis.title.y = element_text(size = 10),
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

p2 <- ggplot(toplo2, aes(x = yr, y = precip_mm / 1e3)) + 
  geom_line() + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F, formula = y~x, color = 'blue') +
  thm +
  theme(
    axis.text.x = element_blank()
  ) +
  labs(
    x = NULL, 
    y = 'Jun-Aug prec. (m)'
  )

# p3 <- ggplot(hydrodat, aes(x = year, y = hy_load)) +
#   geom_line() + 
#   geom_point() + 
#   geom_smooth(method = 'lm', se = F, formula = y~x, color = 'blue') +
#   coord_cartesian(xlim = range(toplo1$yr)) + 
#   thm +
#   theme(
#     axis.text.x = element_blank()
#   ) +
#   labs(
#     x = NULL, 
#     y = expression(paste('Hy. load (', 10^3, ' t/y)'))
#   )

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

toplo <- epcdat %>% 
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
    se = sd(val, na.rm = T)^2,
    lov = t.test(val, na.rm = T)$conf.int[1],
    hiv = t.test(val, na.rm = T)$conf.int[2],
    .by = c('bay_segment', 'yr', 'var', 'loc') 
  )

# # get mixed-effect meta-analysis linear preds
# mixmets <- toplo %>% 
#   group_nest(bay_segment, var, loc) %>% 
#   mutate(
#     data = purrr::map(data, function(x){
#       
#       mod <- mixmeta(avev ~ yr, S = se, random = ~1|yr, data = x, method = 'reml')
#       
#       out <- x %>% 
#         mutate(
#           prd = predict(mod)
#         ) %>% 
#         select(yr, prd)
#       
#       return(out)
#       
#     }) 
#   ) %>% 
#   unnest(data)

wd <- 0.5

toplo5 <- toplo %>% 
  filter(var == 'Water temp. (\u00B0C)') %>% 
  filter(!(yr %in% c(1982, 1985) & bay_segment == 'OTB')) %>%  # missing months create outliers
  filter(!(yr %in% 1975 & bay_segment == 'HB')) # missing months create outliers
p5 <- ggplot(toplo5, aes(x = yr, y = avev, group = loc, color = loc)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), position = position_dodge2(width = wd), show.legend = F, alpha = 0.7) + 
  geom_point(position = position_dodge2(width = wd), size = 0.5) +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  # geom_line(data = mixmets %>% filter(var == 'Water temp. (\u00B0C)'), aes(y = prd)) + 
  geom_smooth(method = 'lm', se = F, formula = y ~ x) +
  facet_grid(~ bay_segment) +
  scale_color_manual(values = c( 'red2', 'red4')) +
  thm +
  theme(
    axis.text.x = element_blank(),
    strip.text = element_text(size = 12)
    ) +
  labs(
    x = NULL, 
    y = '\u00B0C', 
    color = 'Water temp.', 
    shape = NULL
  )
toplo6 <- toplo %>% 
  filter(var == 'Salinity (ppt)')
p6 <- ggplot(toplo6, aes(x = yr, y = avev, group = loc, color = loc)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), position = position_dodge2(width = wd), show.legend = F, alpha = 0.7) + 
  geom_point(position = position_dodge2(width = wd), size = 0.5) +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  # geom_line(data = mixmets %>% filter(var == 'Salinity (ppt)'), aes(y = prd)) +
  geom_smooth(method = 'lm', se = F, formula = y ~ x) +
  facet_grid(~ bay_segment) +
  scale_color_manual(values = c('dodgerblue1', 'dodgerblue4')) +
  thm +
  theme(
    strip.text = element_blank()
  ) +
  labs(
    x = NULL, 
    y = 'ppt', 
    color = 'Salinity', 
    shape = NULL
  )

p <- p1 + p2 + p4 + (p5 + p6 + plot_layout(ncol = 1, guides = 'collect')) + plot_layout(ncol = 1, heights = c(1, 1, 1, 2.75)) & theme(legend.position = 'top')

png(here('figs/suppmeteowqraw.png'), height = 6.5, width = 7, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()

# supp gam temp, sali example -----------------------------------------------------------------

fl <- 'OTB_66_tempbot'
load(file = here(paste0('data/', fl, '.RData')))

toplo <- get(fl)

p1 <- show_prdseries(toplo, ylab = 'Bottom temperature (\u00B0C)', col = 'red2') + 
  labs(subtitle ='Station 66, OTB')

fl <- 'OTB_66_salibot'
load(file = here(paste0('data/', fl, '.RData')))

toplo <- get(fl)

p2 <- show_prdseries(toplo, ylab = 'Bottom salinity (ppt)', col = 'dodgerblue2')

p <- p1 + p2 + plot_layout(ncol = 1)

png(here('figs/suppgamex.png'), height = 5, width = 6, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# kendall -------------------------------------------------------------------------------------

leglab <- expression(paste(yr^{-1}))

# kendall all years
sktres <- epcdat %>% 
  select(bay_segment, station = epchc_station, SampleTime, lon = Longitude, lat = Latitude, yr, 
         mo, matches('Bottom')) %>% 
  pivot_longer(matches('Bottom'), names_to = 'var', values_to = 'val') %>% 
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
ktres <- epcdat %>% 
  select(bay_segment, station = epchc_station, SampleTime, lon = Longitude, lat = Latitude, yr, 
         mo, matches('Bottom')) %>% 
  pivot_longer(matches('Bottom'), names_to = 'var', values_to = 'val') %>% 
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
    varsimp = gsub('^(.*?)_.*$', '\\1', var, perl = T),
    varsimp = factor(varsimp, levels = c('Temp', 'Sal'), labels = c('Temperature', 'Salinity')),
  ) %>%
  nest(.by = c('bay_segment', 'mo', 'varsimp', 'var')) %>% 
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
  unnest('sum')

p2 <- ggplot(ktresplo, aes(x = mo, y = bay_segment, fill = nsigper)) + 
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
    y = 'Bay segment',
    title = '(b) % stations with significant trends by month'
  )
    
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
    var = gsub('^(.*?)_.*$', '\\1', var, perl = T),
    var = factor(var, levels = c('Temp', 'Sal'), labels = c('Temperature', 'Salinity')),
  ) %>% 
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

pthm <- theme_bw(base_family = 'serif') +
  theme(
    legend.position = c(0.93, 0.15), 
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

p1 <- ggplot() + 
  annotation_map_tile('cartolight', zoom = 11) +
  geom_sf(data = tbseg, inherit.aes = F) +
  geom_sf(data = toplo, aes(size = abs(slos), fill = slos, shape = coefsgn, color = pvalcol), stroke = 1) +
  facet_grid( ~ var) +
  # scale_fill_gradient2(leglab, low = 'blue', mid = 'grey',  high = 'tomato1', midpoint = 0) +
  scale_fill_gradientn(leglab, limits = colrng, colors = c('blue', 'grey', 'tomato1')) +
  scale_color_manual(values = c(scales::alpha('black', 0), 'black'), guide = 'none', drop = FALSE) +
  # coord_map() +
  scale_x_continuous(breaks = seq(-82.7, -82.4, length = 4)) +
  scale_shape_manual(values = c(24, 25), drop = FALSE, guide = 'none') +
  pthm +
  scale_size(range = c(0.75, 4), guide = 'none') +
  guides(fill = guide_colourbar(barwidth = 0.4, barheight = 2.5)) + 
  labs(
    title = paste0('(a) Change per year, ', yrsel1[1], '-', yrsel1[2]),
    caption = 'Outlines indicate p < 0.05'
  )

p <- p1 + p2 + plot_layout(ncol = 2, width = c(0.9, 1))

png(here('figs/kendall.png'), height = 4.5, width = 9.5, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# supp kendall --------------------------------------------------------------------------------

leglab <- expression(paste(yr^{-1}))

# kendall all years
sktres <- epcdat %>% 
  select(bay_segment, station = epchc_station, SampleTime, lon = Longitude, lat = Latitude, yr, 
         mo, matches('Bottom')) %>% 
  filter(yr >= yrsel2[1] & yr <= yrsel2[2]) %>% 
  pivot_longer(matches('Bottom'), names_to = 'var', values_to = 'val') %>% 
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
ktres <- epcdat %>% 
  select(bay_segment, station = epchc_station, SampleTime, lon = Longitude, lat = Latitude, yr, 
         mo, matches('Bottom')) %>% 
  filter(yr >= yrsel2[1] & yr <= yrsel2[2]) %>% 
  pivot_longer(matches('Bottom'), names_to = 'var', values_to = 'val') %>% 
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
    varsimp = gsub('^(.*?)_.*$', '\\1', var, perl = T),
    varsimp = factor(varsimp, levels = c('Temp', 'Sal'), labels = c('Temperature', 'Salinity'))
  ) %>%
  nest(.by = c('bay_segment', 'mo', 'varsimp', 'var')) %>% 
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
  unnest('sum')

p2 <- ggplot(ktresplo, aes(x = mo, y = bay_segment, fill = nsigper)) + 
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
    y = 'Bay segment',
    title = '(b) % stations with significant trends by month'
  )

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
    var = gsub('^(.*?)_.*$', '\\1', var, perl = T),
    var = factor(var, levels = c('Temp', 'Sal'), labels = c('Temperature', 'Salinity'))
  ) %>% 
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

pthm <- theme_bw(base_family = 'serif') +
  theme(
    legend.position = c(0.93, 0.15), 
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

p1 <- ggplot() +
  annotation_map_tile('cartolight', zoom = 11) +
  geom_sf(data = tbseg, inherit.aes = F) +
  geom_sf(data = toplo, aes(fill = slos, shape = coefsgn, color = pvalcol), stroke = 1) +
  facet_grid( ~ var) +
  # scale_fill_gradient2(leglab, low = 'blue', mid = 'grey',  high = 'tomato1', midpoint = 0) +
  scale_fill_gradientn(leglab, limits = colrng, colors = c('blue', 'grey', 'tomato1')) +
  scale_color_manual(values = c(scales::alpha('black', 0), 'black'), guide = 'none', drop = FALSE) +
  # coord_map() +
  scale_x_continuous(breaks = seq(-82.7, -82.4, length = 4)) +
  scale_shape_manual(values = c(24, 25), drop = FALSE, guide = 'none') +
  pthm +
  scale_size(range = c(0.75, 4), guide = 'none') +
  guides(fill = guide_colourbar(barwidth = 0.4, barheight = 2.5)) + 
  labs(
    title = paste0('(a) Change per year, ', yrsel2[1], '-', yrsel2[2]),
    caption = 'Outlines indicate p < 0.05'
  )

p <- p1 + p2 + plot_layout(ncol = 2, width = c(0.9, 1))

png(here('figs/suppkendall.png'), height = 4.5, width = 9.5, family = 'serif', units = 'in', res = 300)
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
    se = sd(value, na.rm = T)^2,
    lov = t.test(value)$conf.int[1], 
    hiv = t.test(value)$conf.int[2], 
    cnt = n(),
    .by = c(bay_segment, yr, name)
  ) %>% 
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')),
    name = factor(name, levels = c('temp', 'sal'), labels = c('Water temp. (\u00B0C)', 'Salinity (ppt)'))
  )

# # get mixed-effect meta-analysis linear preds
# mixmets <- toplo %>% 
#   group_nest(bay_segment, name) %>% 
#   mutate(
#     data = purrr::map(data, function(x){
#       
#       mod <- mixmeta(avev ~ yr, S = se, random = ~1|yr, data = x, method = 'reml')
#       
#       out <- x %>% 
#         mutate(
#           prd = predict(mod)
#         ) %>% 
#         select(yr, prd)
#       
#       return(out)
#       
#     }) 
#   ) %>% 
#   unnest(data)

thm <- theme_minimal() + 
  theme(
    strip.placement = 'outside', 
    panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 11), 
    legend.text = element_text(size= 12), 
    axis.text.x = element_text(size = 7), 
    legend.position = 'none'
  )

p <- ggplot(toplo, aes(x = yr, y = avev, color = name)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), show.legend = F, alpha = 0.7) + 
  geom_point(size = 0.75) +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  # geom_line(data = mixmets, aes(y = prd), color = 'steelblue4') +
  scale_color_manual(values = c('red2', 'dodgerblue2')) +
  geom_smooth(method = 'lm', se = F, formula = y ~ x) +
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

# supp pinco temp, sal trends -----------------------------------------------------------------

load(file = here('data/pincotemp.RData'))

toplo <- pincotemp %>% 
  pivot_longer(temp:sal, names_to = 'var', values_to = 'val') %>% 
  summarise(
    avev = mean(val, na.rm = T), 
    se = sd(val, na.rm = T)^2,
    lov = t.test(val)$conf.int[1], 
    hiv = t.test(val)$conf.int[2], 
    cnt = n(),
    nmo = length(unique(mo)),
    .by = c(yr, var)
  ) %>% 
  mutate(
    var = factor(var, levels = c('temp', 'sal'), labels = c('Water temp. (\u00B0C)', 'Salinity (ppt)'))
  )

# # get mixed-effect meta-analysis linear preds
# mixmets <- toplo %>% 
#   group_nest(var) %>% 
#   mutate(
#     data = purrr::map(data, function(x){
#       
#       mod <- mixmeta(avev ~ yr, S = se, random = ~1|yr, data = x, method = 'reml')
#       
#       out <- x %>% 
#         mutate(
#           prd = predict(mod)
#         ) %>% 
#         select(yr, prd)
#       
#       return(out)
#       
#     }) 
#   ) %>% 
#   unnest(data)

thm <- theme_minimal() + 
  theme(
    strip.placement = 'outside', 
    panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 10), 
    legend.text = element_text(size= 12), 
    axis.text.x = element_text(size = 10), 
    legend.position = 'none'
  )

# MTB sampling gap? month coverage? change in sample design?
p <- ggplot(toplo, aes(x = yr, y = avev, color = var)) + 
  geom_linerange(aes(ymin = lov, ymax = hiv), show.legend = F, alpha = 0.7) + 
  geom_point(size = 0.75) +
  # scale_x_continuous(breaks = seq(min(toplo$yr), max(toplo$yr), by = 3)) +
  # geom_line(data = mixmets, aes(y = prd), color = 'steelblue4', method = 'lm', se =) +
  scale_color_manual(values = c('red2', 'dodgerblue2')) +
  geom_smooth(formula = y ~ x, method = 'lm', se = F) +
  facet_grid(var ~ ., switch = 'y', scales = 'free_y') +
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

png(here('figs/supppincotrnds.png'), height = 3.5, width = 3.5, family = 'serif', units = 'in', res = 400)
print(p)
dev.off()

# mixeff example plot -------------------------------------------------------------------------

load(file = here('data/mixmodprds.RData'))

tmpthr <- '30'

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
  geom_text(data = subset(toplo3, thrtyp == 'Temperature > 30'), aes(x = 1975, y = -10, label = slo, col = thrtyp), show.legend = F, 
             vjust = 0, hjust = 0, fontface = 'italic', size = 3) +
  geom_text(data = subset(toplo3, thrtyp == 'Salinity < 25'), aes(x = 1975, y = -30, label = slo, col = thrtyp), show.legend = F, 
            vjust = 0, hjust = 0, fontface = 'italic', size = 3) +
  geom_text(data = subset(toplo3, thrtyp == 'Both'), aes(x = 1975, y = -12, label = slo, col = thrtyp), show.legend = F, 
            vjust = 0, hjust = 0, fontface = 'italic', size = 3) +
  
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

tmpthr <- '30'

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
  geom_text(data = subset(toplo3, thrtyp == 'Temperature > 30'), aes(x = 1998, y = -10, label = slo, col = thrtyp), show.legend = F, 
            vjust = 0, hjust = 0, fontface = 'italic', size = 3) +
  geom_text(data = subset(toplo3, thrtyp == 'Salinity < 25'), aes(x = 1998, y = -30, label = slo, col = thrtyp), show.legend = F, 
            vjust = 0, hjust = 0, fontface = 'italic', size = 3) +
  geom_text(data = subset(toplo3, thrtyp == 'Both'), aes(x = 1998, y = -10, label = slo, col = thrtyp), show.legend = F, 
            vjust = 0, hjust = 0, fontface = 'italic', size = 3) +
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

# seagrass change plots -----------------------------------------------------------------------

load(file = here('data/sgmods.RData'))

epcmod1 <- sgmods$epcmod1
fimmod <- sgmods$fimmod
pincomod <- sgmods$pincomod

##
# epc model plot

epcplo <- modplo_fun(
  mod = epcmod1, 
  xlab1 = '# days temperature > 30 \u00B0C', 
  ylab1 = 'Frequency occurrence', 
  xlab2 = '# days salinity < 25 ppt',
  title1 = '(a) EPC', 
  ismetric = T
)

## 
# fim model plot

fimplo <- modplo_fun(
  mod = fimmod, 
  xlab1 = 'Temperature (\u00B0C)', 
  ylab1 = '% cover', 
  xlab2 = 'Salinity (ppt)',
  title1 = '(b) FIM'
)

# ## 
# # pinco model plot, only yrcat significant
# 
# pincoplo <- modplo_fun(
#   mod = pincomod, 
#   xlab1 = 'Temperature (\u00B0C)', 
#   ylab1 = 'Probability of occurrence', 
#   xlab2 = 'Salinity (ppt)', 
#   isbinom = T,
#   title1 = '(c) PDEM'
# )

p <- epcplo / fimplo 

png(here('figs/sgmods.png'), height = 6, width = 7, family = 'serif', units = 'in', res = 300)
print(p)
dev.off()

# supp hillsborough flow v temp over time -----------------------------------------------------

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

png(here('figs/suppflowtemp.png'), height = 4, width = 6, family = 'serif', units = 'in', res = 600)
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
