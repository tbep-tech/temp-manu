library(tidyverse)
library(lubridate)
library(EnvStats)
library(tbeptools)


# kendall -------------------------------------------------------------------------------------

sktres <- epcdata %>% 
  select(bay_segment, station = epchc_station, SampleTime, lon = Longitude, lat = Latitude, yr, 
         mo, matches('Temp')) %>% 
  pivot_longer(matches('Temp'), names_to = 'var', values_to = 'val') %>% 
  filter(yr < 2023) %>% 
  nest(.by = c('bay_segment', 'var', 'station', 'lon', 'lat')) %>% 
  mutate(
    skt = purrr::pmap(list(station, var, data), function(station, var, data){
      
      cat(station, var, '\n')

      tomod <- data %>% 
        arrange(yr, mo) %>% 
        na.omit()
      
      # run sk test
      ests_sk <- kendallSeasonalTrendTest(val ~ mo + yr, data = tomod)
      
      out <- tibble(
        pval = ests_sk$p.value[2], 
        slos = ests_sk$estimate[2],
        n = nrow(tomod)
      )
      
      return(out)
      
    })
  ) %>% 
  select(-data) %>% 
  unnest(skt)

# # basemap
# dat_ext <- unname(st_bbox(tbseg))
# bsmap1 <- get_stamenmap(bbox = dat_ext, maptype = 'terrain-background', zoom = 11)
# 
# # change opacity of basemap
# mapatt <- attributes(bsmap1)
# bsmap1_transparent <- matrix(adjustcolor(bsmap1, 
#                                          alpha.f = 0.4), 
#                              nrow = nrow(bsmap1))
# attributes(bsmap1_transparent) <- mapatt
# 
# chlleglab <- expression(paste(log[10], " chl-a change (", italic(mu), "g ", L^-1, yr^-1, ")"))
# tnleglab <- expression(paste(log[10], " TN change (mg ", L^-1, yr^-1, ")"))
# 
# pthm <- theme_bw(base_family = 'serif', base_size = 14) +
#   theme(
#     legend.position = 'top',
#     legend.box = 'vertical', 
#     strip.background = element_blank(),
#     axis.title = element_blank(), 
#     axis.text = element_text(size = 9),
#     strip.text = element_text(size = 14)
#   )
# 
# toplo <- modcmp %>%
#   mutate(
#     pvalcol = ifelse(pval < 0.05, T, F), 
#     coefsgn = sign(slos), 
#     coefsgn = factor(coefsgn, levels = c('1', '-1'), labels = c('inc', 'dec')), 
#     test = factor(test, levels = c('mt', 'sk'), labels = c('GAM/Mixed', 'Seasonal Kendall'))
#   ) %>% 
#   left_join(stations, by = c('segment' = 'bay_segment', 'station' = 'epchc_station'))
# 
# toplo1 <- toplo %>% 
#   filter(param == 'tn')
# 
# toplo2 <- toplo %>% 
#   filter(param == 'chla')
# 
# p1 <- ggmap(bsmap1_transparent) +
#   geom_point(data = toplo1, aes(x = Longitude, y = Latitude, size = abs(slos), fill = slos, shape = coefsgn, color = pvalcol), stroke = 1) +
#   facet_grid(~ test) + 
#   scale_fill_gradient2(tnleglab, low = 'green', mid = 'grey',  high = 'tomato1', midpoint = 0) +
#   scale_color_manual(values = c(scales::alpha('black', 0), 'black'), guide = 'none', drop = FALSE) +
#   coord_map() + 
#   scale_shape_manual(values = c(24, 25), guide = 'none', drop = FALSE) +
#   pthm +
#   scale_size(range = c(1, 6), guide = F) +
#   guides(fill = guide_colourbar(barheight = 0.7, barwidth = 16, title.position = 'top', title.hjust = 0.5)) 
# 
# p2 <- ggmap(bsmap1_transparent) +
#   geom_point(data = toplo2, aes(x = Longitude, y = Latitude, size = abs(slos), fill = slos, shape = coefsgn, color = pvalcol), stroke = 1) +
#   facet_grid(~ test) + 
#   scale_fill_gradient2(chlleglab, low = 'green', mid = 'grey',  high = 'tomato1', midpoint = 0) +
#   scale_color_manual(values = c(scales::alpha('black', 0), 'black'), guide = 'none', drop = FALSE) +
#   coord_map() + 
#   scale_shape_manual(values = c(24, 25), guide = 'none', drop = FALSE) +
#   pthm +
#   scale_size(range = c(1, 6), guide = F) +
#   guides(fill = guide_colourbar(barheight = 0.7, barwidth = 16, title.position = 'top', title.hjust = 0.5)) + 
#   labs(
#     caption = 'Outlines indicate p < 0.05'
#   )
# 
# p <- p1 + p2 + plot_layout(ncol = 2)
# png('figs/trndcmp.png', height = 6, width = 11, family = 'serif', units = 'in', res = 300)
# p
# dev.off()
# 
