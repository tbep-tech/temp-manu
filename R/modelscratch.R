library(tidyverse)
library(lmerTest)
library(tbeptools)
library(here)
library(vegan)
library(scales)

load(file = here('data/thrdat.RData'))

seglng <- c('Old Tampa Bay', 'Hillsborough Bay', 'Middle Tampa Bay', 'Lower Tampa Bay')
segshr <- c('OTB', 'HB', 'MTB', 'LTB')

##
# seagrass fo
transectocc <- anlz_transectocc(transect)
transectavespp <- anlz_transectavespp(transectocc, by_seg = TRUE)

fodat <- transectavespp %>% 
  select(-nsites) %>% 
  filter(bay_segment %in% segshr) %>% 
  filter(Savspecies %in% c('Halodule', 'Syringodium', 'Thalassia')) %>% 
  mutate(
    Savspecies = droplevels(Savspecies)
  ) %>% 
  pivot_wider(names_from = 'Savspecies', values_from = 'foest', values_fill = 0)

##
# threshold count data

# these had strong models
salthr <- '25'
tmpthr <- '30'

cntdat <- thrdat %>% 
  filter(salithr == paste0('sali_', salthr)) %>% 
  filter(tempthr == paste0('temp_', tmpthr)) %>% 
  summarise(
    cnt = mean(cnt),
    .by = c('thrtyp', 'bay_segment', 'yr')
  ) %>% 
  filter(yr >= min(fodat$yr) & yr <= max(fodat$yr)) %>% 
  pivot_wider(names_from = thrtyp, values_from = cnt)

##
# load data

load(file = url('https://github.com/tbep-tech/load-estimates/raw/main/data/totanndat.RData'))

tndat <- totanndat %>% 
  select(bay_segment, yr = year, tn_load) %>% 
  filter(yr >= min(fodat$yr) & yr <= max(fodat$yr)) %>% 
  filter(bay_segment %in% seglng) %>% 
  mutate(
    bay_segment = factor(bay_segment, levels = seglng, labels = segshr)
  )

##
# chl data

chldat <- epcdata %>% 
  select(bay_segment, yr, chla) %>% 
  summarise(
    chla = mean(chla, na.rm = T),
    .by = c('bay_segment', 'yr')
  ) %>% 
  filter(yr >= min(fodat$yr) & yr <= max(fodat$yr)) %>% 
  filter(bay_segment %in% segshr) %>% 
  mutate(
    bay_segment = factor(bay_segment, levels = segshr)
  )

##
# combine

# NMS
alldat <- fodat %>% 
  inner_join(cntdat, by = c('yr', 'bay_segment')) %>% 
  inner_join(tndat, by = c('yr', 'bay_segment')) %>% 
  inner_join(chldat, by = c('yr', 'bay_segment'))
  # filter(bay_segment == 'HB')

toord <- alldat %>% 
  select(-yr, -bay_segment) %>% 
  mutate_all(rescale, to = c(0, 1))

vec_ext <- 1
coord_fix <- F
size <- 3
repel <- F
arrow <- 0.2
txt <- 3
alpha <- 0.8
ext <- 1.1
exp <- 0.1
parse <- T
ellipse <- F

orddat <- metaMDS(toord, k = 3, trymax = 50)
grps <- alldat$yr

ggord(orddat, axes = c('1', '2'), ellipse = ellipse, grp_in = grps,
      parse = parse, vec_ext = vec_ext, coord_fix = coord_fix, size = size, 
      repel = repel, arrow = arrow, txt = txt, alpha = alpha, ext = ext, 
      exp = exp) +
  theme(
    legend.title = element_blank(), 
    legend.position = 'top'
  )


alldat <- fodat %>% 
  filter(!(Halodule == 0 & Syringodium == 0 & Thalassia == 0)) %>%
  inner_join(cntdat, by = c('yr', 'bay_segment')) %>% 
  inner_join(tndat, by = c('yr', 'bay_segment'))# %>% 
# filter(bay_segment == 'HB')

# RDA
toord <- tmp %>% 
  select(-yr, -bay_segment)

fotmp <- toord[, c('Halodule', 'Syringodium', 'Thalassia')]
envtmp <- toord[, c('salicnt', 'tempcnt', 'bothcnt', 'tn_load')] %>% 
  mutate_all(rescale, to = c(0, 1))
ord <- rda(fotmp, envtmp, scale = F)

vec_ext <- 0.7
coord_fix <- T
addsize <- 4
size <- 3
repel <- F
arrow <- 0.2
txt <- 4
alpha <- 0.8
ext <- 1.2
exp <- 0.07
parse <- T
ellipse <- F

grps <- alldat$bay_segment

ggord(ord, axes = c('1', '2'),  ptslab = T, ellipse = ellipse, grp_in = grps, addsize = addsize,
      parse = parse, vec_ext = vec_ext, coord_fix = coord_fix, size = size, 
      repel = repel, arrow = arrow, txt = txt, alpha = alpha, ext = ext, 
      exp = exp)


tomod <- alldat %>% 
  filter(bay_segment == 'OTB')
mod <- gam(Thalassia ~ s(tn_load) + s(salicnt), data = tomod)

mod <- lm(Thalassia ~ tn_load + salicnt + tempcnt + bothcnt, data = alldat)
mod1 <- lm(Thalassia ~ tn_load, data = alldat)

mods <- alldat %>% 
  pivot_longer(names_to = 'spp', values_to = 'fo', Halodule:Thalassia) %>% 
  group_by(bay_segment, spp) %>% 
  nest() %>% 
  mutate(
    mod = purrr::map(data, function(x){

      mod <- lm(fo ~ chla, data = x)
      res <- residuals(mod)
      x$resid <- res
      mod <- lm(resid ~ 0 + bothcnt, data = x)
      # mod <- step(mod)
      return(mod)
      
    })
  )
  