library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)

# EPC eval ------------------------------------------------------------------------------------

data(epccmbdat)

tomod <- epccmbdat %>% 
  filter(!bay_segment %in% c('LTB'))

##
# single smoothers only, temp, sal

mod <- gam(total ~ s(yr, by = bay_segment) + s(la, by = bay_segment) + s(temp, by = bay_segment) + s(sal, by = bay_segment), data = tomod, method = 'REML')

toplo <- smooth_estimates(mod) %>% 
  add_confint() %>% 
  # add_residuals(data = ., model = mod)
  select(.smooth, .estimate, yr, la, temp, sal, .lower_ci, .upper_ci) %>% 
  pivot_longer(cols = c(yr, la, temp, sal), names_to = 'var', values_to = 'value') %>% 
  filter(!is.na(value)) %>% 
  mutate(
    bay_segment = gsub('^.*segment', '', .smooth)
  ) %>% 
  filter(!((value < min(tomod[tomod$bay_segment == 'HB', 'sal']) | value > max(tomod[tomod$bay_segment == 'HB', 'sal'])) & bay_segment == 'HB' & var == 'sal')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'HB', 'temp']) | value > max(tomod[tomod$bay_segment == 'HB', 'temp'])) & bay_segment == 'HB' & var == 'temp')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'HB', 'la']) | value > max(tomod[tomod$bay_segment == 'HB', 'la'])) & bay_segment == 'HB' & var == 'la')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'OTB', 'sal']) | value > max(tomod[tomod$bay_segment == 'OTB', 'sal'])) & bay_segment == 'OTB' & var == 'sal')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'OTB', 'temp']) | value > max(tomod[tomod$bay_segment == 'OTB', 'temp'])) & bay_segment == 'OTB' & var == 'temp')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'OTB', 'la']) | value > max(tomod[tomod$bay_segment == 'OTB', 'la'])) & bay_segment == 'OTB' & var == 'la')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'MTB', 'sal']) | value > max(tomod[tomod$bay_segment == 'MTB', 'sal'])) & bay_segment == 'MTB' & var == 'sal')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'MTB', 'temp']) | value > max(tomod[tomod$bay_segment == 'MTB', 'temp'])) & bay_segment == 'MTB' & var == 'temp')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'MTB', 'la']) | value > max(tomod[tomod$bay_segment == 'MTB', 'la'])) & bay_segment == 'MTB' & var == 'la'))


ggplot(toplo, aes(x = value, y = .estimate, color = bay_segment, fill = bay_segment)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  facet_wrap(~var, scales = 'free')

##
# single smoothers only, both

mod <- gam(total ~ s(yr, by = bay_segment) + s(la, by = bay_segment) + s(both, by = bay_segment), data = tomod, method = 'REML')

toplo <- data_slice(mod, 
                    both = evenly(both), 
                    la = evenly(la),
                    bay_segment = c('HB', 'OTB', 'MTB'), 
                    yr = c(2000, 2020)
) %>% 
  fitted_values(mod, data = ., scale = 'response') %>% 
  mutate(
    .fitted = ifelse(.fitted < 0 | .fitted > 1, NA, .fitted)
  )

ggplot(toplo, aes(x = both, y = .fitted, color = factor(yr), fill = factor(yr), group = yr)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line() + 
  facet_wrap(~bay_segment)

ggplot(toplo, aes(x = la, y = both, z = .fitted)) + 
  geom_tile(aes(fill = .fitted)) + 
  facet_wrap(yr~bay_segment) + 
  scale_fill_viridis_c() + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))

# FIM eval------------------------------------------------------------------------------------

data(fimsgtempdat)

tomod <- fimsgtempdat %>% 
  filter(!bay_segment %in% c('LTB')) %>% 
  select(date, sgcov, sgpres, yr, temp, sal, secchi_m, secchi_on_bottom, bay_segment) %>% 
  na.omit() %>%
  mutate(mo = month(date)) %>% 
  # filter(month(date) %in% c(6:9)) %>%
  summarise(
    sgcov = mean(sgcov),
    sgpres = mean(sgpres),
    temp = mean(temp),
    sal = mean(sal),
    secchi_m = mean(secchi_m),
    .by = c(bay_segment, yr, mo)
  ) %>% 
  mutate(
    la = dplyr::case_when(
      bay_segment %in% "OTB" ~ 1.49 / secchi_m,
      bay_segment %in% "HB" ~ 1.61 / secchi_m,
      bay_segment %in% "MTB" ~ 1.49 / secchi_m,
      bay_segment %in% "LTB" ~ 1.84 / secchi_m
    )
  )
  # filter(!secchi_on_bottom)

##
# single smoothers only

mod <- gam(sgcov ~ s(yr, by = bay_segment) + s(la, by = bay_segment) + s(temp, by = bay_segment) + s(sal, by = bay_segment) + s(mo, bs = 'cc') + ti(yr, mo, bs = c('tp', 'cc'), by = bay_segment) + ti(la, mo, bs = c('tp', 'cc'), by = bay_segment) + ti(temp, mo, bs = c('tp', 'cc'), by = bay_segment) + ti(sal, mo, bs = c('tp', 'cc'), by = bay_segment), knots = list(doy = c(1, 12)), data = tomod, method = 'REML')

toplo <- smooth_estimates(mod) %>% 
  add_confint() %>% 
  # add_residuals(data = ., model = mod)
  select(.smooth, .estimate, yr, secchi_m, temp, sal, .lower_ci, .upper_ci) %>% 
  pivot_longer(cols = c(yr, secchi_m, temp, sal), names_to = 'var', values_to = 'value') %>% 
  filter(!is.na(value)) %>% 
  mutate(
    bay_segment = gsub('^.*segment', '', .smooth)
  ) %>% 
  filter(!((value < min(tomod[tomod$bay_segment == 'HB', 'sal']) | value > max(tomod[tomod$bay_segment == 'HB', 'sal'])) & bay_segment == 'HB' & var == 'sal')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'HB', 'temp']) | value > max(tomod[tomod$bay_segment == 'HB', 'temp'])) & bay_segment == 'HB' & var == 'temp')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'HB', 'secchi_m']) | value > max(tomod[tomod$bay_segment == 'HB', 'secchi_m'])) & bay_segment == 'HB' & var == 'secchi_m')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'OTB', 'sal']) | value > max(tomod[tomod$bay_segment == 'OTB', 'sal'])) & bay_segment == 'OTB' & var == 'sal')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'OTB', 'temp']) | value > max(tomod[tomod$bay_segment == 'OTB', 'temp'])) & bay_segment == 'OTB' & var == 'temp')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'OTB', 'secchi_m']) | value > max(tomod[tomod$bay_segment == 'OTB', 'secchi_m'])) & bay_segment == 'OTB' & var == 'secchi_m')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'MTB', 'sal']) | value > max(tomod[tomod$bay_segment == 'MTB', 'sal'])) & bay_segment == 'MTB' & var == 'sal')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'MTB', 'temp']) | value > max(tomod[tomod$bay_segment == 'MTB', 'temp'])) & bay_segment == 'MTB' & var == 'temp')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'MTB', 'secchi_m']) | value > max(tomod[tomod$bay_segment == 'MTB', 'secchi_m'])) & bay_segment == 'MTB' & var == 'secchi_m'))


ggplot(toplo, aes(x = value, y = .estimate, color = bay_segment, fill = bay_segment)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  facet_wrap(bay_segment~var, scales = 'free')

# PINCO eval----------------------------------------------------------------------------------

data("pincotemp")

tomod <- pincotemp %>% 
  filter(!bay_segment %in% c('LTB')) %>% 
  select(date, allsg, yr, site, temp, sal, secchi_m, secchi_on_bottom, bay_segment) %>% 
  na.omit() %>%
  mutate(
    mo = month(date),
    site = factor(site)
    ) %>% 
  summarise(
    allsg = sum(allsg) / length(allsg),
    temp = mean(temp),
    sal = mean(sal),
    secchi_m = mean(secchi_m),
    .by = c(yr, mo, site)
  ) %>% 
  mutate(
    la = 1.49 / secchi_m
  )
# filter(!secchi_on_bottom)

##
# single smoothers only

mod <- gam(allsg ~ s(yr, by = site) + s(la, by = site) + s(temp, by = site) + s(sal, by = site) + s(mo, bs = 'cc') + ti(yr, mo, bs = c('tp', 'cc'), by = site) + ti(la, mo, bs = c('tp', 'cc'), by = site) + ti(temp, mo, bs = c('tp', 'cc'), by = site) + ti(sal, mo, bs = c('tp', 'cc'), by = site), knots = list(doy = c(1, 12)), data = tomod, method = 'REML')

toplo <- smooth_estimates(mod) %>% 
  add_confint() %>% 
  # add_residuals(data = ., model = mod)
  select(.smooth, .estimate, yr, secchi_m, temp, sal, .lower_ci, .upper_ci) %>% 
  pivot_longer(cols = c(yr, secchi_m, temp, sal), names_to = 'var', values_to = 'value') %>% 
  filter(!is.na(value))

ggplot(toplo, aes(x = value, y = .estimate)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  facet_wrap(~var, scales = 'free')

