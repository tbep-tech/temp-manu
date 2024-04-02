library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)

data(epccmbdat)

tomod <- epccmbdat %>% 
  filter(!bay_segment %in% c('LTB')) #%>% 
  # .[-22, ]

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

mod <- gam(total ~ te(la, yr, by = bay_segment) + te(temp, yr, by = bay_segment) + te(sal, yr, by = bay_segment), data = tomod, method = 'REML')

salslc <- data_slice(mod,
                    yr = c(2000, 2020),
                    sal = evenly(sal, n = 100), 
                    bay_segment = evenly(bay_segment)
                    )
toploslc <- smooth_estimates(mod, select = 'te(sal,yr)', data = salslc, partial_match = T) %>% 
  add_confint() %>% 
  filter(!((value < min(tomod[tomod$bay_segment == 'HB', 'sal']) | value > max(tomod[tomod$bay_segment == 'HB', 'sal'])) & bay_segment == 'HB' & var == 'sal')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'OTB', 'sal']) | value > max(tomod[tomod$bay_segment == 'OTB', 'sal'])) & bay_segment == 'OTB' & var == 'sal')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'MTB', 'sal']) | value > max(tomod[tomod$bay_segment == 'MTB', 'sal'])) & bay_segment == 'MTB' & var == 'sal')) 

p1 <- ggplot(toploslc, aes(x = sal, y = .estimate, color = factor(yr), fill = factor(yr), group = yr)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  facet_wrap(~bay_segment, scales = 'free')
tempslc <- data_slice(mod,
                    yr = c(2000, 2020),
                    temp = evenly(temp, n = 100), 
                    bay_segment = evenly(bay_segment)
                    )
toploslc <- smooth_estimates(mod, select = 'te(temp,yr)', data = tempslc, partial_match = T) %>% 
  add_confint() %>% 
  filter(!((.estimate < min(tomod[tomod$bay_segment == 'HB', 'temp']) | .estimate > max(tomod[tomod$bay_segment == 'HB', 'temp'])) & bay_segment == 'HB' & var == 'temp')) %>%
  filter(!((.estimate < min(tomod[tomod$bay_segment == 'OTB', 'temp']) | .estimate > max(tomod[tomod$bay_segment == 'OTB', 'temp'])) & bay_segment == 'OTB' & var == 'temp')) %>%
  filter(!((.estimate < min(tomod[tomod$bay_segment == 'MTB', 'temp']) | .estimate > max(tomod[tomod$bay_segment == 'MTB', 'temp'])) & bay_segment == 'MTB' & var == 'temp'))

p2 <- ggplot(toploslc, aes(x = temp, y = .estimate, color = factor(yr), fill = factor(yr), group = yr)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  facet_wrap(~bay_segment, scales = 'free')
laslc <- data_slice(mod,
                    yr = c(2000, 2020),
                    la = evenly(la, n = 100), 
                    bay_segment = evenly(bay_segment)
                    )
toploslc <- smooth_estimates(mod, select = 'te(la,yr)', data = laslc, partial_match = T) %>% 
  add_confint() %>% 
  filter(!((value < min(tomod[tomod$bay_segment == 'HB', 'la']) | value > max(tomod[tomod$bay_segment == 'HB', 'la'])) & bay_segment == 'HB' & var == 'la')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'OTB', 'la']) | value > max(tomod[tomod$bay_segment == 'OTB', 'la'])) & bay_segment == 'OTB' & var == 'la')) %>%
  filter(!((value < min(tomod[tomod$bay_segment == 'MTB', 'la']) | value > max(tomod[tomod$bay_segment == 'MTB', 'la'])) & bay_segment == 'MTB' & var == 'la'))

p3 <- ggplot(toploslc, aes(x = la, y = .estimate, color = factor(yr), fill = factor(yr), group = yr)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  facet_wrap(~bay_segment, scales = 'free')


p1 + p2 + p3 + plot_layout(ncol = 1, guides = 'collect')


toplo <- smooth_estimates(mod, smooth = data = toprd) %>% 
  add_confint() %>% 
  filter(.type == 'TPRS') %>% 
  select(.smooth, .estimate, yr, la, temp, sal, .lower_ci, .upper_ci) %>% 
  pivot_longer(cols = c(yr, la, temp, sal), names_to = 'var', values_to = 'value') %>% 
  filter(!is.na(value)) %>% 
  mutate(
    bay_segment = gsub('^.*segment', '', .smooth)
  )

ggplot(toplo, aes(x = sal, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = yr, color = yr) +
  geom_line(aes(color = yr))+
  facet_wrap(bay_segment~yr, scales = 'free')
