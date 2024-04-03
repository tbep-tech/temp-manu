library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)

# EPC eval ------------------------------------------------------------------------------------

data(epccmbdat)

tomod <- epccmbdat %>% 
  filter(!bay_segment %in% c('LTB'))

##
# single smoothers only

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
# interaction terms

mod <- gam(total ~ te(la, yr, by = bay_segment) + te(temp, yr, by = bay_segment) + te(sal, yr, by = bay_segment), data = tomod, method = 'REML')

salslc <- data_slice(mod,
                    yr = c(2000, 2020),
                    sal = evenly(sal, n = 100), 
                    bay_segment = evenly(bay_segment)
                    )
toploslc <- smooth_estimates(mod, select = 'te(sal,yr)', data = salslc, partial_match = T) %>% 
  add_confint() %>% 
  filter(!((sal < min(tomod[tomod$bay_segment == 'HB', 'sal']) | sal > max(tomod[tomod$bay_segment == 'HB', 'sal'])) & bay_segment == 'HB')) %>%
  filter(!((sal < min(tomod[tomod$bay_segment == 'OTB', 'sal']) | sal > max(tomod[tomod$bay_segment == 'OTB', 'sal'])) & bay_segment == 'OTB')) %>%
  filter(!((sal < min(tomod[tomod$bay_segment == 'MTB', 'sal']) | sal > max(tomod[tomod$bay_segment == 'MTB', 'sal'])) & bay_segment == 'MTB')) 

p1 <- ggplot(toploslc, aes(x = sal, y = .estimate, color = factor(yr), fill = factor(yr), group = yr)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  facet_wrap(~bay_segment) +#, scales = 'free_y') +
  labs(
    x = NULL, 
    y = 'Sal partial'
  )
tempslc <- data_slice(mod,
                    yr = c(2000, 2020),
                    temp = evenly(temp, n = 100), 
                    bay_segment = evenly(bay_segment)
                    )
toploslc <- smooth_estimates(mod, select = 'te(temp,yr)', data = tempslc, partial_match = T) %>% 
  add_confint() %>% 
  filter(!((temp < min(tomod[tomod$bay_segment == 'HB', 'temp']) | temp > max(tomod[tomod$bay_segment == 'HB', 'temp'])) & bay_segment == 'HB')) %>%
  filter(!((temp < min(tomod[tomod$bay_segment == 'OTB', 'temp']) | temp > max(tomod[tomod$bay_segment == 'OTB', 'temp'])) & bay_segment == 'OTB')) %>%
  filter(!((temp < min(tomod[tomod$bay_segment == 'MTB', 'temp']) | temp > max(tomod[tomod$bay_segment == 'MTB', 'temp'])) & bay_segment == 'MTB'))

p2 <- ggplot(toploslc, aes(x = temp, y = .estimate, color = factor(yr), fill = factor(yr), group = yr)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  facet_wrap(~bay_segment) +#, scales = 'free_y') + 
  labs(
    x = NULL, 
    y = 'Temp partial'
  )
laslc <- data_slice(mod,
                    yr = c(2000, 2020),
                    la = evenly(la, n = 100), 
                    bay_segment = evenly(bay_segment)
                    )
toploslc <- smooth_estimates(mod, select = 'te(la,yr)', data = laslc, partial_match = T) %>% 
  add_confint() %>% 
  filter(!((la < min(tomod[tomod$bay_segment == 'HB', 'la']) | la > max(tomod[tomod$bay_segment == 'HB', 'la'])) & bay_segment == 'HB')) %>%
  filter(!((la < min(tomod[tomod$bay_segment == 'OTB', 'la']) | la > max(tomod[tomod$bay_segment == 'OTB', 'la'])) & bay_segment == 'OTB')) %>%
  filter(!((la < min(tomod[tomod$bay_segment == 'MTB', 'la']) | la > max(tomod[tomod$bay_segment == 'MTB', 'la'])) & bay_segment == 'MTB'))

p3 <- ggplot(toploslc, aes(x = la, y = .estimate, color = factor(yr), fill = factor(yr), group = yr)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  facet_wrap(~bay_segment) +#, scales = 'free_y') + 
  labs(
    x = NULL, 
    y = 'LA partial'
  )

p1 + p2 + p3 + plot_layout(ncol = 1, guides = 'collect')

# FIM eval------------------------------------------------------------------------------------

data(fimsgtempdat)

tomod <- fimsgtempdat %>% 
  filter(!bay_segment %in% c('LTB')) %>% 
  select(date, sgcov, sgpres, yr, la, temp, sal, secchi_m, secchi_on_bottom, bay_segment) %>% 
  na.omit() %>%
  # filter(month(date) %in% c(7:10)) %>% 
  summarise(
    sgpres = mean(sgpres),
    la = mean(la),
    temp = mean(temp),
    sal = mean(sal),
    secchi_m = mean(secchi_m),
    .by = c(bay_segment, yr)
  )
  # filter(!secchi_on_bottom)

##
# single smoothers only

mod <- gam(sgpres ~ s(yr, by = bay_segment, k = 20) + s(secchi_m, by = bay_segment) + s(temp, by = bay_segment) + s(sal, by = bay_segment), data = tomod, method = 'REML')

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
  facet_wrap(~var, scales = 'free')

## 
# interaction terms

mod <- gam(sgcov ~ te(secchi_m, yr, by = bay_segment) + te(temp, yr, by = bay_segment) + te(sal, yr, by = bay_segment), data = tomod, method = 'REML')

salslc <- data_slice(mod,
                     yr = c(2000, 2020),
                     sal = evenly(sal, n = 100), 
                     bay_segment = evenly(bay_segment)
)
toploslc <- smooth_estimates(mod, select = 'te(sal,yr)', data = salslc, partial_match = T) %>% 
  add_confint() %>% 
  filter(!((sal < min(tomod[tomod$bay_segment == 'HB', 'sal']) | sal > max(tomod[tomod$bay_segment == 'HB', 'sal'])) & bay_segment == 'HB')) %>%
  filter(!((sal < min(tomod[tomod$bay_segment == 'OTB', 'sal']) | sal > max(tomod[tomod$bay_segment == 'OTB', 'sal'])) & bay_segment == 'OTB')) %>%
  filter(!((sal < min(tomod[tomod$bay_segment == 'MTB', 'sal']) | sal > max(tomod[tomod$bay_segment == 'MTB', 'sal'])) & bay_segment == 'MTB')) 

p1 <- ggplot(toploslc, aes(x = sal, y = .estimate, color = factor(yr), fill = factor(yr), group = yr)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  facet_wrap(~bay_segment) +#, scales = 'free_y') +
  labs(
    x = NULL, 
    y = 'Sal partial'
  )
tempslc <- data_slice(mod,
                      yr = c(2000, 2020),
                      temp = evenly(temp, n = 100), 
                      bay_segment = evenly(bay_segment)
)
toploslc <- smooth_estimates(mod, select = 'te(temp,yr)', data = tempslc, partial_match = T) %>% 
  add_confint() %>% 
  filter(!((temp < min(tomod[tomod$bay_segment == 'HB', 'temp']) | temp > max(tomod[tomod$bay_segment == 'HB', 'temp'])) & bay_segment == 'HB')) %>%
  filter(!((temp < min(tomod[tomod$bay_segment == 'OTB', 'temp']) | temp > max(tomod[tomod$bay_segment == 'OTB', 'temp'])) & bay_segment == 'OTB')) %>%
  filter(!((temp < min(tomod[tomod$bay_segment == 'MTB', 'temp']) | temp > max(tomod[tomod$bay_segment == 'MTB', 'temp'])) & bay_segment == 'MTB'))

p2 <- ggplot(toploslc, aes(x = temp, y = .estimate, color = factor(yr), fill = factor(yr), group = yr)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  facet_wrap(~bay_segment) +#, scales = 'free_y') + 
  labs(
    x = NULL, 
    y = 'Temp partial'
  )
secchislc <- data_slice(mod,
                    yr = c(2000, 2020),
                    secchi_m = evenly(secchi_m, n = 100), 
                    bay_segment = evenly(bay_segment)
)
toploslc <- smooth_estimates(mod, select = 'te(secchi_m,yr)', data = secchislc, partial_match = T) %>% 
  add_confint() %>% 
  filter(!((secchi_m < min(tomod[tomod$bay_segment == 'HB', 'secchi_m']) | secchi_m > max(tomod[tomod$bay_segment == 'HB', 'secchi_m'])) & bay_segment == 'HB')) %>%
  filter(!((secchi_m < min(tomod[tomod$bay_segment == 'OTB', 'secchi_m']) | secchi_m > max(tomod[tomod$bay_segment == 'OTB', 'secchi_m'])) & bay_segment == 'OTB')) %>%
  filter(!((secchi_m < min(tomod[tomod$bay_segment == 'MTB', 'secchi_m']) | secchi_m > max(tomod[tomod$bay_segment == 'MTB', 'secchi_m'])) & bay_segment == 'MTB'))

p3 <- ggplot(toploslc, aes(x = secchi_m, y = .estimate, color = factor(yr), fill = factor(yr), group = yr)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  facet_wrap(~bay_segment) +#, scales = 'free_y') + 
  labs(
    x = NULL, 
    y = 'Secchi partial'
  )

p1 + p2 + p3 + plot_layout(ncol = 1, guides = 'collect')

# PINCO eval----------------------------------------------------------------------------------

data("pincotemp")

tomod <- pincotemp %>% 
  filter(!bay_segment %in% c('LTB')) %>% 
  select(date, allsg, yr, la, temp, sal, secchi_m, secchi_on_bottom, bay_segment) %>% 
  na.omit() %>%
  # filter(month(date) %in% c(7:10)) %>%
  summarise(
    allsg = sum(allsg) / length(allsg),
    la = mean(la),
    temp = mean(temp),
    sal = mean(sal),
    secchi_m = mean(secchi_m),
    .by = c(bay_segment, yr)
  )
# filter(!secchi_on_bottom)

##
# single smoothers only

mod <- gam(allsg ~ s(yr) + s(secchi_m) + s(temp) + s(sal), data = tomod, method = 'REML')

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

## 
# interaction terms

mod <- gam(allsg ~ te(secchi_m, yr) + te(temp, yr) + te(sal, yr), data = tomod, method = 'REML')

salslc <- data_slice(mod,
                     yr = c(2000, 2020),
                     sal = evenly(sal, n = 100)
)
toploslc <- smooth_estimates(mod, select = 'te(sal,yr)', data = salslc, partial_match = T) %>% 
  add_confint()

p1 <- ggplot(toploslc, aes(x = sal, y = .estimate, color = factor(yr), fill = factor(yr), group = yr)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  labs(
    x = NULL, 
    y = 'Sal partial'
  )
tempslc <- data_slice(mod,
                      yr = c(2000, 2020),
                      temp = evenly(temp, n = 100)
)
toploslc <- smooth_estimates(mod, select = 'te(temp,yr)', data = tempslc, partial_match = T) %>% 
  add_confint()

p2 <- ggplot(toploslc, aes(x = temp, y = .estimate, color = factor(yr), fill = factor(yr), group = yr)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+  
  labs(
    x = NULL, 
    y = 'Temp partial'
  )
secchislc <- data_slice(mod,
                        yr = c(2000, 2020),
                        secchi_m = evenly(secchi_m, n = 100)
)
toploslc <- smooth_estimates(mod, select = 'te(secchi_m,yr)', data = secchislc, partial_match = T) %>% 
  add_confint()

p3 <- ggplot(toploslc, aes(x = secchi_m, y = .estimate, color = factor(yr), fill = factor(yr), group = yr)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line()+ 
  labs(
    x = NULL, 
    y = 'Secchi partial'
  )

p1 + p2 + p3 + plot_layout(ncol = 1, guides = 'collect')