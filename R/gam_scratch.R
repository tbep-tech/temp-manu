library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)
library(here)
library(broom)
library(knitr)

source(here('R/funcs.R'))

load(file = here('data/epccmbdat.RData'))
load(file = here('data/fimsgtempdat.RData'))
load(file = here('data/pincotemp.RData'))

# EPC -----------------------------------------------------------------------------------------

tomod <- epccmbdat %>% 
  filter(!bay_segment %in% c('LTB'))

mod <- gam(total ~ s(yr, by = bay_segment) + s(la, by = bay_segment) + s(temp, by = bay_segment) + s(sal, by = bay_segment), data = tomod, method = 'REML')

mod %>% 
  tidy() %>% 
  mutate(
    star = p_ast2(p.value)
  ) %>% 
  kable(digits = 3)

smmod <- summary(mod)
smtab <- tibble(
  R.sq = smmod$r.sq,
  `Deviance explained` = paste0(round(100*smmod$dev.expl, 0), '%')
)
kable(smtab, digits = 2, align = 'l')

draw(mod, ncol = 3)

mod <- gam(total ~ s(yr, by = bay_segment) + s(la, by = bay_segment) + s(both, by = bay_segment), data = tomod, method = 'REML')

mod %>% 
  tidy() %>% 
  mutate(
    star = p_ast2(p.value)
  ) %>% 
  kable(digits = 3)

smmod <- summary(mod)
smtab <- tibble(
  R.sq = smmod$r.sq,
  `Deviance explained` = paste0(round(100*smmod$dev.expl, 0), '%')
)
kable(smtab, digits = 2, align = 'l')

draw(mod, ncol = 3)

# FIM -----------------------------------------------------------------------------------------

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

mod <- gam(sgcov ~ s(yr, by = bay_segment) + s(la, by = bay_segment) + s(temp, by = bay_segment) + s(sal, by = bay_segment) + s(mo, bs = 'cc', by = bay_segment) + ti(yr, mo, bs = c('tp', 'cc'), by = bay_segment) + ti(la, mo, bs = c('tp', 'cc'), by = bay_segment) + ti(temp, mo, bs = c('tp', 'cc'), by = bay_segment) + ti(sal, mo, bs = c('tp', 'cc'), by = bay_segment), knots = list(doy = c(1, 12)), data = tomod, method = 'REML')

mod %>% 
  tidy() %>% 
  mutate(
    star = p_ast2(p.value)
  ) %>% 
  kable(digits = 3)

smmod <- summary(mod)
smtab <- tibble(
  R.sq = smmod$r.sq,
  `Deviance explained` = paste0(round(100*smmod$dev.expl, 0), '%')
)
kable(smtab, digits = 2, align = 'l')

draw(mod, select = 's(temp)', ncol = 3, partial_match = TRUE)
draw(mod, select = 's(sal)', ncol = 3, partial_match = TRUE)
draw(mod, select = 'ti(temp,mo)', ncol = 3, partial_match = TRUE)
draw(mod, select = 'ti(sal,mo)', ncol = 3, partial_match = TRUE)

tempslc <- smooth_estimates(mod, select = 'ti(temp,mo)', partial_match = T) %>% 
  add_confint() %>% 
  filter(mo %in% c(1, 7)) %>% 
  filter(bay_segment %in% c('HB', 'MTB'))
ggplot(tempslc, aes(x = temp, y = .estimate, color = factor(mo), fill = factor(mo), group = mo)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line() +
  facet_wrap(~bay_segment, scales = 'free_y') + 
  labs(
    x = 'Temperature',
    y = 'Partial effect',
    color = 'Month',
    fill = 'Month'
  )

salslc <- smooth_estimates(mod, select = 'ti(sal,mo)', partial_match = T) %>% 
  add_confint() %>% 
  filter(mo %in% c(3, 9)) %>% 
  filter(bay_segment %in% c('MTB'))
ggplot(salslc, aes(x = sal, y = .estimate, color = factor(mo), fill = factor(mo), group = mo)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line() +
  facet_wrap(~bay_segment, scales = 'free_y') + 
  labs(
    x = 'Salinity',
    y = 'Partial effect',
    color = 'Month',
    fill = 'Month'
  )

# PDEM ----------------------------------------------------------------------------------------

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

mod <- gam(allsg ~ s(yr, by = site) + s(la, by = site) + s(temp, by = site) + s(sal, by = site) + s(mo, bs = 'cc', by = site) + ti(yr, mo, bs = c('tp', 'cc'), by = site) + ti(la, mo, bs = c('tp', 'cc'), by = site) + ti(temp, mo, bs = c('tp', 'cc'), by = site) + ti(sal, mo, bs = c('tp', 'cc'), by = site), knots = list(doy = c(1, 12)), data = tomod, method = 'REML')

mod %>% 
  tidy() %>% 
  mutate(
    star = p_ast2(p.value)
  ) %>% 
  kable(digits = 3)

smmod <- summary(mod)
smtab <- tibble(
  R.sq = smmod$r.sq,
  `Deviance explained` = paste0(round(100*smmod$dev.expl, 0), '%')
)
kable(smtab, digits = 2, align = 'l')

draw(mod, select = 's(temp)', ncol = 4, partial_match = TRUE)
draw(mod, select = 's(sal)', ncol = 4, partial_match = TRUE)
draw(mod, select = 'ti(temp,mo)', ncol = 4, partial_match = TRUE)
draw(mod, select = 'ti(sal,mo)', ncol = 4, partial_match = TRUE)

tempslc <- smooth_estimates(mod, select = 'ti(temp,mo)', partial_match = T) %>% 
  add_confint() %>% 
  filter(mo %in% c(1, 7)) %>% 
  filter(site %in% c('E1', 'E4'))
ggplot(tempslc, aes(x = temp, y = .estimate, color = factor(mo), fill = factor(mo), group = mo)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line() +
  facet_wrap(~site, scales = 'free_y') + 
  labs(
    x = 'Temperature',
    y = 'Partial effect',
    color = 'Month',
    fill = 'Month'
  )

salslc <- smooth_estimates(mod, select = 'ti(sal,mo)', partial_match = T) %>% 
  add_confint() %>% 
  filter(mo %in% c(3, 9)) %>% 
  filter(site %in% c('E4'))
ggplot(salslc, aes(x = sal, y = .estimate, color = factor(mo), fill = factor(mo), group = mo)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line() +
  facet_wrap(~site, scales = 'free_y') + 
  labs(
    x = 'Salinity',
    y = 'Partial effect',
    color = 'Month',
    fill = 'Month'
  )
