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

mod <- gam(total ~ s(yr, by = bay_segment) + s(la, by = bay_segment) + s(temp, by = bay_segment) + s(sal, by = bay_segment) + ti(la, yr, by = bay_segment) + ti(temp, yr, by = bay_segment) + ti(sal, yr, by = bay_segment), data = tomod, method = 'REML')

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
draw(mod, select = 'ti(sal,yr)', ncol = 3, partial_match = TRUE)

mod <- gam(total ~ s(yr, by = bay_segment) + s(la, by = bay_segment) + s(both, by = bay_segment) + ti(la, yr, by = bay_segment) + ti(both, yr, by = bay_segment), data = tomod, method = 'REML')

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

draw(mod, select = 's(both)', ncol = 3, partial_match = TRUE)
draw(mod, select = 's(la)', ncol = 3, partial_match = TRUE)
draw(mod, select = 'ti(both,yr)', ncol = 3, partial_match = TRUE)

# FIM -----------------------------------------------------------------------------------------

tomod <- fimsgtempdat %>% 
  filter(!bay_segment %in% c('LTB')) %>%
  select(date, sgcov, sgpres, yr, temp, sal, secchi_m, secchi_on_bottom, bay_segment) %>% 
  na.omit() %>%
  mutate(mo = month(date)) %>% 
  filter(month(date) %in% c(7:11)) %>%
  summarise(
    sgcov = mean(sgcov),
    sgpres = mean(sgpres),
    temp = mean(temp),
    sal = mean(sal),
    secchi_m = mean(secchi_m),
    .by = c(bay_segment, yr)
  ) %>%
  mutate(
    la = dplyr::case_when(
      bay_segment %in% "OTB" ~ 1.49 / secchi_m,
      bay_segment %in% "HB" ~ 1.61 / secchi_m,
      bay_segment %in% "MTB" ~ 1.49 / secchi_m,
      bay_segment %in% "LTB" ~ 1.84 / secchi_m
    )
  )

mod <- gam(sgcov ~ s(yr, by = bay_segment) + s(temp, by = bay_segment) + s(sal, by = bay_segment) + 
             ti(sal, yr, by = bay_segment) + ti(temp, yr, by = bay_segment), data = tomod, method = 'REML')

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
draw(mod, select = 'ti(sal,yr)', ncol = 3, partial_match = TRUE)

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
  filter(mo >= 7 & mo <= 11) %>% 
  summarise(
    allsg = sum(allsg) / length(allsg),
    temp = mean(temp),
    sal = mean(sal),
    secchi_m = mean(secchi_m),
    .by = c(yr)
  ) %>%
  mutate(
    la = 1.49 / secchi_m
  )

mod <- gam(allsg ~ s(yr) + s(temp) + s(sal) + ti(yr, temp) + ti(yr, sal), 
            data = tomod, method = 'REML')

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

draw(mod, select = 's(la)', partial_match = TRUE, constant = coef(mod)[1], fun = inv_link(mod))

draw(mod, select = 's(yr)', partial_match = TRUE, constant = coef(mod)[1], fun = inv_link(mod))

draw(mod, select = 's(temp)', partial_match = TRUE)
draw(mod, select = 's()', ncol = 4, partial_match = TRUE)
draw(mod, select = 'ti(yr,temp)', partial_match = TRUE)
draw(mod, select = 'ti(yr,sal)', partial_match = TRUE)

tempslc <- smooth_estimates(mod, select = 'ti(yr,sal)', partial_match = T) %>% 
  add_confint() %>% 
  filter(yr == 2022)
ggplot(tempslc, aes(x = sal, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, color = NA) +
  geom_line() + 
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
