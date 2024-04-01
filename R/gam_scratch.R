library(mgcv)

data(epccmbdat)

tomod <- epccmbdat %>% 
  # select(-bbave, -chla) %>% 
  filter(!bay_segment %in% c('LTB', 'MTB', 'OTB'))

mod <- gam(total ~ s(yr) + s(la) + s(temp) + s(sal), data = tomod, method = 'REML')

draw(mod)
toplo <- smooth_estimates(mod) %>% 
  add_confint() %>% 
  select(.smooth, .estimate, yr, la, temp, sal, .lower_ci, .upper_ci) %>% 
  pivot_longer(cols = c(yr, la, temp, sal), names_to = 'var', values_to = 'value') %>% 
  filter(!is.na(value))

ggplot(toplo, aes(x = value, y = .estimate)) + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
  geom_line()+ 
  facet_wrap(~var, scales = 'free')

### repeat down to here but for each bay segment model

# slc <- data_slice(mod, 
#                     yr = evenly(yr, n = 3), 
#                     temp = evenly(temp, n = 100)
#                     )
# toplo <- fitted_values(mod, data = slc, scale = 'response')
# 
# ggplot(toplo, aes(x = temp, y = .fitted)) + 
#   geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
#   geom_line()+ 
#   facet_wrap(~yr)
