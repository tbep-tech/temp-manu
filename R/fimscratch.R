
seglng <- c('Old Tampa Bay', 'Hillsborough Bay', 'Middle Tampa Bay', 'Lower Tampa Bay')
segshr <- c('OTB', 'HB', 'MTB', 'LTB')

load(file = here('data/fimsgtempdat.RData'))

# relating sg presence to temp/sal ------------------------------------------------------------

toplo1 <- fimsgtempdat %>% 
  st_set_geometry(NULL) %>% 
  mutate(
    year = year(date),
    yrgroup = ifelse(year < 2000, '1996 - 1999', '2000 - 2022'),
    mo = month(date, label = T, abbr = F)
  ) %>% 
  filter(depth < 2) %>% 
  filter(mo %in% c('February', 'August')) %>% 
  pivot_longer(names_to = 'var', values_to = 'val', c(temp, sal)) %>% 
  summarise(
    avev = mean(val),
    hiv = t.test(val)$conf.int[2],
    lov = t.test(val)$conf.int[1], 
    .by = c('yrgroup', 'bay_segment', 'mo', 'var', 'sgpres')
  )

wd <- 0.5
ggplot(toplo1[toplo1$var == 'temp', ], aes(x = yrgroup, color = sgpres)) + 
  geom_point(aes(y = avev), position = position_dodge(width = wd), size = 1) + 
  geom_errorbar(aes(ymin = lov, ymax = hiv), position = position_dodge(width = wd), width = 0) +
  facet_grid(mo~bay_segment, scales = 'free_y')

toplo2 <- fimsgtempdat %>% 
  st_set_geometry(NULL) %>% 
  filter(depth < 2) %>%
  rename(
    Temp = temp, 
    Sal = sal
  ) %>% 
  mutate(
    yr = year(date),
    mo = month(date), 
    sgpres = ifelse(sgpres, 1, 0),
    yrcat = cut(yr, breaks = c(-Inf, 2016, Inf), labels = c('Recovery (pre - 2016)', 'Decline (2016 - present)'), right = F)
  ) %>% 
  filter(mo == 8)

# ggplot(toplo2[toplo2$mo %in% c(2, 4, 6, 8, 10, 12), ], aes(x = date, y = temp, col = factor(sgpres))) + 
#   geom_point() + 
#   facet_grid(mo~bay_segment, scales = 'free_y') + 
#   geom_smooth(method = 'lm', se = F)

tomod <- toplo2
mod <- glm(sgpres ~ Temp*yrcat*bay_segment + Sal*yrcat*bay_segment, data = tomod, family = binomial('logit')) %>% 
  step()


toplo2 <- getprd_fun(mod, depvar = 'Temp')
toplo3 <- getprd_fun(mod, depvar = 'Sal')
p2 <- ggplot(toplo2, aes(x = Temp)) + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = 'lightgrey', alpha = 0.5) + 
  geom_line(aes(y = visregFit)) + 
  geom_rug(data = tomod[tomod$sgpres == 0,], aes(x = Temp, y = sgpres), sides = 'b', linewidth = 1, color = 'blue') +
  geom_rug(data = tomod[tomod$sgpres == 1,], aes(x = Temp, y = sgpres), sides = 't', linewidth = 1, color = 'red') +
  facet_grid(bay_segment ~ yrcat) + 
  scale_y_continuous(labels = scales::percent) +
  labs(
    y = 'Probability of seagrass present', 
    x = 'Temperature'
  ) + 
  theme_bw() + 
  theme(
    strip.background = element_blank(), 
    panel.grid.minor = element_blank()
  )

p3 <- ggplot(toplo3, aes(x = Sal)) + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = 'lightgrey', alpha = 0.5) + 
  geom_line(aes(y = visregFit)) + 
  geom_rug(data = tomod[tomod$sgpres == 0,], aes(x = Sal, y = sgpres), sides = 'b', linewidth = 1, color = 'blue') +
  geom_rug(data = tomod[tomod$sgpres == 1,], aes(x = Sal, y = sgpres), sides = 't', linewidth = 1, color = 'red') +
  facet_grid(bay_segment ~ yrcat) + 
  scale_y_continuous(labels = scales::percent) +
  labs(
    y = 'Probability of seagrass present', 
    x = 'Salinity'
  ) + 
  theme_bw() + 
  theme(
    strip.background = element_blank(), 
    panel.grid.minor = element_blank()
  )



toplo2 <- fimsgtempdat %>% 
  st_transform(crs = st_crs(sgmanagement)) %>% 
  .[sgmanagement[sgmanagement$areas == 21, ], ] %>% 
  st_set_geometry(NULL) %>% 
  # filter(depth < 2) %>% 
  rename(
    Temp = temp, 
    Sal = sal
  ) %>% 
  mutate(
    yr = year(date),
    mo = month(date), 
    sgpres = ifelse(sgpres, 1, 0),
    yrcat = cut(yr, breaks = c(-Inf, 2016, Inf), labels = c('Recovery (pre - 2016)', 'Decline (2016 - present)'), right = F)
  ) %>% 
  filter(mo %in% c(7, 8, 9))

tomod <- toplo2
mod <- glm(sgpres ~ Temp*yrcat + Sal*yrcat, data = tomod, family = binomial('logit'))

toplo2 <- visreg(mod, 'Temp', by = 'yrcat', scale = 'response', plot = F)$fit
toplo3 <- visreg(mod, 'Sal', by = 'yrcat', scale = 'response', plot = F)$fit

p2 <- ggplot(toplo2, aes(x = Temp)) + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = 'lightgrey', alpha = 0.5) + 
  geom_line(aes(y = visregFit)) + 
  geom_rug(data = tomod[tomod$sgpres == 0,], aes(x = Temp, y = sgpres), sides = 'b', linewidth = 1, color = 'blue') +
  geom_rug(data = tomod[tomod$sgpres == 1,], aes(x = Temp, y = sgpres), sides = 't', linewidth = 1, color = 'red') +
  facet_grid(~ yrcat) + 
  scale_y_continuous(labels = scales::percent) +
  labs(
    y = 'Probability of seagrass present', 
    x = 'Temperature'
  ) + 
  theme_bw() + 
  theme(
    strip.background = element_blank(), 
    panel.grid.minor = element_blank()
  )

p3 <- ggplot(toplo3, aes(x = Sal)) + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = 'lightgrey', alpha = 0.5) + 
  geom_line(aes(y = visregFit)) + 
  geom_rug(data = tomod[tomod$sgpres == 0,], aes(x = Sal, y = sgpres), sides = 'b', linewidth = 1, color = 'blue') +
  geom_rug(data = tomod[tomod$sgpres == 1,], aes(x = Sal, y = sgpres), sides = 't', linewidth = 1, color = 'red') +
  facet_grid(~ yrcat) + 
  scale_y_continuous(labels = scales::percent) +
  labs(
    y = 'Probability of seagrass present', 
    x = 'Salinity'
  ) + 
  theme_bw() + 
  theme(
    strip.background = element_blank(), 
    panel.grid.minor = element_blank()
  )


# temp/sal trends -----------------------------------------------------------------------------

load(file = here('data/fimsgtempdat.RData'))

col_fun <- function(data, min_val, max_val, n = 100) {
  # Normalize the data to be between -1 and 1
  normalized_data <- (data - min_val) / (max_val - min_val) * 2 - 1
  
  # Create a divergent color palette
  palette <- colorRampPalette(c('blue', 'grey', 'tomato1'))(n)
  
  # Map normalized data to colors
  color_values <- ifelse(normalized_data < 0,
                         palette[ceiling((normalized_data + 1) * (n / 2))],
                         palette[ceiling(normalized_data * (n / 2) + (n / 2))])
  
  return(color_values)
}

# aggregate by bay segment
sumfimdat <- fimsgtempdat %>% 
  st_set_geometry(NULL) %>% 
  mutate(
    yr = year(date), 
    mo = month(date)
  ) %>% 
  summarise(
    temp = mean(temp), 
    sal = mean(sal),
    .by = c(yr, mo, bay_segment)
  ) %>% 
  pivot_longer(names_to = 'var', values_to = 'val', temp:sal)

# kendall all years
sktres <- sumfimdat %>% 
  nest(.by = c('bay_segment', 'var')) %>% 
  mutate(
    skt = purrr::pmap(list(bay_segment, var, data), function(bay_segment, var, data){
      
      cat(bay_segment, var, '\n')
      
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

totab <- sktres %>% 
  mutate(
    pval = p_ast(pval), 
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) %>% 
  arrange(var, bay_segment)

reactable(
  totab,
  columns = list(
    slos = colDef(
      style =  function(value) {
        color <- col_fun(value, min_val = -0.201, max_val = 0.201)
        list(background = color)
      }, 
      format = colFormat(digits = 2)
    )
  )
)

##
# by sg pres/abs

# aggregate by bay segment, uneven month distribution across sg management areas
sumfimdat <- fimdat %>% 
  st_set_geometry(NULL) %>% 
  mutate(
    yr = year(date), 
    mo = month(date)
  ) %>% 
  summarise(
    temp = mean(temp), 
    sal = mean(sal),
    .by = c(yr, mo, bay_segment, sgpres)
  ) %>% 
  pivot_longer(names_to = 'var', values_to = 'val', temp:sal)

# kendall all years by sg pres/abs
sktres <- sumfimdat %>% 
  nest(.by = c('bay_segment', 'sgpres', 'var')) %>% 
  mutate(
    skt = purrr::pmap(list(bay_segment, sgpres, var, data), function(bay_segment, sgpres,var, data){
      
      cat(bay_segment, sgpres, var, '\n')
      
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

totab <- sktres %>% 
  mutate(
    pval = p_ast(pval), 
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
  ) %>% 
  arrange(var, sgpres, bay_segment)

reactable(
  totab,
  columns = list(
    slos = colDef(
      style =  function(value) {
        color <- col_fun(value, min_val = -0.201, max_val = 0.201)
        list(background = color)
      }, 
      format = colFormat(digits = 2)
    )
  )
)
