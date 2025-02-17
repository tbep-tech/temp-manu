library(tidyverse)
library(tbeptools)

data(epcdata)
data(saliprd)
data(tempprd)

sta <- 66

epcdat <- epcdata %>% 
  filter(epchc_station == sta) %>% 
  mutate(
    date = as.Date(SampleTime)
  ) %>% 
  filter(yr < 2023)

toplo1 <- tempprd %>% 
  filter(station == sta) %>% 
  filter(param == 'tempbot') %>% 
  unnest(prd)
toplo2 <- saliprd %>% 
  filter(station == sta) %>% 
  filter(param == 'salibot') %>% 
  unnest(prd)

p1a <- ggplot() + 
  geom_point(data = epcdat, aes(x = date, y = Temp_Water_Bottom_degC)) + 
  # geom_line(data = toplo1, aes(x = date, y = value), color = '#FF2A2A', linewidth = 1) + 
  scale_x_date(limits = c(as.Date('1975-01-01'), as.Date('2023-01-01'))) +
  theme_minimal() + 
  theme(
    panel.grid.minor = element_blank()
  ) + 
  labs(
    x = NULL, 
    y = expression('Bottom temp. ('*degree*C*')')
  )
p1b <- ggplot() + 
  geom_line(data = toplo1, aes(x = date, y = value), color = '#FF2A2A', linewidth = 1) + 
  scale_x_date(limits = c(as.Date('1975-01-01'), as.Date('2023-01-01'))) +
  theme_void()

p2a <- ggplot() + 
  geom_point(data = epcdat, aes(x = date, y = Sal_Bottom_ppth)) + 
  # geom_line(data = toplo2, aes(x = date, y = value), color = '#7EAEFE', linewidth = 1) + 
  scale_x_date(limits = c(as.Date('1975-01-01'), as.Date('2023-01-01'))) +
  theme_minimal() + 
  theme(
    panel.grid.minor = element_blank()
  ) + 
  labs(
    x = NULL, 
    y = expression('Bottom salinity (psu)')
  )
p2b <- ggplot() + 
  geom_line(data = toplo2, aes(x = date, y = value), color = '#7EAEFE', linewidth = 1) +
  scale_x_date(limits = c(as.Date('1975-01-01'), as.Date('2023-01-01'))) +
  theme_void() 

png('~/Desktop/tempprda.png', width = 6, height = 2.5, units = 'in', res = 300, bg = 'transparent')
print(p1a)
dev.off()

png('~/Desktop/tempprdb.png', width = 6, height = 2.5, units = 'in', res = 300, bg = 'transparent')
print(p1b)
dev.off()

png('~/Desktop/saliprda.png', width = 6, height = 2.5, units = 'in', res = 300, bg = 'transparent')
print(p2a)
dev.off()

png('~/Desktop/saliprdb.png', width = 6, height = 2.5, units = 'in', res = 300, bg = 'transparent')
print(p2b)
dev.off()
