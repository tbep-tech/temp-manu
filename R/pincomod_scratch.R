data(pincotemp)

tomod <- pincotemp %>%
  mutate(
    yrcat = case_when(
      yr <= 2015 ~ 'Recovery (pre - 2016)',
      T ~ 'Decline (2016 - present)'
    ), 
    yrcat = factor(yrcat, levels = c('Recovery (pre - 2016)', 'Decline (2016 - present)'))
  ) #%>% 
  # filter(yr > 2008) %>%
  # filter(
  # month(date) >= 7 & month(date) <= 11
  # ) #%>%
  # filter(site == 'E1') %>%
  # filter(hr %in% c(9:12))

# models
pincosalmod <- gam(allsg ~ yrcat + s(sal, by = yrcat), data = tomod, family = binomial(link = 'logit'),
                 method = 'REML')

pincotempmod <- gam(allsg ~ yrcat + s(temp, by = yrcat), data = tomod, family = binomial(link = 'logit'),
                  method = 'REML')

toplo1 <- getprd_fun3(pincosalmod, 'sal')
toplo2 <- getprd_fun3(pincotempmod, 'temp')

thm <- theme_bw() + 
  theme(
    panel.grid.minor = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(size = 12), 
    axis.text = element_text(size = 11), 
    axis.title = element_text(size = 12)
  )

# separate plots by bay segment
p1 <- ggplot(toplo1, aes(x = sal, y = visregFit)) + 
  geom_rug(data = tomod[tomod$allsg == 0,], aes(x = sal, y = allsg), sides = 'b', linewidth = 1, color = 'red') +
  geom_rug(data = tomod[tomod$allsg == 1,], aes(x = sal, y = allsg), sides = 't', linewidth = 1, color = 'blue') +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = 'dodgerblue2', alpha = 0.2) +
  geom_line(color = 'dodgerblue2', size = 1) +
  coord_cartesian(
    ylim = c(0, max(toplo1$visregUpr))#,
    # xlim = c(10, 30)
  ) +
  facet_grid(~ yrcat, scales = 'free') +
  thm + 
  labs(
    x = 'Salinity (psu)', 
    y = 'Probability of seagrass present', 
    title = '(a) Seagrass presence vs salinity'
  )

p2 <- ggplot(toplo2, aes(x = temp, y = visregFit)) + 
  geom_rug(data = tomod[tomod$allsg == 0,], aes(x = temp, y = allsg), sides = 'b', linewidth = 1, color = 'red') +
  geom_rug(data = tomod[tomod$allsg == 1,], aes(x = temp, y = allsg), sides = 't', linewidth = 1, color = 'blue') +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = 'red2', alpha = 0.2) +
  geom_line(color = 'red2', size = 1) +
  coord_cartesian(ylim = c(0, max(toplo2$visregUpr))) +
  facet_grid(~ yrcat, scales = 'free') +
  coord_cartesian(
    ylim = c(0, 0.35)#,
    # xlim = c(10, max(pincotemp$temp, na.rm = T))
  ) +
  thm + 
  labs(
    x = 'Water temp. (\u00B0 C)', 
    y = 'Prob of seagrass present',
    title = '(b) Seagrass presence vs temperature'
  )

p <- p1 + p2 + plot_layout(ncol = 1)
p

# data(pincotemp)
# 
# tomod <- pincotemp %>%
#   mutate(
#     yrcat = case_when(
#       yr <= 2015 ~ 'Recovery (pre - 2016)',
#       T ~ 'Decline (2016 - present)'
#     ), 
#     yrcat = factor(yrcat, levels = c('Recovery (pre - 2016)', 'Decline (2016 - present)'))#, 
#     # yr = factor(yr, ordered = T)
#   )
# 
# # models
# pincomod <- glm(allsg ~ yrcat*temp*sal, data = tomod, family = binomial(link = 'logit')) %>% 
#   step()
# 
# visreg(pincomod, xvar = 'temp', by = 'sal',  breaks = c(15,20, 30), cond = list(yrcat = 'Decline (2016 - present)'), scale = 'response')
# visreg(pincomod, xvar = 'temp', by = 'sal',  breaks = c(15,20, 30), cond = list(yrcat = 'Recovery (pre - 2016)'), scale = 'response')
