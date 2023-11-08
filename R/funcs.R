# function for formatting p-values in tables
p_ast <- function(x){
  
  sig_cats <- c('**', '*', '')
  sig_vals <- c(-Inf, 0.005, 0.05, Inf)
  
  out <- cut(x, breaks = sig_vals, labels = sig_cats, right = FALSE)
  out <- as.character(out)
  
  return(out)
  
}

# function for formatting p-values in tables
p_ast2 <- function(x){
  
  sig_cats <- c('**', '*', '$^.$', '')
  sig_vals <- c(-Inf, 0.005, 0.05, 0.1, Inf)
  
  out <- cut(x, breaks = sig_vals, labels = sig_cats, right = FALSE)
  out <- as.character(out)
  
  return(out)
  
}

# function for formatting p-values in text
p_txt <- function(x){
  
  sig_cats <- c('p < 0.005', 'p < 0.05', 'ns')
  sig_vals <- c(-Inf, 0.005, 0.05, Inf)
  
  out <- cut(x, breaks = sig_vals, labels = sig_cats, right = FALSE)
  out <- as.character(out)
  
  return(out)
  
}

#define function to extract overall p-value of model
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# summarize lm for inline text
modtxt_fun <- function(mod){
  r2 <- paste0('Adj. R$^2$ = ', round(summary(mod)$adj.r.squared, 2))
  fv <- round(summary(mod)$fstatistic, 2)
  fv <- paste0('F = ', fv[1], ', df = ', fv[2], ', ', fv[3])
  pv <- p_txt(overall_p(mod))
  out <- paste(r2, fv, pv, sep = ', ')
  return(out)
}

# run length encoding for vector of TF, return longest run of T
runfunc <- function(cnt){
  lngth <- length(na.omit(cnt))
  # if(lngth < 300) {
  #   out <- NA
  # } else {
    x <- rle(cnt)
    rn <- x$lengths[x$values]
    if(length(rn) == 0)
      out <- 0
    else
      out <- max(rn)
    # }
  return(out)
}

# prep epc, fim, pinco data prior to modelign 
modprep <- function(dat){
  
  dat %>% 
    filter(!bay_segment %in% c('LTB')) %>% 
    mutate(
      yrcat = case_when(
        yr <= 2016 ~ 'pre', 
        yr > 2016 ~ 'post'
      ), 
      yrcat = factor(yrcat, levels = c('pre', 'post'), 
                     labels = c('Recovery (pre - 2016)', 'Decline (2016 - present)')), 
      bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB'))
    )
  
}

# plotting function for glm mods of seagrass change
modplo_fun <- function(mod, xlab1, ylab1, subtitle1 = NULL, title1 = NULL, 
                       xlab2, ylab2 = NULL, subtitle2 = NULL, title2 = NULL, isbinom = F, ismetric = F) {

  if(ismetric){
    salcond <- quantile(mod$model$sal, 0.5)
    tempcond <- quantile(mod$model$temp, 0.5)
  }
  
  if(!ismetric){
    salcond <- quantile(mod$model$sal, 0.5)
    tempcond <- quantile(mod$model$temp, 0.5)
  }
  
  toplo1 <- visreg(mod, xvar = 'temp', by = 'yrcat', plot = F, scale = 'response', cond = list(sal = salcond), data = mod$model)
  toplo1a <- toplo1$fit
  
  toplo2 <- visreg(mod, xvar = 'sal', by = 'yrcat', plot = F, scale = 'response', cond = list(temp = tempcond), data = mod$model)
  toplo2a <- toplo2$fit
  
  ylim <- c(min(toplo1a$visregLwr, toplo2a$visregLwr), max(toplo1a$visregUpr, toplo2a$visregUpr))
  
  p1 <- ggplot(toplo1a, aes(x = temp)) + 
    geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = 'red2', alpha = 0.5) + 
    geom_line(aes(y = visregFit), color = 'red2') + 
    facet_grid(~ yrcat, scales = 'free_y') + 
    coord_cartesian(ylim = ylim) +
    labs(
      y = ylab1,
      x = xlab1,
      subtitle = subtitle1, 
      title = title1
    ) + 
    theme_bw() + 
    theme(
      strip.background = element_blank(), 
      panel.grid.minor = element_blank()
    )
  
  p2 <- ggplot(toplo2a, aes(x = sal)) + 
    geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = 'dodgerblue2', alpha = 0.5) + 
    geom_line(aes(y = visregFit), color = 'dodgerblue2') + 
    facet_grid(~ yrcat, scales = 'free_y') + 
    coord_cartesian(ylim = ylim) +
    labs(
      y = ylab2,
      x = xlab2,
      subtitle = subtitle2, 
      title = title2
    ) + 
    theme_bw() + 
    theme(
      strip.background = element_blank(), 
      panel.grid.minor = element_blank()
    )
  
  if(!isbinom){
    toplo1b <- toplo1$res
    toplo2b <- toplo2$res
    p1 <- p1 +
      geom_point(data = toplo1b, aes(x = temp, y = visregRes), size = 0.5, alpha = 0.5)
    p2 <- p2 +
      geom_point(data = toplo2b, aes(x = sal, y = visregRes), size = 0.5, alpha = 0.5)
  }
  
  if(isbinom){
    rawdat <- mod$model
    p1 <- p1 +
      geom_rug(data = rawdat[rawdat$allsg == 0, ], aes(x = temp), sides = 'b', alpha = 0.5) +
      geom_rug(data = rawdat[rawdat$allsg == 1, ], aes(x = temp), sides = 't', alpha = 0.5) 
    p2 <- p2 +
      geom_rug(data = rawdat[rawdat$allsg == 0, ], aes(x = sal), sides = 'b', alpha = 0.5) +
      geom_rug(data = rawdat[rawdat$allsg == 1, ], aes(x = sal), sides = 't', alpha = 0.5)
  }
  
  p <- p1 + p2 + plot_layout(ncol = 2)
  
  return(p)
  
}
