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
p_txt <- function(x, addp = TRUE){
  
  if(x < 0.001){
    out <- 'p < 0.001'
  } else {
    out <- paste('p =', sprintf('%.3f', round(x, 3)))
  }
  
  if(!addp)
    out <- gsub('^p\\s|\\=\\s', '', out)
  
  return(out)
  
}

# gam table function
gam_table <- function(mod, cap = NULL){
  
  out <- mod %>% 
    tidy() %>% 
    rowwise() %>% 
    mutate(
      p.value = p_txt(p.value, addp = F)
    ) %>% 
    ungroup() %>% 
    mutate_if(is.numeric, round, 3) %>%
    rename(
      'Smoother' = term,
      'Ref.df' = ref.df,
      'F' = statistic,
      'p' = p.value
    ) %>%
    kable(digits = 3, align = 'lrrrr', caption = cap)
  
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
modtxt_fun <- function(modsum){
  
  lapply(modsum, function(x){
    n <- paste0('n = ', x$n)
    r2 <- paste0('Adj. R$^2$ = ', x$rsq)
    dev <- paste0('Deviance explained = ', x$dev, '\\%')
    out <- paste(n, r2, dev, sep = ', ')
    return(out)
  })
  
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

# gam plot for fim and pdem, main text
gamplo_fun <- function(sgmods){
  
  mod <- sgmods$fimmod
  smths <- c('s(yr)', 's(temp)', 's(sal)')
  labels <- c('Year', 'Temperature (\u00B0C)', 'Salinity (ppt)')
  cols <- c('grey20', 'red2', 'dodgerblue2')
  
  toplo1 <- tibble(smooths = factor(smths, levels = smths)) %>% 
    group_nest(smooths) %>% 
    mutate(
      data = purrr::map(smooths, function(x) smooth_estimates(mod, as.character(x), partial_match = T) %>% add_confint()),
      xlab = factor(smooths, levels = smths,
                    labels = labels
      ),
      xvar = gsub('s\\((.*)\\)', '\\1', smooths), 
      cols = cols
    ) %>% 
    unnest('data') %>% 
    group_by(smooths, bay_segment, xlab, xvar, cols) %>% 
    nest() %>% 
    mutate(
      plos = purrr::pmap(list(data, xvar, xlab, cols, bay_segment), function(data, xvar, xlab, cols, bay_segment){
        
        # xlab <- ifelse(bay_segment == 'HB', as.character(xlab), "")
        xlab <- as.character(xlab)
        ylab <- ifelse(bay_segment == 'OTB', 'Partial effect', '')
        xsz <- ifelse(xvar == 'yr', 8, 8)
        
        yrng <- 0.9 * max(abs(range(data$.estimate, data$.lower_ci, data$.upper_ci)))
        
        if(bay_segment == 'HB' & xvar == 'sal')
          yrng <-100
        
        ttl <- ifelse(bay_segment == 'OTB' & xvar == 'yr', 
                      '(a) FIM', '')
        
        bay_segment <- ifelse(xvar == 'yr', as.character(bay_segment), '')
        
        p <- ggplot(data, aes(x = get(xvar), y = .estimate)) +
          geom_hline(yintercept = 0, linetype = 'dashed') +
          geom_ribbon(aes(ymin  = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = cols) +
          geom_line(color = cols) +
          coord_cartesian(ylim = c(-yrng, yrng)) +
          theme_minimal() +
          theme(
            panel.grid.minor = element_blank(), 
            plot.subtitle = element_text(hjust = 0.5, size = 12), 
            axis.text.x = element_text(size = xsz),
            plot.margin = margin(0, 0, 0, 0)
          ) +
          labs(
            subtitle = bay_segment, 
            title = ttl,
            x = xlab,
            y = ylab
          )
        
        return(p)
        
      })
    )
  
  mod <- sgmods$pincomod
  smths <- c('s(yr)', 's(temp)', 's(sal)')
  labels <- c('Year', 'Temperature (\u00B0C)', 'Salinity (ppt)')
  cols <- c('grey20', 'red2', 'dodgerblue2')
  
  toplo2 <- tibble(smooths = factor(smths, levels = smths)) %>% 
    group_nest(smooths) %>% 
    mutate(
      data = purrr::map(smooths, function(x) smooth_estimates(mod, as.character(x), partial_match = T) %>% add_confint()),
      xlab = factor(smooths, levels = smths,
                    labels = labels
      ),
      xvar = gsub('s\\((.*)\\)', '\\1', smooths), 
      cols = cols,
      plos = purrr::pmap(list(data, xvar, xlab, cols), function(data, xvar, xlab, cols){
        
        xlab <- as.character(xlab)
        # ylab <- ifelse(xvar == 'yr', 'Partial effect', '')
        ylab <- ''
        xsz <- ifelse(xvar == 'yr', 8, 8)
        
        yrng <- 0.9 * max(abs(range(data$.estimate, data$.lower_ci, data$.upper_ci)))
        
        ttl <- ifelse(xvar == 'yr', '(b) PDEM', '')
        bay_segment <- ifelse(xvar == 'yr', 'OTB', '')
        
        p <- ggplot(data, aes(x = get(xvar), y = .estimate)) +
          geom_hline(yintercept = 0, linetype = 'dashed') +
          geom_ribbon(aes(ymin  = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = cols) +
          geom_line(color = cols) +
          coord_cartesian(ylim = c(-yrng, yrng)) +
          theme_minimal() +
          theme(
            panel.grid.minor = element_blank(),
            plot.subtitle = element_text(hjust = 0.5, size = 12), 
            axis.text.x = element_text(size = xsz), 
            plot.margin = margin(0, 0, 0, 0)
          ) +
          labs(
            title = ttl,
            x = xlab,
            subtitle = bay_segment,
            y = ylab
          )
        
        return(p)
        
      })
    )
  
  p <- (toplo1$plos[[1]] + toplo1$plos[[2]] + toplo1$plos[[3]] + toplo2$plos[[1]] + plot_layout(ncol = 4)) / 
    (toplo1$plos[[4]] + toplo1$plos[[5]] + toplo1$plos[[6]] + toplo2$plos[[2]] + plot_layout(ncol = 4)) /
    (toplo1$plos[[7]] + toplo1$plos[[8]] + toplo1$plos[[9]] + toplo2$plos[[3]] + plot_layout(ncol = 4))
  
  return(p)
  
}

# plot all s() for gam, supplement
epcgamplo_fun <- function(mod, smths, labels, cols){

  out <- tibble(smooths = factor(smths, levels = smths)) %>% 
    group_nest(smooths) %>% 
    mutate(
      data = purrr::map(smooths, function(x) smooth_estimates(mod, as.character(x), partial_match = T) %>% add_confint()),
      xlab = factor(smooths, levels = smths,
                    labels = labels
      ),
      xvar = gsub('s\\((.*)\\)', '\\1', smooths), 
      cols = cols
    ) %>% 
    unnest('data') %>% 
    group_by(smooths, bay_segment, xlab, xvar, cols) %>% 
    nest() %>% 
    mutate(
      plos = purrr::pmap(list(data, xvar, xlab, cols, bay_segment), function(data, xvar, xlab, cols, bay_segment){
        
        xlab <- as.character(xlab)#ifelse(bay_segment == 'HB', as.character(xlab), "")
        ylab <- ifelse(bay_segment == 'OTB', 'Partial effect', '')
        bay_segment <- ifelse(xvar == 'yr', as.character(bay_segment), '')
        xsz <- ifelse(xvar == 'yr', 8, 8)
        
        yrng <- 0.9 * max(abs(range(data$.estimate, data$.lower_ci, data$.upper_ci)))
        
        p <- ggplot(data, aes(x = get(xvar), y = .estimate)) +
          geom_hline(yintercept = 0, linetype = 'dashed') +
          geom_ribbon(aes(ymin  = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = cols) +
          geom_line(color = cols) +
          coord_cartesian(ylim = c(-yrng, yrng)) +
          theme_minimal() +
          theme(
            panel.grid.minor = element_blank(), 
            plot.subtitle = element_text(hjust = 0.5, size = 12), 
            axis.text.x = element_text(size = xsz)#,
            # plot.margin = margin(0, 0, 0, 0)
          ) +
          labs(
            subtitle = bay_segment, 
            x = xlab,
            y = ylab
          )
        
        return(p)
        
      })
    )
  
    return(out)
  
}

