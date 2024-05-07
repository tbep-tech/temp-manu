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
 
  n <- paste('n =', nrow(mod$model))
  smmod <- summary(mod)
  rsq <- paste0('R.sq. ', sprintf('%.2f', round(smmod$r.sq, 2)))
  dev <- paste0('Deviance explained ', round(100 * smmod$dev.expl, 0), '%')
  cap <- paste0(cap, ' ', n, ', ', rsq, ', ', dev)
  
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

# plot all s() for gam
gamplo_fun <- function(mod, smths, labels, cols, bayseg = T){

  out <- tibble(smooths = factor(smths, levels = smths)) %>% 
    group_nest(smooths) %>% 
    mutate(
      data = purrr::map(smooths, function(x) smooth_estimates(mod, as.character(x), partial_match = T) %>% add_confint()), 
      xlab = factor(smooths, levels = smths,
                    labels = labels
      ),
      xvar = gsub('s\\((.*)\\)', '\\1', smooths),
      cols = cols,
      plos = purrr::pmap(list(data, xvar, xlab, cols), function(data, xvar, xlab, cols){

        p <- ggplot(data, aes(x = get(xvar), y = .estimate)) +
          geom_hline(yintercept = 0, linetype = 'dashed') +
          geom_ribbon(aes(ymin  = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = cols) +
          geom_line(color = cols) +
          theme_minimal() +
          theme(
            panel.grid.minor = element_blank()
          ) +
          labs(
            x = xlab,
            y = 'Partial effect'
          )

        if(bayseg)
          p <- p +
            facet_wrap(~ bay_segment, ncol = 3, scales = 'free_x')

        return(p)

      })
    )
  
  return(out)
  
}

