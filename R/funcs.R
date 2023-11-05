# function for formatting p-values in tables
p_ast <- function(x){
  
  sig_cats <- c('**', '*', '')
  sig_vals <- c(-Inf, 0.005, 0.05, Inf)
  
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

# get model predictions glm
getprd_fun <- function(modin, depvar = 'Sal'){
  
  fitotb <- visreg(modin, depvar, by = 'yrcat', cond = list(bay_segment = 'OTB'), scale = 'response', plot = F) %>% 
    .$fit %>%  
    mutate(
      bay_segment = 'OTB'
    )
  fithb <- visreg(modin, depvar, by = 'yrcat', cond = list(bay_segment = 'HB'), scale = 'response', plot = F) %>% 
    .$fit %>% 
    mutate(
      bay_segment = 'HB'
    ) 
  fitmtb <- visreg(modin, depvar, by = 'yrcat', cond = list(bay_segment = 'MTB'), scale = 'response', plot = F) %>% 
    .$fit %>% 
    mutate(
      bay_segment = 'MTB'
    )
  
  out <- bind_rows(fitotb, fithb, fitmtb) %>% 
    mutate(
      bay_segment = factor(bay_segment, levels = segshr)
    )
  
  return(out)
  
}

# get model predictions gam, fim
getprd_fun2 <- function(mod, indvar = c('sal', 'temp')){
  
  indvar <- match.arg(indvar)
  
  out <- unique(mod$model[, c('bay_segment', 'yrcat')]) %>% 
    group_by(bay_segment, yrcat) %>% 
    nest() %>% 
    mutate(
      data = purrr::pmap(list(bay_segment, yrcat), function(bay_segment, yrcat){
        
        condls <- list(bay_segment = bay_segment, yrcat = yrcat, fctcrs = paste(yrcat, bay_segment, sep = ':'))
        prd <- visreg(mod, xvar = indvar, cond = condls, plot = F, scale = 'response')
        
        out <- prd$fit %>% 
          select(-bay_segment, -yrcat)
        
        return(out)
        
      })
    ) %>% 
    unnest(data)

  return(out)
  
}


# get model predictions gam, pinco
getprd_fun3 <- function(mod, indvar = c('sal', 'temp')){
  
  out <- unique(mod$model) %>% 
    group_by(yrcat) %>% 
    nest() %>% 
    mutate(
      data = purrr::pmap(list(yrcat), function(yrcat){
        
        condls <- list(yrcat = yrcat)
        prd <- visreg(mod, xvar = indvar, cond = condls, plot = F, scale = 'response')
        
        out <- prd$fit %>% 
          select(-yrcat)
        
        return(out)
        
      })
    ) %>% 
    unnest(data)
  
  return(out)
  
}