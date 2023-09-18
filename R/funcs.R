# function for formatting p-values in tables
p_ast <- function(x){
  
  sig_cats <- c('**', '*', '')
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

# get model predictions for 
getprd_fun <- function(modin){
  
  salfitotb <- visreg(modin, 'Sal', by = 'yrcat', cond = list(bay_segment = 'OTB'), scale = 'response', plot = F) %>% 
    .$fit %>%  
    mutate(
      bay_segment = 'OTB'
    )
  salfithb <- visreg(modin, 'Sal', by = 'yrcat', cond = list(bay_segment = 'HB'), scale = 'response', plot = F) %>% 
    .$fit %>% 
    mutate(
      bay_segment = 'HB'
    ) 
  salfitmtb <- visreg(modin, 'Sal', by = 'yrcat', cond = list(bay_segment = 'MTB'), scale = 'response', plot = F) %>% 
    .$fit %>% 
    mutate(
      bay_segment = 'MTB'
    )
  
  out <- bind_rows(salfitotb, salfithb, salfitmtb) %>% 
    mutate(
      bay_segment = factor(bay_segment, levels = segshr)
    )
  
  return(out)
  
}