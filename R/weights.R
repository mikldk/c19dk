#' @importFrom magrittr "%>%" 
#' @importFrom dplyr filter select inner_join mutate
calc_obs_pred <- function(samples, data_indlagte) {
  obs_int <- data_indlagte %>% 
    filter(inc_deviation) %>% 
    select(Dato, Intensiv)
  res_int <- samples$sim_indlagt_int %>% 
    inner_join(obs_int, by = "Dato")
  res_int
}

#' @importFrom magrittr "%>%" 
#' @importFrom dplyr mutate
normalise_residuals <- function(d, min_weight = 0.1) {
  if (min_weight < 1e-10 | min_weight > (1-1e-10)) {
    stop("min_weight must be between 1e-10 and 1-1e-10")
  }
  
  x <- 1 - min_weight
  
  # residuals high: bad parameters
  # weights low: bad parameters
  # weights scaled from 0.1 to 1 such that worst model gets 0.1 and best 1.
  d %>% 
    mutate(weight = (1-x*(weight_raw-min(weight_raw))/(max(weight_raw)-min(weight_raw))))
}


#' @importFrom dplyr group_by summarise
#' @export
weights_mse <- function(samples, data_indlagte, min_weight = 0.1) {
  # (y - yhat)^2
  
  res_int <- calc_obs_pred(samples, data_indlagte) %>% 
    mutate(res = Antal - Intensiv)

  res_int_w <- res_int %>% 
    group_by(MC) %>% 
    summarise(weight_raw = sum(res^2)) %>% 
    normalise_residuals(min_weight = min_weight)
  
  res_int_w
}

#' @importFrom dplyr group_by summarise
#' @export
weights_msre <- function(samples, data_indlagte, min_weight = 0.1) {
  # (log(y+1) - log(yhat+1))^2
  
  res_int <- calc_obs_pred(samples, data_indlagte) %>% 
    mutate(res = log(Antal + 1) - log(Intensiv + 1))
  
  res_int_w <- res_int %>% 
    group_by(MC) %>% 
    summarise(weight_raw = sum(res^2)) %>% 
    normalise_residuals(min_weight = min_weight)
  
  res_int_w
}

