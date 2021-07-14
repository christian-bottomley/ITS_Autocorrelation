#------------------------------------------------------------------------------#
# Evaluation of methods (Prais-Winsten, ARMA and Newey-West) in terms of bias, #
# mean square error and coverage using parameter estimates from simulated      #
# data. Performance measures and Monte Carlo error calculated using rsim       #
# package. Results outputted as a table (csv file) and as a matrix of plots.   #
#------------------------------------------------------------------------------#

  library(tidyverse)
  library(rsimsum)
	
# Plotting function
  
  plotRes <- function(dat, y) {
    ggplot(dat, aes(x = n, y = .data[[y]], color = method)) +
      geom_line(size = 1.7) +
      scale_color_manual(values = c("Linear Reg" = "#D47DC5",
                                    "Prais Winsten" = "#B8643D",
                                    "ARMA" = "#8CD9D2",
                                    "Newey West" = "#37A447",
                                    "MA(1)" = "#D1D075",
                                    "MA(2)" = "#85C9D6",
                                    "MA(3)" = "#000000")) +
      theme_bw(base_size = 20) +
      theme(legend.position = "top") +
      labs(y = y) +
      facet_wrap(~autocorr)
  }	

# Load data and reshape to long format
  
  load("est.rda")
  est_long <- 
    est %>%
    as.data.frame() %>%
    select(-b1, -b2, -b3, -sd) %>%
    pivot_longer(-c("simno","scenario", "ma1", "ma2", "ma3", "n"), 
                         names_to = c(".value", "method"),
                         names_pattern = "(.+)(.+)"
                         ) %>%
   rename(estimate = est) %>%
   mutate(dfres = ifelse(method %in% 1:2, n-3, 1000)) 
  
# Add ACF
  
  ARMAacf2 <- function(mavalues) {ARMAacf(ma = mavalues)[2:4]} 
  mavalues <- est_long[, c("ma1", "ma2", "ma3")]
  rvalues <- t(apply(mavalues, 1, ARMAacf2))
  colnames(rvalues) <- c("r1", "r2", "r3")
  est_long <- cbind(est_long, rvalues)
  
# Use simsum to calculate performance of each method for each scenario
    
   prfm <- 
     simsum(est_long, estvarname = "estimate", se = "se", df = "dfres",
            true = -1, methodvar = "method", by = "scenario", ref = "1"
           )  
  
# Dataframe of results
   
  dfRes <- 
    prfm$sum %>%
    filter(stat == "bias" | stat == "mse" | stat == "cover")  %>% 
    pivot_wider(names_from = stat, values_from = c("est", "mcse")) %>%
    select("scenario", everything())
  
  param <- 
    est_long[, c( "scenario", "method", "r1", "r2", "r3", "n")] %>%
    distinct() %>%
    mutate(r1 = format(round(r1, 2))) %>%
    mutate(r2 = format(round(r2, 2))) %>%
    mutate(r3 = format(round(r3, 2))) %>%
    unite("autocorr", "r1", "r2", "r3", sep = "|")  %>%
    mutate(scenario = as.factor(scenario)) %>%
    mutate(method = as.factor(method))
  
  dfRes <- 
    param %>%
    inner_join(dfRes, by = c("scenario", "method")) %>%
    mutate(method = fct_recode(method, 
                               "Linear Reg" = "1",
                               "Prais Winsten" = "2",
                               "ARMA" = "3",
                               "Newey West" = "4",
                               "MA(1)" = "5",
                               "MA(2)" = "6",
                               "MA(3)" = "7",
                               )
           ) %>%
    rename(bias = est_bias) %>%
    rename(coverage = est_cover) %>%
    rename(mse = est_mse) %>%
    mutate_if(is.numeric, ~round(.,3)) 

  write.csv(dfRes, file = "tabl_perform.csv")
  
# Summarise results
  
  dfRes %>%
    mutate(lag1corr = as.numeric(str_sub(autocorr, 1, 4))) %>%  
    mutate(lb_bias = bias - qnorm(0.975) * mcse_bias) %>%
    mutate(ub_bias = bias + qnorm(0.975) * mcse_bias) %>%
    #filter(n<50) %>%
    #filter(lag1corr>=0.6) %>%
    group_by(method) %>%
    summarise(
      bias_mean = mean(bias),
      coverage_mean = mean(coverage),
      rmse_mean = mean(sqrt(mse)),
      )
  
# Graphs: bias, coverage, root mse	
    
  width <- 10
  height <- 1.2 * width
  
  dfRes %>%
    filter(!method %in% c("MA(1)", "MA(2)", "Newey West")) %>%
    plotRes("bias") %>%
    ggsave(file = "fig_bias.pdf", height = height, width = width)
  
  dfRes %>%
    filter(!method %in% c("MA(1)", "MA(2)", "Newey West")) %>%
    mutate(rmse = sqrt(mse)) %>%
    plotRes("rmse") %>%
    ggsave(file = "fig_rmse.pdf", height = height, width = width)
 
  dfRes %>%
    filter(!method %in% c("MA(1)", "MA(2)")) %>%
    plotRes("coverage") %>%
    ggsave(file = "fig_coverage.pdf", height = height, width = width)
  
  dfRes %>%
    filter(method %in% c("MA(1)", "MA(2)", "MA(3)")) %>%
    plotRes("coverage") %>%
    ggsave(file = "fig_coverage_ma.pdf", height = height, width = width)
    

    
  