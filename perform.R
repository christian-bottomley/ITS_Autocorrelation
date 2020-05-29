#------------------------------------------------------------------------------#
# Evaluation of methods (Prais-Winsten, ARMA and Newey-West) in terms of bias, #
# mean square error and coverage using parameter estimates from simualated     #
# data. Performance measures and Monte Carlo error calculated using rsim       #
# package. Results outputted as a table and as a matrix of plots.              #
#------------------------------------------------------------------------------#

  library(tidyverse)
  library(wesanderson)
  library(rsimsum)
	
# Plotting function
  
  plotRes <- function(dat, y) {
    ggplot(dat, aes(x = n, y = .data[[y]], color = method)) +
      geom_line(size = 1.7) +
      scale_color_manual(values = wes_palette("IsleofDogs1")) +
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
    mutate(dfres = n-3) %>%
    pivot_longer(-c("simno","scenario", "r1", "r2", "r3", "n", "dfres"), 
                         names_to = c(".value", "method"),
                         names_pattern = "(.+)(.+)"
                         ) %>%
   rename(estimate = est) 
    
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
                               "Newey West" = "4"
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
    filter(method != "Newey West") %>%
    plotRes("bias") %>%
    ggsave(file = "plot_bias.tiff", height = height, width = width)
 
  dfRes %>%
    plotRes("coverage") %>%
    ggsave(file = "plot_coverage.tiff", height = height, width = width)
 
  dfRes %>%
    filter(method != "Newey West") %>%
    mutate(rmse = sqrt(mse)) %>%
    plotRes("rmse") %>%
    ggsave(file = "plot_rmse.tiff", height = height, width = width)
