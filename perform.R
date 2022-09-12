#------------------------------------------------------------------------------#
# Evaluation of methods (Prais-Winsten, ARMA and Newey-West) in terms of bias, #
# mean square error and coverage using parameter estimates from simulated      #
# data. Performance measures and Monte Carlo error calculated using simsum     #
# in Stata (R package simsum does not use degrees of freedom to calculate      #
# power). Results outputted as a table (csv file) and as plots.                #
#------------------------------------------------------------------------------#

  library(tidyverse)
  library(gridExtra)
  library(gtable)
  library(grid)
  library(rsimsum)
  library(ggridges)
  library(RStata)

  options("RStata.StataPath" 
          = "/Applications/Stata/StataSE.app/Contents/MacOS/stata-se")
  
  options("RStata.StataVersion" = 17)
	
# Plotting function

  plotRes <- function(dat, y, colors) {
    ggplot(dat, aes(x = n, y = .data[[y]], color = method)) +
      geom_line(size = 1.7) +
      scale_color_manual(values = colors) +
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
    dplyr::select(-b1, -b2, -b3, -sd) %>%
    pivot_longer(-c("simno","scenario", "ma1", "ma2", "ma3", "n"), 
                         names_to = c(".value", "method"),
                         names_sep = "_"
                         ) %>% 
    mutate(
      df = ifelse(method %in% c("lr","pw"), n-3, df),
      df = ifelse(method %in% c("arma", "nw", "ma1", "ma2", "ma3"), 1000, df)
      ) 

# Use simsum to calculate performance of each method for each scenario
  
  stata_cmd <- 
    "simsum coef, true(-1) id(simno) methodvar(method) se(se) df(df) ///
     by(scenario n) mcse format(%4.3f) clear
    "
  dfres <- stata(stata_cmd, data.in = est_long, data.out = T, stata.echo = F)
  
# Reshape dataframe of results
  
  dfres<-
    dfres %>%   
    select(-statnum) %>% 
    rename_with(~ paste0("est_", .x), 
                -c("scenario","statcode", "n", ends_with("mcse"))) %>% 
    rename_with(~ paste0("mcse_", .x), ends_with("mcse"))
  
  colnames(dfres)<-gsub("coef","", colnames(dfres))
  colnames(dfres)<-gsub("_mcse","", colnames(dfres))
  
  dfres <-
    dfres %>% 
    filter(statcode == "bias" | statcode == "empse" 
           | statcode == "cover" | statcode == "power") %>% 
    pivot_longer(cols = starts_with(c("est", "mcse")),
                 names_to  = c(".value", "method"),
                 names_sep = "_"
                 ) %>% 
    mutate(statcode = factor(statcode)) %>% 
    pivot_wider(names_from = "statcode", 
                values_from = c("est", "mcse")) 
   
# Compute ACF for each scenario 
  
  ARMAacf2 <- function(mavalues) {ARMAacf(ma = mavalues)[2:4]} 
  
  mavalues <- 
    est %>% 
    as.data.frame() %>% 
    select(scenario, ma1, ma2, ma3) %>% 
    distinct()
  
  rvalues <- 
    t(apply(mavalues[,c("ma1","ma2","ma3")], 1, ARMAacf2)) %>% 
    as.data.frame() %>% 
    add_column(mavalues$scenario) 

  colnames(rvalues) <- c("r1", "r2", "r3", "scenario")

# Edit dataframe and export
  
  dfres <-
    dfres %>% 
    inner_join(rvalues, by="scenario") %>% 
    mutate(r1 = format(round(r1, 2)),
           r2 = format(round(r2, 2)),
           r3 = format(round(r3, 2))) %>%
    unite("autocorr", "r1", "r2", "r3", sep = "|")  %>%
    mutate(scenario = as.factor(scenario),
           method = factor(method,
                           levels = c("lr","pw","arma","nw",
                                      "ma1","ma2","ma3","ma3kr")),
           method = fct_recode(method, 
                        "Linear Reg" = "lr",
                        "Prais Winsten" = "pw",
                        "ARMA" = "arma",
                        "Newey West" = "nw",
                        "MA(1)" = "ma1",
                        "MA(2)" = "ma2",
                        "MA(3)" = "ma3",
                        "MA(3) K-R" = "ma3kr")
           ) %>% 
    rename(bias = est_bias) %>%
    rename(coverage = est_cover) %>%
    rename(power = est_power) %>% 
    rename(empse = est_empse) %>%
    mutate(coverage = coverage/100,
           mcse_cover = mcse_cover/100,
           power = power/100,
           mcse_power = mcse_power/100,
           mse = bias^2 + empse^2,
           mcse_mse = mcse_bias^2 + mcse_empse^2
    ) %>% 
    select(-c(empse, mcse_empse)) %>% 
    relocate(autocorr, .after = "scenario") %>% 
    relocate(mse, .after = "coverage") %>% 
    relocate(mcse_mse, .after = "mcse_cover") %>% 
    group_by(scenario) %>% 
    arrange(method,.by_group = T) %>% 
    ungroup()
    
  dfres %>% 
    mutate_if(is.numeric, ~sprintf(.,fmt = "%.3f")) %>% 
    write.csv(file = "table_perform.csv")
  
# Summarise results
  
  dfres %>%
    mutate(lag1corr = as.numeric(str_sub(autocorr, 1, 4))) %>%  
    mutate(lb_bias = bias - qnorm(0.975) * mcse_bias) %>%
    mutate(ub_bias = bias + qnorm(0.975) * mcse_bias) %>%
    #filter(n<=50) %>%
    #filter(lag1corr<0.3) %>%
    group_by(method) %>%
    summarise(
      bias_mean = mean(bias),
      coverage_mean = mean(coverage),
      rmse_mean = mean(sqrt(mse)),
      )
  
# Plot: bias, coverage, root mse	
    
  width <- 10
  height <- 1.2 * width

  colors = c("#D47DC5", "#B8643D", "#8CD9D2", "#37A447", "#D1D075",
             "#85C9D6","#000000", "#137735")
  
  dfres %>%
    filter(method %in% c("Linear Reg", "Prais Winsten", "ARMA", "MA(3)")) %>%
    plotRes("bias", colors = colors[c(1,2,3,7)]) %>%
    ggsave(file = "fig_bias.pdf", height = height, width = width)
  
  dfres %>%
    filter(method %in% c("Linear Reg", "Prais Winsten", "ARMA", "MA(3)")) %>%
    mutate(rmse = sqrt(mse)) %>%
    plotRes("rmse", colors = colors[c(1,2,3,7)]) %>%
    ggsave(file = "fig_rmse.pdf", height = height, width = width)
 
  dfres %>%
    filter(method %in% c("Linear Reg", "Prais Winsten", "ARMA", 
                         "Newey West", "MA(3)")) %>%
    plotRes("coverage", colors = colors[c(1,2,3,4,7)]) %>%
    ggsave(file = "fig_coverage.pdf", height = height, width = width)
  
  dfres %>%
    filter(method %in% c("MA(1)", "MA(2)", "MA(3)", "MA(3) K-R")) %>%
    plotRes("coverage", colors = colors[c(5,6,7,8)]) %>%
    ggsave(file = "fig_coverage_ma.pdf", height = height, width = width)
  
  dfres %>%
    filter(method %in% c("Linear Reg", "Prais Winsten", "ARMA", 
                         "Newey West", "MA(3) K-R")) %>% 
    filter(coverage > 0.85) %>% 
    filter(n <= 50) %>% 
    plotRes("power", colors = colors[c(1,2,3,4,8)]) %>%
    ggsave(file = "fig_power.pdf", height = height, width = width)
  
# Plot: distribution of SEs and coverage
  
  empse2_cover <-
    dfres %>%
    filter(scenario == 30) %>% 
    mutate(empse2 = mse - bias^2) %>% 
    select(method, empse2, coverage) 
  
  dftmp <-
  est %>% 
    as.data.frame() %>% 
    filter(scenario == 30) %>% 
    select(starts_with("se")) %>% 
    pivot_longer(cols = starts_with("se"),
                 names_to = "method",
                 values_to = "se",
                 names_prefix = "se_") %>%
    mutate(method = factor(method, 
                            levels = c("lr","pw","arma","nw",
                                       "ma1","ma2","ma3","ma3kr")),
           method = fct_recode(method, 
                        "Linear Reg" = "lr",
                        "Prais Winsten" = "pw",
                        "ARMA" = "arma",
                        "Newey West" = "nw",
                        "MA(1)" = "ma1",
                        "MA(2)" = "ma2",
                        "MA(3)" = "ma3",
                        "MA(3) K-R" = "ma3kr",
    )) %>% 
    left_join(empse2_cover, by = "method") %>% 
    mutate(method = fct_recode(method, 
                               "LR" = "Linear Reg",
                               "PW" = "Prais Winsten",
                               "ARMA" = "ARMA",
                               "NW" = "Newey West",
                               "MA(1)" = "MA(1)",
                               "MA(2)" = "MA(2)",
                               "MA(3)" = "MA(3)",
                               "MA(3)KR" = "MA(3) K-R"),
            modse2 = se^2) 
    
 p1 <-
   ggplot(dftmp, aes(x = method, y = modse2, color = method)) + 
    scale_color_manual(values = colors) +
    theme_bw(base_size = 20) +
    ylim(0,1.5) +
    labs(x = "", y = ~se^2) +
    geom_jitter(shape = 1, 
                position = position_jitter(width = 0.3, seed = 111), 
                show.legend = F,
                alpha = 0.5,
                size = 1) +
    stat_summary(fun = mean, 
                 geom = "crossbar", 
                 color = "darkorange1",
                 lwd = 0.5,
                 width = 0.65) +
    geom_point(aes(y=empse2), 
               color = "darkorange1", 
               shape = 4,
               stroke = 1,
               size = 5) 
  
 p2 <-   
    ggplot(dftmp[1:8,], aes(x = method, y = coverage,  fill = method)) + 
    scale_fill_manual(values = colors) +  
    theme_bw(base_size = 20) +
    geom_bar(stat = "identity", width = 0.5, show.legend = F) +
    labs(x = "", y = "coverage") +
    scale_y_continuous(limits = c(0.7, 1.0), 
                       breaks = c(0.7, 0.8, 0.9, 0.95, 1),
                       oob = scales::squish) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") 

 g1 <- ggplotGrob(p1)
 g2 <- ggplotGrob(p2)
 g <- rbind(g2, g1, size = "first")

 ggsave(g, file = "fig_se2_dist.pdf")

 
# Plot: distribution of MA model parameter estimates
  
  dftmp <-
   est %>% 
    as.data.frame() %>% 
    filter(scenario == 30) %>% 
    mutate(sig2_ma3 = sig_ma3^2,
           sig2_ma3kr = sig_ma3kr^2) %>% 
    select(ma1, ma2, ma3, sd,
           theta1_ma3, theta2_ma3, theta3_ma3, sig2_ma3,
           theta1_ma3kr, theta2_ma3kr, theta3_ma3kr, sig2_ma3kr) %>% 
    pivot_longer(cols = starts_with("theta") | starts_with("sig"),
                 names_to  = "param_method",
                 values_to = "est")  %>% 
    separate(param_method, into = c("param", "method"), 
             remove = F,
             sep="_") %>% 
    mutate(true = ma1,
           true = ifelse(param == "theta2", ma2, true),
           true = ifelse(param == "theta3", ma3, true),
           true = ifelse(param == "sig2", sd^2, true),
           param = factor(param),
           method = factor(method)
           ) %>% 
    filter(est < 2 & est > -1)

  ggplot(dftmp, aes(x = param, y = est, fill = method, color = method)) + 
    geom_jitter(position = position_jitterdodge(
                  jitter.width = 0.3,
                  dodge.width = 0.5,
                  seed = 111), 
                size = 1,
                pch = 21,
                alpha = 0.5,
                show.legend = T) +
    theme_bw(base_size = 24) +
    scale_fill_manual(
      labels = c("ML", "REML"), 
      values = colors[7:8]
      ) +
    scale_color_manual(
      labels = c("ML", "REML"), 
      values =  colors[7:8]
    ) +
    labs(x = "", y ="estimate") +
    stat_summary(fun = median, 
                geom = "crossbar", 
                position = position_dodge(),
                width = 0.4,
                lwd = 1,
                color = "darkorange1",
                show.legend = F) +
    geom_point(aes(y=true), 
               color = "darkorange1", 
               shape = 4,
               size = 8,
               stroke = 3,
               show.legend = F) +
    scale_x_discrete(labels = c('sig2' = expression(sigma^{2}),
                                'theta1'   = expression(theta[1]),
                                'theta2'   = expression(theta[2]),
                                'theta3'   = expression(theta[3])
                                )
                     ) +
   guides(colour = guide_legend(override.aes = list(size=6))) 
 
  ggsave(file = "fig_arma_param_dist.pdf", height = height, width = width)
 

