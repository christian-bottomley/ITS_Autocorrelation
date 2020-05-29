#------------------------------------------------------------------------------#
# Data Example: Analysis of the impact of pneumoccocal vaccine introduction on # 
# the incidence of severe and very severe pneumonia in children 2-59 months in #
# Kilifi, Kenya. Analysis illustrates the use of differnt methods for adjusting#
# autocorrelation (Prais-Winsten method, ARMA modelling and Newey-West method) #                                                               #
#------------------------------------------------------------------------------#  

  library(tidyverse)
  library(lubridate)
  library(cowplot)
  library(prais)
  library(sandwich)
  library(forecast)

# Data management  
  
  dat <-
    read_csv("data_pneumo_ex.csv") %>%
    mutate(rate_10000 = 10000 * rate) %>%
    mutate(log2rate = log(rate_10000, base = 2)) %>%
    mutate(day = 14) %>%
    mutate(month = match(month, month.abb)) %>%
    mutate(date = make_date(year, month, day)) 
    dat$month <- factor(dat$month) 
    contrasts(dat$month) <- contr.sum(12)
    
# Fit models  
  
  lr <- lm(log2rate ~ vac + roll + studymo + month 
           + strike1 + strike2, data = dat)

  pw <- prais_winsten(log2rate ~ vac + roll + studymo + month 
                      + strike1 + strike2, data = dat) 
  
  aa <- auto.arima(dat$log2rate, 
                   xreg = model.matrix.lm(lr)[,-1], 
                   stepwise = T,
                   d = 0, 
                   D = 0, 
                   ic = 'aicc', 
                   seasonal = F
                   )
  
# Table of trend and vaccine impact estimates
  
  par <- c("studymo", "vac")
  
  est <- rbind(coef(lr)[par], coef(pw)[par], coef(aa)[par], coef(lr)[par])
  
  se <- rbind(sqrt(diag(vcov(lr))[par]),
              sqrt(diag(vcovHC.prais(pw))[par]),
              sqrt(diag(vcov(aa))[par]),
              sqrt(diag(NeweyWest(lr, 
                                  lag = NULL, prewhite = F, 
                                  adjust = T)[par, par]))
              )
  
  lb <- 100 * (1 - 2^(est + 1.96 * se))
  ub <- 100 * (1 - 2^(est - 1.96 * se))
  est <- 100 * (1 - 2^est)
  
  res_table <- cbind(
                    format(est[ ,"studymo"], digits = 2),
                    paste(
                     format(lb[ ,"studymo"], digits = 2),
                     format(ub[ ,"studymo"], digits = 2), sep =","),
                     format(est[ ,"vac"], digits = 2, nsmall = 1),
                     paste(
                      format(lb[ ,"vac"], digits = 2, nsmall = 1),
                      format(ub[ ,"vac"], digits = 2, nsmall = 1), sep =",")
                    )
  colnames(res_table) <- c("trend","ci trend", "vac", "ci vac")
  floor(bwNeweyWest(lr))
  view(res_table)              
  
# Add predicted values from linear reg to dataset
  
  dat <-
  dat %>% 
    mutate(pred_lr =   coef(lr)["(Intercept)"] 
         + coef(lr)["studymo"] * studymo
         + coef(lr)["vac"] * vac
    ) 
    
# Plot data and fitted model
 
  p1 <- ggplot(dat, aes(x = date)) +
    theme_bw() +
    geom_area(aes(y = log2rate),  fill = "gray", alpha = 0.75) +
    geom_line(aes(y = pred_lr), color = "black", size = 0.75,
              data = filter(dat, date < dmy("01-01-2011")) 
               ) +
    geom_line(aes(y = pred_lr), color = "black", size = 0.75,
              data = filter(dat, date > dmy("01-04-2011")) 
    ) +
    geom_rect(aes(xmin = dmy("01-01-2011"), 
                  xmax = dmy("31-03-2011"), 
                  ymin = 0, ymax = 4.8), fill = "#CC99CC", alpha = 0.5) + 
    geom_text(x = dmy("07-01-2011"), y = 5.0, label="roll-out") + 
    scale_x_date(breaks = c(dmy("01-01-2002"), dmy("01-01-2006"), 
                            dmy("01-01-2010"), dmy("01-01-2014"), 
                            dmy("01-01-2018")), 
                 date_minor_breaks = "1 year", date_labels = "%b %y") +
    ylab("monthly rate per 10,000 (log scale)") +
    scale_y_continuous(breaks = 1:5, labels = c("2", "4", "8", "16", "32"))

# Plot ACF
  
  ciline <- 1.96/sqrt(length(lr$residuals))
  acf_res <- acf(lr$residuals, plot = FALSE)
  acfdf <-  with(acf_res, data.frame(lag, acf)) 
 
  p2 <- ggplot(data = acfdf) +
    theme_bw() +
    geom_segment(aes(x = lag, xend = lag, y = acf, yend = 0)) +
    geom_hline(yintercept = ciline, linetype = "dashed") +
    geom_hline(yintercept = -ciline, linetype = "dashed") +
    scale_y_continuous(breaks = seq(-0.2, 1, 0.2)) 
    
# Combine plots
  
  p <- plot_grid(p1, p2, labels = "AUTO")
  save_plot("pneumo_ex_fig.tiff", p)
    
  
  
  
  