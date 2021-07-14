#------------------------------------------------------------------------------#
# Simualtion study comparing different methods (Prais-Wimnsten, Newey West,    #
# ARMA) for handling autocorrelation in interrupted time series.               #
# Code generates simulated datasets and produces estimates and standard errors #
# for subsequent analysis of performance                                       #
#------------------------------------------------------------------------------#

  library(MASS)
  library(prais)
  library(sandwich)
  library(forecast)
  library(sfsmisc)

  set.seed(1066)

# Function to generate data
# Input: vector of parameter values (ma1, ma2, ma3, n, b1, b2, b3, sd)
# Output: dataframe (y=ITS, x=intervention indicator, t=time)

  genDat <- function(param) {
    ma <- param[1:3]
    n <- param[4]
    b <- param[5:7]
    sdv <- param[8]
    x <- c( rep( 0, n/2 ), rep( 1, n/2 ) )
    t <- 1:n
    mu <- b[1] + b[2] * x + b[3] * t
    e <- arima.sim(model = list(ma = ma), n = n, sd = sdv)
    y <- mu + e
    return(data.frame(y, x, t))
   }

# Function to suppress output (Hadley Wickham)

  quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }

# Function to fit models
# Input: dataframe (y, x, t)
# Output: dataframe of parameter estimates and standard errors

  fitMods <- function(df) {

  # Linear regression
    lr <- lm(y ~ x + t, data = df)
    est_lr <- coef(lr)["x"]
    se_lr <- sqrt(vcov(lr)["x", "x"])

  # Prais Winsten
   	pw <- quiet(prais_winsten(y ~ x + t, data = df))
    est_pw <- coef(pw)["x"]
    se_pw <- sqrt(vcovHC.prais(pw)["x", "x"])

  # ARMA
   	aa <- auto.arima(df$y,
   	                 xreg = cbind(df$x, df$t),
   	                 stepwise = T,
   	                 d = 0,
   	                 D = 0,
   	                 ic = 'aicc',
   	                 seasonal = F)
   	converge_aa <- ifelse(sum(diag(vcov(aa)) <= 0) > 0, 0, 1)
   	est_aa <- ifelse(converge_aa == 1,
   	                 coef(aa)["xreg1"],
   	                 NaN
   	                 )
   	se_aa <- ifelse(converge_aa == 1,
   	                sqrt(vcov(aa)["xreg1","xreg1"]),
   	                NaN
   	                )

  # Newey West
   	se_nw <- sqrt(NeweyWest(lr,
   	                        lag = NULL,
   	                        prewhite = F,
   	                        adjust = T)["x", "x"]
   	              )

  # MA1
   	ma1 <- arima(df$y, xreg = cbind(df$x, df$t), order = c(0, 0, 1))
   	est_ma1 <- coef(ma1)["cbind(df$x, df$t)1"]
   	se_ma1 <- sqrt(vcov(ma1)["cbind(df$x, df$t)1", "cbind(df$x, df$t)1"])

  # MA2
   	ma2 <- arima(df$y, xreg = cbind(df$x, df$t), order = c(0, 0, 2))
   	est_ma2 <- coef(ma2)["cbind(df$x, df$t)1"]
   	se_ma2 <- sqrt(vcov(ma2)["cbind(df$x, df$t)1", "cbind(df$x, df$t)1"])

  # MA3
   	ma3 <- arima(df$y, xreg = cbind(df$x, df$t), order = c(0, 0, 3))
   	est_ma3 <- coef(ma3)["cbind(df$x, df$t)1"]
   	se_ma3 <- sqrt(vcov(ma3)["cbind(df$x, df$t)1", "cbind(df$x, df$t)1"])

  # Combine and output estimates and standard errors
    est <- c(est_lr, est_pw, est_aa, est_lr, est_ma1, est_ma2, est_ma3)
    se <- c(se_lr, se_pw, se_aa, se_nw, se_ma1, se_ma2, se_ma3)
    cbind(est, se)
  }

# Matrix of parameter values (ma1-ma3, n, b1-b3, sd)
# Each row represents a different scenario

  nma <- 20
  ma1values <- sort(runif(nma, min = 0, max = 1))
  ma2values <- runif(nma, min = 0, max = ma1values)
  ma3values <- runif(nma, min = 0, max = ma2values)
  nvalues  <- c(20, 50, 100, 300)
  nn <- length(nvalues)
  b1 <- 4
  b2 <- -1
  b3 <- -1
  sdv <- 1

  param <- cbind(
    rep(ma1values, each = nn),
    rep(ma2values, each = nn),
    rep(ma3values, each = nn),
    rep(nvalues, nma),
    rep(b1, nma * nn),
    rep(b2, nma * nn),
    rep(b3, nma * nn),
    rep(sdv, nma * nn)
  )

  colnames(param) <- c("ma1", "ma2", "ma3","n", "b1", "b2", "b3", "sd")
  param[ , "b3"] <-  param[ , "b3"] / param[ , "n"]
  
# Models fitted to n=2000 simulated datasets for each scenario
# Parameter values and standard errors stored in matrix "est"

  nsim <- 2000
  est <- matrix(0, nsim * nrow(param), 2 + ncol(param) + 14)
  colnames(est) <- c(
                     "simno", "scenario", "ma1", "ma2", "ma3",
                     "n", "b1", "b2", "b3", "sd",
                     "est1", "est2", "est3", "est4", "est5", "est6", "est7", 
                     "se1", "se2", "se3", "se4", "se5","se6", "se7"
                     )
  
  for(i in 1:nrow(param)) {
    start_time <- Sys.time() 
    simRes <- replicate(nsim, fitMods(
                                      genDat(param = param[i, ])
                                      )
                        )
    rows <- ((i - 1) * nsim + 1):(i * nsim)
    est[rows, 1] <- 1:nsim
    est[rows, 2] <- i
    est[rows, 3:10] <- matrix(param[i, ], ncol = 8, nrow = nsim, byrow = T)
    est[rows, 11:17] <-  t(simRes[ , 1, ])
    est[rows, 18:24] <-  t(simRes[ , 2, ])
   	end_time <- Sys.time() 
   	print(i)
   	print(end_time - start_time)
   }   
     
 # Save results

  save(est, file="est.rda")
  
