#------------------------------------------------------------------------------#
# Simulation study comparing different methods (Prais-Winsten, Newey West,     #
# ARMA) for handling autocorrelation in interrupted time series.               #
# Code generates simulated datasets and produces estimates and standard errors #
# for subsequent analysis of performance                                       #
#------------------------------------------------------------------------------#

  library(MASS)
  library(prais)
  library(sandwich)
  library(forecast)
  library(sfsmisc)
  library(RStata)
  
  options("RStata.StataPath" 
          = "/Applications/Stata/StataSE.app/Contents/MacOS/stata-se")
  
  options("RStata.StataVersion" = 17)

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
    coef_lr <- coef(lr)["x"]
    se_lr <- sqrt(vcov(lr)["x", "x"])

  # Prais Winsten
   	pw <- quiet(prais_winsten(y ~ x + t, data = df))
    coef_pw <- coef(pw)["x"]
    se_pw <- sqrt(vcovHC.prais(pw)["x", "x"])

  # ARMA
   	arma <- auto.arima(df$y,
   	                 xreg = cbind(df$x, df$t),
   	                 stepwise = T,
   	                 d = 0,
   	                 D = 0,
   	                 ic = 'aicc',
   	                 seasonal = F)
   	converge_arma <- ifelse(sum(diag(vcov(arma)) <= 0) > 0, 0, 1)
   	coef_arma <- ifelse(converge_arma == 1,
   	                 coef(arma)["xreg1"],
   	                 NaN
   	                 )
   	se_arma <- ifelse(converge_arma == 1,
   	                sqrt(vcov(arma)["xreg1","xreg1"]),
   	                NaN
   	                )

  # Newey West
   	se_nw <- sqrt(NeweyWest(lr,
   	                        lag = NULL,
   	                        prewhite = F,
   	                        adjust = T)["x", "x"]
   	              )

  # MA1 maximum likelihood
   	ma1 <- arima(df$y, xreg = cbind(df$x, df$t), order = c(0, 0, 1))
   	coef_ma1 <- coef(ma1)["cbind(df$x, df$t)1"]
   	se_ma1 <- sqrt(vcov(ma1)["cbind(df$x, df$t)1", "cbind(df$x, df$t)1"])
   	sig_ma1 <-sqrt(ma1$sigma2)
   	theta1_ma1 <- coef(ma1)["ma1"]
   	
   
  # MA2 maximum likelihood
   	ma2 <- arima(df$y, xreg = cbind(df$x, df$t), order = c(0, 0, 2))
   	coef_ma2 <- coef(ma2)["cbind(df$x, df$t)1"]
   	se_ma2 <- sqrt(vcov(ma2)["cbind(df$x, df$t)1", "cbind(df$x, df$t)1"])
   	sig_ma2 <-sqrt(ma2$sigma2)
   	theta1_ma2 <- coef(ma2)["ma1"]
   	theta2_ma2 <- coef(ma2)["ma2"]

  # MA3 maximum likelihood
   	ma3 <- arima(df$y, xreg = cbind(df$x, df$t), order = c(0, 0, 3))
   	coef_ma3 <- coef(ma3)["cbind(df$x, df$t)1"]
   	se_ma3 <- sqrt(vcov(ma3)["cbind(df$x, df$t)1", "cbind(df$x, df$t)1"])
   	sig_ma3 <-sqrt(ma3$sigma2)
   	theta1_ma3 <- coef(ma3)["ma1"]
   	theta2_ma3 <- coef(ma3)["ma2"]
   	theta3_ma3 <- coef(ma3)["ma3"]
   
  # MA3 REML + Kenward-Roger correction (only for n<=50)
   	if(nrow(df) <= 50) {
   	  stata_cmd <- 
   	  "mixed y  x  t, residuals(ma 3, t(t)) reml dfmethod(kroger) 
      matrix A = r(table) 
      matrix list A
      gen coef = A[1,1] 
      gen se = A[2,1]
      gen df = A[7,1]
      gen logsigma = A[1,4]
      gen theta1 = A[1,5]
      gen theta2 = A[1,6]
      gen theta3 = A[1,7]
      "
   	  ma3kr <- stata(stata_cmd, data.in = df, data.out = T, stata.echo = F)
   	  coef_ma3kr <- ma3kr[1,"coef"]
   	  se_ma3kr <- ma3kr[1,"se"]
    	df_ma3kr <- ma3kr[1,"df"]
    	sig_ma3kr <- exp(ma3kr[1,"logsigma"])
    	theta1_ma3kr <- ma3kr[1,"theta1"]
    	theta2_ma3kr <- ma3kr[1,"theta2"]
    	theta3_ma3kr <- ma3kr[1,"theta3"]
    	} else {
   	  coef_ma3kr <- NA
   	  se_ma3kr <- NA
   	  df_ma3kr <- NA
   	  sig_ma3kr <- NA
   	  theta1_ma3kr <- NA
   	  theta2_ma3kr <- NA
   	  theta3_ma3kr <- NA
      }
   	
  # Output estimates and standard errors
    coef <- c(coef_lr, coef_pw, coef_arma, coef_lr, 
              coef_ma1, coef_ma2, coef_ma3, coef_ma3kr)
    se <- c(se_lr, se_pw, se_arma, se_nw, 
            se_ma1, se_ma2, se_ma3, se_ma3kr)
    sig <- c(rep(NA, 4), sig_ma1, sig_ma2, sig_ma3, sig_ma3kr)
    df_ma3kr <- c(rep(NA, 7), df_ma3kr)
    theta1 <- rep(NA, 8)
    theta1[5:8] <- c(theta1_ma1, theta1_ma2, theta1_ma3, theta1_ma3kr)
    theta2 <- rep(NA, 8)
    theta2[6:8] <- c(theta2_ma2, theta2_ma3, theta2_ma3kr)
    theta3 <- rep(NA, 8)
    theta3[7:8] <- c(theta3_ma3, theta3_ma3kr)
    
    out <- cbind(coef, se, sig, theta1, theta2, theta3, df_ma3kr)
    rownames(out) <- c("lr", "pw", "arma", "nw", "ma1", "ma2", "ma3", "ma3kr")
    out  
  }

# Matrix of parameter values (ma1-ma3, n, b1-b3, sd)
# Each row represents a different scenario

  nma <- 20
  ma1values <- sort(runif(nma, min = 0, max = 1))
  ma2values <- runif(nma, min = 0, max = ma1values)
  ma3values <- runif(nma, min = 0, max = ma2values)
  ma1values[1] <- 0
  ma2values[1] <- 0
  ma3values[1] <- 0
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
# Estimates and standard errors stored in matrix "est"

  nsim <- 2000
  columnnames <-
    c("simno", "scenario", 
      "ma1", "ma2", "ma3", "n", "b1", "b2", "b3", "sd",
      "coef_lr", "coef_pw", "coef_arma", "coef_nw", 
      "coef_ma1", "coef_ma2", "coef_ma3", "coef_ma3kr",
      "se_lr", "se_pw", "se_arma", "se_nw", 
      "se_ma1", "se_ma2", "se_ma3", "se_ma3kr",
      "sig_ma1", "sig_ma2", "sig_ma3", "sig_ma3kr",
      "theta1_ma1", "theta1_ma2",  "theta1_ma3", "theta1_ma3kr", 
      "theta2_ma2", "theta2_ma3", "theta2_ma3kr",
      "theta3_ma3", "theta3_ma3kr",
      "df_ma3kr")
  
  est <- matrix(0, nsim * nrow(param), length(columnnames))
  colnames(est) <- columnnames
  
  for(i in 1:nrow(param)) {
    
    start_time <- Sys.time() 
    
    simRes <- replicate(nsim, 
                        fitMods(
                                genDat(param = param[i, ])
                            )
                        )
    
    rows <- ((i - 1) * nsim + 1):(i * nsim)
    est[rows, 1] <- 1:nsim
    est[rows, 2] <- i
    est[rows, 3:10] <- matrix(param[i, ], ncol = 8, nrow = nsim, byrow = T)
    est[rows, 11:18] <-  t(simRes[ , "coef", ])
    est[rows, 19:26] <-  t(simRes[ , "se", ])
    est[rows, 27:30] <-  t(simRes[ c("ma1", "ma2", "ma3", "ma3kr"), "sig", ])
    est[rows, 31:34] <-  t(simRes[ c("ma1", "ma2", "ma3", "ma3kr"), "theta1", ])
    est[rows, 35:37] <-  t(simRes[ c("ma2", "ma3", "ma3kr"), "theta2", ])
    est[rows, 38:39]    <-  simRes[c("ma3", "ma3kr"), "theta3", ]
    est[rows, 40]    <-  simRes["ma3kr", "df_ma3kr", ]
  
   	end_time <- Sys.time() 
   	print(i)
   	print(end_time - start_time)
   
   	}   
     
 # Save results

  save(est, file="est.rda")
  
  
