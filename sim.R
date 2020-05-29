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

# Function to generate correlation matrix                                      
# Input:                                                                       
# 1) ITS length (n)                                                            
# 2) vector of lag1-lag3 correlations (r)                                      
# Output: n x n correlation matrix                                             

  genCorrMat <- function(n, r) {
    rmat <- diag(1, n)
    for(i in 1:n) {
	    for(j in 1:n) {
	      if (abs(i-j) == 1) {
	        rmat[i,j] <- r[1]
	      }
	      if (abs(i-j) == 2) {
	        rmat[i,j] <- r[2]
	      }
	      if (abs(i-j) == 3) {
	        rmat[i,j] <- r[3]
	       }	     
	     }
    }
    rmat <- posdefify(rmat) 
	  return(rmat)
  }
  
# Function to generate data 
# Input: vector of parameter values (r1, r2, r3, n, b1, b2, b3, sd)
# Output: dataframe (y=ITS, x=intervention indicator, t=time)  
         
  genDat <- function(param) {
    r <- param[1:3]
    n <- param[4]
    b <- param[5:7]
    sdv <- param[8]
    x <- c( rep( 0, n/2 ), rep( 1, n/2 ) ) 
    t <- 1:n    
    mu <- b[1] + b[2] * x + b[3] * t 
    rmat <- genCorrMat(n, r)
    e <- mvrnorm(1, rep(0, n), sdv^2*rmat )
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
  
  # Combine and output estimates and standard errors 
    est <- c(est_lr, est_pw, est_aa, est_lr) 		
    se <- c(se_lr, se_pw, se_aa, se_nw) 
    cbind(est, se)  	
  }
  
# Matrix of parameter values (r1-r3, n, b1-b3, sd)
# Each row represents a different scenario
  
  nr <- 20
  r1values <- sort(runif(nr, min = 0, max = 0.8))
  r2values <- runif(nr, min = 0, max = r1values)
  r3values <- runif(nr, min = 0, max = r2values)
  nvalues  <- c(20, 50, 100, 300) 
  nn <- length(nvalues)
  b1 <- 4 
  b2 <- -1
  b3 <- -1
  sdv <- 1
  
  param <- cbind(
    rep(r1values, each = nn),
    rep(r2values, each = nn),
    rep(r3values, each = nn),
    rep(nvalues, nr),
    rep(b1, nr * nn),
    rep(b2, nr * nn),
    rep(b3, nr * nn),
    rep(sdv, nr * nn)
  )
  
  colnames(param) <- c("r1", "r2", "r3", "n", "b1", "b2", "b3", "sd")
  param[ , "b3"] <-  param[ , "b3"] / param[ , "n"]     
    
# Models fitted to n=2000 simulated datasets for each scenario
# Parameters and standard errors stored in matrix "est"

  nsim <- 2000
  est <- matrix(0, nsim * nrow(param), ncol(param) + 2 + 8)
  colnames(est) <- c(
                     "simno", "scenario",
                     "r1", "r2", "r3", "n", "b1", "b2", "b3", "sd",
                     "est1", "est2", "est3", "est4",
                     "se1", "se2", "se3", "se4"
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
    est[rows, 11:14] <-  t(simRes[ , 1, ])
    est[rows, 15:18] <-  t(simRes[ , 2, ])
   	end_time <- Sys.time() 
   	print(i)
   	print(end_time - start_time)
   }   
     
 # Save results

  save(est, file="est.rda")
  
  
