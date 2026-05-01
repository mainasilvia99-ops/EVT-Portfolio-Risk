###############################################
# Extreme Financial Risk Meaurement EVT MODELS 
###############################################

par(mar = c(4,4,2,1))

#Comparing tails of distributions

compare.densities <- function(x.range = c(0, 5), y.range = c(0, 3/2)){
  f1 <- function(x){dnorm(x, mean = 1, sd = sqrt(2))}
  f2 <- function(x, alpha = 4, lambda = 3){alpha * lambda^alpha/(lambda+x)^(alpha + 1)}
  
  curve(f1, from = x.range[1], to = x.range[2], lwd = 3, lty = 1, col = 'brown2', 
        ylim = y.range, ylab = 'densities')
  curve(f2, from = x.range[1], to = x.range[2], lwd = 3, lty = 2, 
        col = 'cornflowerblue', add = TRUE )
  
  legend("topright", ncol = 1, cex = 2, legend = c("f1", "f2"),
         col = c('brown2', 'cornflowerblue'), lwd = 4, lty = c(1, 2))
}


compare.densities()

#zoom on tails
compare.densities(c(3,6),c(0,0.12))


#Creating functions for the GEV density
#case "xi" is non-zero

dGEV <- function(x, m, s, xi){
  z <- 1 + xi*(x - m)/s
  ifelse(z <= 0, 0,
         (1/s) * z^(-1/xi - 1) * exp(-z^(-1/xi)))
}

#case "xi" is zero
dGumbel <- function(x, m, s){
  (1/s) * exp(-(x - m)/s) * exp(-exp(-(x - m)/s))
}

pGEV <- function(x, m, s, xi){
  z <- 1 + xi*(x - m)/s
  ifelse(z <= 0, 0, exp(-z^(-1/xi)))
}

#Generating sample of maximum values -Block maxima generator

generate.max <- function(B = 50, n = 10, rdist, param = c(0,1), a = 1, b = 0){
  my.max <- numeric(B)
  
  for(i in 1:B){
    sample <- rdist(n, param[1], param[2])
    my.max[i] <- (max(sample) - b)/a
  }
  
  hist(my.max, breaks = "FD", prob = TRUE,
       col = "springgreen2",
       main = paste(" Maxima (n = ", n,") "),
       xlab = "")
  
  return(invisible(my.max))
}


#B is the no of replicates of the simulations
#n is the sample size of the extreme values from the simulations
#Because we want a sample of simulated maxima we replicate it B times to see the distribution of the sample n

###########################################
# (NEGATIVE) WEIBULL (Uniform → ξ < 0)
###########################################

#We try with F ~ Uniform is in MDA of negative weibull with psi < 0 n=5
generate.max(B = 50000, n = 5, rdist = runif, param = c(0,1), a = 1/5, b = 1)

curve(dGEV(x, -1, 1, -1), type = "l", lwd = 2, lty = 2,
      col = "brown3", n = 500, add = TRUE)

#We try with F ~ Uniform is in MDA of negative weibull with psi < 0 n=10
generate.max(B = 50000, n = 10, rdist = runif, param = c(0,1), a = 1/10, b = 1)

curve(dGEV(x, -1, 1, -1), type = "l", lwd = 2, lty = 2,
      col = "brown3", n = 500, add = TRUE)

#We try with F ~ Uniform is in MDA of negative weibull with psi < 0 n=100
generate.max(B = 50000, n = 100, rdist = runif, param = c(0,1), a = 1/100, b = 1)

curve(dGEV(x, -1, 1, -1), type = "l", lwd = 2, lty = 2,
      col = "brown3", n = 500, add = TRUE)

#We try with F ~ Uniform is in MDA of negative weibull with psi < 0 n=500
generate.max(B = 50000, n = 500, rdist = runif, param = c(0,1), a = 1/500, b = 1)

curve(dGEV(x, -1, 1, -1), type = "l", lwd = 2, lty = 2,
      col = "brown3", n = 500, add = TRUE)

###########################################
# GUMBEL (Exponential → ξ = 0)
###########################################

#We try with F ~ Exponential, n = 5 is in MDA of to gumbel with psi = 0
generate.max(B = 50000, n = 5, rdist = rgamma, param = c(1, 1), a = 1, b = log(5)) #n=5,10,100,500

curve(dGumbel(x, 0, 1), type = "l", lwd = 4, lty = 2,
      col = "brown3", n = 500, add = TRUE)

#We try with F ~ Exponential, n = 10 is in MDA of to gumbel with psi = 0
generate.max(B = 50000, n = 10, rdist = rgamma, param = c(1, 1), a = 1, b = log(10)) #n=5,10,100,500

curve(dGumbel(x, 0, 1), type = "l", lwd = 3, lty = 2,
      col = "brown3", n = 500, add = TRUE)

#We try with F ~ Exponential, n = 100 is in MDA of to gumbel with psi = 0
generate.max(B = 50000, n = 100, rdist = rgamma, param = c(1, 1), a = 1, b = log(100)) #n=5,10,100,500

curve(dGumbel(x, 0, 1), type = "l", lwd = 3, lty = 2,
      col = "brown3", n = 500, add = TRUE)

#We try with F ~ Exponential, n = 500 is in MDA of to gumbel with psi = 0
generate.max(B = 50000, n = 500, rdist = rgamma, param = c(1, 1), a = 1, b = log(500)) #n=5,10,100,500

curve(dGumbel(x, 0, 1), type = "l", lwd = 3, lty = 2,
      col = "brown3", n = 500, add = TRUE)

#We try with F ~ Exponential, n = 1000 is in MDA of to gumbel with psi = 0
generate.max(B = 50000, n = 1000, rdist = rgamma, param = c(1, 1), a = 1, b = log(1000)) #n=5,10,100,500

curve(dGumbel(x, 0, 1), type = "l", lwd = 3, lty = 2,
      col = "brown3", n = 500, add = TRUE)

#We try with F ~ Exponential, n = 5000 is in MDA of to gumbel with psi = 0
generate.max(B = 50000, n = 5000, rdist = rgamma, param = c(1, 1), a = 1, b = log(5000)) #n=5,10,100,500

curve(dGumbel(x, 0, 1), type = "l", lwd = 2, lty = 2,
      col = "brown3", n = 500, add = TRUE)

###########################################
# PARETO (Fréchet → ξ > 0)
###########################################

#We first create a Pareto 'random generator'

#pareto random generator function

rpareto <- function(n, alpha, lambda){
  u <- runif(n)
  lambda*((1 - u)^(-1/alpha)-1)
}

alpha <- 5
xi <- 1/alpha


# F ~ Pareto, n = 10 is in MDA of to frechet with psi > 0

generate.max(B = 50000, n = 10, 
             rdist = rpareto, 
             param = c(alpha, 1), 
             a = 10^(1/alpha) - 1, 
             b = 0) #n&a=10,100,1000,3000

curve(dGEV(x, 1, 1/5, xi), type = "l", lwd = 3, lty = 2,
      col = "brown3", n = 500, add = TRUE)

# F ~ Pareto, n = 100 is in MDA of to frechet with psi > 0

generate.max(B = 50000, n = 100, rdist = rpareto, param = c(alpha, 1), a = 100^(1/alpha) - 1, b = 0) #n&a=10,100,1000,3000

curve(dGEV(x, 1, 1/5, xi), type = "l", lwd = 3, lty = 2,
      col = "brown3", n = 500, add = TRUE)


# F ~ Pareto, n = 1000 is in MDA of to frechet with psi > 0

generate.max(B = 50000, n = 1000, rdist = rpareto, param = c(alpha, 1), a = 1000^(1/alpha) - 1, b = 0) #n&a=10,100,1000,3000

curve(dGEV(x, 1, 1/5, xi), type = "l", lwd = 3, lty = 2,
      col = "brown3", n = 500, add = TRUE)


# F ~ Pareto, n = 5000 is in MDA of to frechet with psi > 0

generate.max(B = 50000, n = 5000, rdist = rpareto, param = c(alpha, 1), a = 5000^(1/alpha) - 1, b = 0) #n&a=10,100,1000,3000

curve(dGEV(x, 1, 1/5, xi), type = "l", lwd = 3, lty = 2,
      col = "brown3", n = 500, add = TRUE)


#Convergence is much slower for the Pareto...

####################################
#FITTING THE GEV TO REAL MARKET DATA
####################################


#LOAD DATA + PREP DATA FROM EXCEL
 library(dplyr)

return.data <- stock_returns_data
return.data$date <- as.Date(return.data$date)

head(return.data)

portfolio <- return.data %>%
  group_by(date) %>%
  summarise(portfolio_return = mean(returns, na.rm = TRUE))

#Extract MONTHLY worst loss (BLOCK MAXIMA)
library(lubridate)

min_returns <- portfolio %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(min_return = min(portfolio_return, na.rm = TRUE))

head(min_returns)

unique(substr(min_returns$month, 1, 4))

range(return.data$date)

#Convert to EVT losses
losses <- -min_returns$min_return

summary(losses)

hist(losses, breaks = "FD", 
     col = 'lightgoldenrod', prob = TRUE,
     main = "Histogram of the largest monthly loss", 
     xlab = '')

#The shape suggests a GEV with psi >= 0 i.e. Gumbel or Frechet

#We now fit the GEV to our data, by Maximum Likelihood MLE

#negative log-likelihood of our data (function of 3 parameters )
nll <- function(p){
  -sum(log(dGEV(losses, p[1],p[2],p[3])))
  }

#initial parameters for optimization
iv <- c(0.0, 0.01, 0.1)

#minimize f
fit <- nlm(nll, iv)

#display parameter estimates
(par <- fit$estimate)


#ASSESS THE FIT WITH DIAGNOSTIC PLOTS
#--- plot 1 ---
#Histogram vs fitted density

hist(losses, breaks = "FD",
     col = 'lightgoldenrod', 
     prob = TRUE,
     main = "Histogram of the largest monthly loss", 
     xlab = '', ylim = c(0, 60))

curve(dGEV(x, m = par[1], s = par[2], xi = par[3]), from =min(losses), to = max(losses), lwd = 3, col = 'brown2', add = TRUE)

#--- plot 2 ---
#COMPARING CDFs (Empirical vs fitted GEV)

# CDF comparison
plot(ecdf(losses), main = "Empirical vs GEV CDF",
     ylab = "CDF", lwd = 2)

curve(pGEV(x, par[1], par[2], par[3]),
      from = min(losses), to = max(losses),
      add = TRUE, col = "brown2", lwd = 2)



#--- plot 3 ---
#QQ plot
#careful qq plot uses package EnvStats with shape parameter as -psi
install.packages("EnvStats")

library(EnvStats)


EnvStats::qqPlot(losses,
                 distribution = "gevd",
                 param.list = list(location = par[1],
                                   scale = par[2],
                                   shape = -par[3]),
                 add.line = TRUE,
                 line.col = "brown3")

#INCREASE DATA FREQUENCY
#Extract WEEKLY worst loss (BLOCK MAXIMA)

library(dplyr)
library(lubridate)

min_returns_weekly <- portfolio %>%
  mutate(week = floor_date(date, "week")) %>%
  group_by(week) %>%
  summarise(min_return = min(portfolio_return, na.rm = TRUE))

class(portfolio$date)

#Convert to EVT losses
weekly_losses <- -min_returns$min_return

summary(weekly_losses)

hist(weekly_losses, breaks = "FD", 
     col = 'lightgoldenrod', prob = TRUE,
     main = "Histogram of the largest weekly loss", 
     xlab = '')

#The shape suggests a GEV with psi >= 0 i.e. Gumbel or Frechet

#We now fit the GEV to our data, by Maximum Likelihood MLE

#negative log-likelihood of our data (function of 3 parameters )
iv <- c(0, 1, 0.1)

nll_safe <- function(p){
  m <- p[1]
  s <- p[2]
  xi <- p[3]
  
  if(s <= 0) return(1e10)
  
  dens <- dGEV(losses, m, s, xi)
  
  if(any(!is.finite(dens)) || any(dens <= 0)) return(1e10)
  
  -sum(log(dens))
}

fit <- optim(iv, nll_safe, method = "L-BFGS-B",
             lower = c(-Inf, 1e-6, -Inf))

par <- fit$par

length(weekly_losses)
summary(weekly_losses)

#ASSESS THE FIT WITH DIAGNOSTIC PLOTS
#--- plot 1 ---
#Histogram vs fitted density

hist(weekly_losses, breaks = "FD",
     col = 'lightgoldenrod', 
     prob = TRUE,
     main = "Histogram of the largest weekly loss", 
     xlab = '', ylim = c(0, 60))

curve(dGEV(x, m = par[1], s = par[2], xi = par[3]), 
      from =min(weekly_losses), 
      to = max(weekly_losses), 
      lwd = 3, col = 'brown2', 
      add = TRUE)

#--- plot 2 ---
#COMPARING CDFs (Empirical vs fitted GEV)

# CDF comparison
plot(ecdf(weekly_losses), main = "Empirical vs GEV CDF",
     ylab = "CDF", lwd = 2)

curve(pGEV(x, par[1], par[2], par[3]),
      from =min(weekly_losses), 
      to = max(weekly_losses), 
      add = TRUE, col = "brown2", lwd = 2)



#--- plot 3 ---
#QQ plot
#careful qq plot uses package EnvStats with shape parameter as -psi

library(EnvStats)


EnvStats::qqPlot(weekly_losses,
                 distribution = "gevd",
                 param.list = list(location = par[1],
                                   scale = par[2],
                                   shape = -par[3]),
                 add.line = TRUE,
                 line.col = "brown3")

###############################
#Compute EVT VaR for fitted GEV
###############################


mu = 0.0235
sigma = 0.0112
xi = 0.15996
portfolio.value = 1000000
par = c(mu, sigma, xi)

#GEV quantile function (VaR engine)
qGEV <- function(p, mu, sigma, xi){
  if(sigma <= 0) stop("scale must be > 0")
  
  if(abs(xi) < 1e-8){
    return(mu - sigma * log(-log(p)))
  } else {
    return(mu + (sigma/xi) * ((-log(p))^(-xi) - 1))
  }
}

#Weekly GEV VaR (95% and 99%)
VaR_95_gev <- qGEV(0.95, par[1], par[2], par[3])
VaR_95_gev

VaR_95_gevvalue <- VaR_95*portfolio.value
VaR_95_gevvalue

VaR_99_gev <- qGEV(0.99, par[1], par[2], par[3])
VaR_99_gev

VaR_99_gevvalue <- VaR_99*portfolio.value
VaR_99_gevvalue

#Expected Shortfall
ES_empirical <- mean(weekly_losses[weekly_losses >= VaR_99])
ES_empirical

ES_value <- ES_empirical * portfolio.value
ES_value

###########################
#EFRM EVT MODELS PART II
###########################

###########################################
# POTS MODEL - EXCEEDANCES METHOD
###########################################

#Generate samples of exceedances
generate.ex <- function(n = 10000, u, rdist, param = c(0, 1), seed = 1){
  
  set.seed(seed)
  
  sample <- rdist(n, param[1], param[2]) #generate sample 
  
  W <- sample[sample > u] - u #retain values over 'u'and subtract 'u'
  
  hist(W, breaks = "FD", 
       col = 'lightgreen', 
       prob = TRUE,
       main = "Histogram of exceedances",
       xlab = ""
       )

  
}

#-----------------------------------------------------------
# --- For a UNIFORM, in MDA of GEV psi = -1 i.e weibull ---

#exceedances of the uniform are exactly uniform regardless of the threshold level

#u=1/2
generate.ex(u = 1/2, rdist = runif)


#u=0.9
generate.ex(u = 0.9, rdist = runif)

#---------------------------------------------------------

# --- For a NORMAL, in MDA of GEV psi = 0 i.e gumbel ---

#u=1
generate.ex(u = 1, rdist = rnorm)

#u=2
generate.ex(u = 2, rdist = rnorm)

#-----------------------------------------------------------

# --- For a STUDENT-T, in MDA of GEV psi > 0 i.e frechet ---

#u=1
generate.ex(u = 1, rdist = rt, param = c(3, 0))

#u=2
generate.ex(u = 2, rdist = rt, param = c(3, 0))


#-----------------------------------------------------------
#FITTING DATA TO GPD
#-----------------------------------------------------------

dGPD1 <- function(x, s, e){ #density for non-negative 'e'
 
   ifelse(x < 0, 0, 1/s*(1 + e*x/s)^(-1/e-1))
  
}


dGPD2 <- function(x, s, e) { #density for negative 'e'
  
  ifelse(x < 0, 0, ifelse (x > -s/e, 0, 1/s*(1 + e*x/s)^(-1/e-1)))
  
}
  

#---------------------------------------
#FUNCTION TO FIND ML ESTIMATES
#---------------------------------------

#Argument 'negative.e' specifies whether we want to fit with negative 'e' or not

fit.GPD <- function(my.data, u, negative.e = F){
  
  W <- my.data[my.data > u] - u #retain values for 'u', and subtract 'u'
  
  #negative log-likelihood of our W data (function of 3 parameters)
  if(negative.e == TRUE){
   
     nll <- function(p){-sum(log(dGPD2(W, p[1], p[2])))}
     
  } else{
    
    nll <- function(p){-sum(log(dGPD1(W, p[1], p[2])))}
    
  }
  
  #method of moments estimates
  e.mm <- 0.5*(mean(W^2) - 2*mean(W)^2)/(mean(W^2)+mean(W)^2)
  
  s.mm <- mean(W)*(1 - e.mm)
  
  #matrix of constraints: this is used in the optimization to force xi to be either >=0 or <=0
  #we also specify an initial 'e' estimate 
  if(negative.e == TRUE){
    
    mat <- rbind(c(1, 0), c(0, -1))
    
    e.i <- -1/1000 #set something negative but close to 0, otherwise initial optimization might fail 
  } else{
    
    mat <- rbind(c(1, 0), c(0, 1))
    
    e.i <- 1/1000 #set something negative but close to 0, otherwise initial optimization might fail
  }
  
  #minimize f (subject to constraints)
  fit <- constrOptim(c(s.mm, e.i), nll, ui = mat, ci = c(0, 0), grad = NULL)
  
  #display parameter estimates
  (par <- fit$par)
}


###########################################
#FITTING GPD TO SIMULATED DATA
###########################################

#Pareto random generator
rpareto <- function(n, alpha, lambda){
  
  Q <- function(u) {
    lambda*((1 - u)^(-1/alpha) -1)
  }
  
  return(Q(runif(n)))
}

set.seed(1)

n <- 5000

unif.data <- runif(n)

norm.data <- rnorm(n)

exp.data <- rexp(n)

gam.data <- rgamma(n, 1/2, 1/2)

pareto.data <- rpareto(n, 3, 4)

t.data <- rt(n, 3, 0)

#--- Fit to UNIFORM DATA ---
fit.GPD(unif.data, u = quantile(unif.data, 0.95), negative.e = TRUE)

#--- Fit to NORMAL DATA ---
fit.GPD(norm.data, u = quantile(norm.data, 0.95), negative.e = TRUE)

#--- Fit to EXPONENTIAL DATA ---
fit.GPD(exp.data, u = quantile(exp.data, 0.95))

#--- Fit to GAMMA DATA ---
fit.GPD(gam.data, u = quantile(gam.data, 0.95))

#--- Fit to PARETO DATA ---
fit.GPD(pareto.data, u = quantile(pareto.data, 0.95))

#--- Fit to STUDENT-T DATA ---
fit.GPD(t.data, u = quantile(t.data, 0.95))


#####################
#MEAN EXCESS PLOT
#####################

e.u <- function(sample, u){
  
  sum((sample - u)*(sample > u))/sum((sample > u))
}

#'vectorize' the function, so we can apply the function to a vector 'u'

e.u <- Vectorize(e.u, vectorize.args = 'u')


plot.e.u <- function(sample, nb.exclude = 0){
  
  #N is the no of points on the plot not sample size
  n <- length(sample) - (nb.exclude + 1)
  
  par(mar = c(5, 5, 2, 2))
  
  plot(sort(sample)[1:n], e.u(sample, sort(sample)[1:n]),
       pch = 19, cex = 1.5, col = 'brown3',
       xlab = "u",
       ylab = "e(u)",
       cex.lab = 1.75)
}

# --- UNIFORM MEAN EXCESS PLOT ---
plot.e.u(unif.data, 2)

# --- NORMAL MEAN EXCESS PLOT ---
plot.e.u(norm.data, 2)

# --- EXPONENTIAL MEAN EXCESS PLOT ---
plot.e.u(exp.data, 2)

# --- GAMMA MEAN EXCESS PLOT ---
plot.e.u(gam.data, 2)

# --- PARETO MEAN EXCESS PLOT ---
plot.e.u(pareto.data, 2)

# --- STUDENT-T MEAN EXCESS PLOT ---
plot.e.u(t.data, 2)


#-------------------------------------
#CHECKING THE FIT: REAL DATA EXAMPLE 
#-------------------------------------

#We use the "Danish fire" insurance data from package QRM

install.packages(c("QRM"), ("moments"))

library("QRM")

library("moments")

danishlosses <- danish$FIRE.LOSSES

summary(danishlosses)

sd(danishlosses)

kurtosis(danishlosses)

quantile(danishlosses, c(0.9, 0.95, 0.99))

hist(log(danishlosses), breaks = "FD",
     col = 'lightgoldenrod', prob = TRUE,
     main = "Histogram of Danish Fire Insurance: Log(losess)",
     xlab = "Log(losses)"
     )


#--- Mean excess plot of Danish fire data ---
par(mar = c(4,4,2,1))

plot.e.u(danishlosses)

#exclude the 2 largest data points
plot.e.u(danishlosses, nb.exclude = 2)

#plot appears fairly linear up to value 10 a kink appears and then it straightens out
#it seem reasonable to set u=10.
#let's see if a GPD fits the exceedances abouve u=10

threshold <- 10
(par <- fit.GPD(losses, u = threshold))

##########################################
#DIAGNOSTIC PLOTS TO ASSESS GPD MODEL FIT
##########################################

#--- Histogram vs fitted GPD plot ---
W <- losses[losses > threshold] - threshold

hist(W, col = 'lightgoldenrod',
     prob = TRUE, breaks = 50,
     main = "Histogram of the largest monthly loss",
     xlab = ""
     )

curve(dGPD1(x, s = par[1], e = par[2]),
      from = 0, to = max(losses),
      lwd = 3, col = 'brown3', add = TRUE
      )



#--- Comparing CDFs: Empirical vs GPD CDF plot ---
pGPD <- function(x, s, e){
  ifelse(x <= 0, 0, 1-(1+e*x/s)^(-1/e)
         ) #assumes positive 'e'
}

plot(ecdf(W), do.points = F,
     main = "Empirical vs GPD CDF",
     xlab = "",
     ylab = "CDF",
     lwd = 2
     )

curve(pGPD(x, s = par[1], e = par[2]), 
           from = 0, to = max(losses),
           col = 'brown2',lwd =2, add = TRUE
      )


#--- QQ plot ---

qGPD <- function(u, s, e){
  s*(-1 + (1 - u)^(-e))/e
}

m <- length(W) #sample size for exceedances

plot(qGPD((1:m)/(m+1), s = par[1], e = par[2]),
     sort(W), xlab =  "Fitted GPD quantiles",
     ylab = "Data quantiles"
     )

abline(0, 1, lwd = 2, col = 'brown2')


#overall the fit is good but form qqplot the 3 largest data point are larger than expected so we can do a KS test

ks.test(W, pGPD, s = par[1], e = par[2])

#############################################
# 27/04/26
# FITTING GPD TO MONTHLY PORTFOLIO LOSSES
#############################################

#---------------------------------------
# 1. DATA PREPARATION
#---------------------------------------

losses <- -min_returns$min_return

summary(losses)

hist(losses, breaks = "FD",
     col = "lightgoldenrod",
     prob = TRUE,
     main = "Extreme Monthly Portfolio Losses",
     xlab = "")

sd(losses)
kurtosis(losses)
quantile(losses, c(0.9, 0.95, 0.99))

hist(log(losses), breaks = "FD",
     col = "lightgoldenrod",
     prob = TRUE,
     main = "Log Extreme Losses",
     xlab = "Log(losses)")

#---------------------------------------
# 2. MEAN EXCESS PLOT
#---------------------------------------

par(mar = c(4,4,2,1))
plot.e.u(losses)

#---------------------------------------
# 3. THRESHOLD SELECTION
#---------------------------------------

u <- quantile(losses, 0.95)

# Fit initial GPD (optional check)
fit <- fit.GPD(losses, u = u)

#---------------------------------------
# 4. EXCEEDANCES
#---------------------------------------

W  <- losses[losses > u] - u
Nu <- length(W)
N  <- length(losses)

#---------------------------------------
# 5. MLE ESTIMATION (STABLE)
#---------------------------------------

gpd_nll <- function(p){
  sigma <- p[1]
  xi    <- p[2]
  
  if(sigma <= 0) return(1e10)
  
  z <- 1 + xi * W / sigma
  
  if(any(z <= 0)) return(1e10)
  
  -sum(log((1/sigma) * z^(-1/xi - 1)))
}

fit <- optim(
  par = c(sd(W), 0.1),
  fn = gpd_nll,
  method = "L-BFGS-B",
  lower = c(1e-6, -1),
  upper = c(Inf, 5)
)

sigma_hat <- fit$par[1]
xi_hat    <- fit$par[2]

#---------------------------------------
# 6. VaR FUNCTION (POT MODEL)
#---------------------------------------

VaR_GPD <- function(p, u, sigma, xi, N, Nu){
  u + (sigma/xi) * (((N/Nu * (1 - p))^(-xi)) - 1)
}

VaR_99  <- VaR_GPD(0.99,  u, sigma_hat, xi_hat, N, Nu)
VaR_995 <- VaR_GPD(0.995, u, sigma_hat, xi_hat, N, Nu)

VaR_99_value  <- VaR_99  * portfolio.value
VaR_995_value <- VaR_995 * portfolio.value

VaR_99
VaR_995

VaR_99_value
VaR_995_value

#---------------------------------------
# 7. ES FUNCTION (CORRECT)
#---------------------------------------

ES_GPD <- function(p, u, sigma, xi, N, Nu){
  VaR <- VaR_GPD(p, u, sigma, xi, N, Nu)
  
  ES <- (VaR + sigma - xi * u) / (1 - xi)
  
  return(ES)
}

ES_99 <- ES_GPD(0.99, u, sigma_hat, xi_hat, N, Nu)

ES_99_value <- ES_99 * portfolio.value

ES_99
ES_99_value

#---------------------------------------
# 8. BACKTESTING
#---------------------------------------

violations <- mean(losses > VaR_99)
violations


#---------------------------------------
# 9. DIAGNOSTIC PLOTS
#---------------------------------------

par(mar = c(4,4,2,1))

# --- Histogram fit ---
hist(W, prob = TRUE, breaks = 30,
     col = "lightgoldenrod",
     main = "GPD Fit to Exceedances",
     xlab = "")

curve((1/sigma_hat) * (1 + xi_hat*x/sigma_hat)^(-1/xi_hat - 1),
      from = 0, to = max(W),
      col = "brown3", lwd = 2, add = TRUE)

# --- CDF ---
pGPD <- function(x, sigma, xi){
  ifelse(x <= 0, 0,
         1 - (1 + xi*x/sigma)^(-1/xi))
}

plot(ecdf(W), do.points = FALSE,
     main = "Empirical vs GPD CDF",
     xlab = "",
     ylab = "CDF",
     lwd = 2)

curve(pGPD(x, sigma_hat, xi_hat),
      from = 0, to = max(W),
      col = "brown2", lwd = 2, add = TRUE)

# --- QQ plot ---
qGPD <- function(p, sigma, xi){
  sigma/xi * ((1 - p)^(-xi) - 1)
}

m <- length(W)

plot(qGPD((1:m)/(m+1), sigma_hat, xi_hat),
     sort(W),
     xlab = "Theoretical Quantiles",
     ylab = "Empirical Quantiles")

abline(0,1, col="red", lwd=2)

#---------------------------------------
# 10. KS TEST (DIAGNOSTIC ONLY)
#---------------------------------------

ks.test(W, pGPD, sigma_hat, xi_hat)
