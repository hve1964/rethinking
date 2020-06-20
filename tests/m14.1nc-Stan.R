################################################################################
# McElreath's (2020) multilevel linear model m14.1
# file: m14.1nc-Stan.R
# date: So, 07.06.2020
################################################################################

graphics.off()
rm(list=ls())

#-------------------------------------------------------------------------------
## Load packages
#-------------------------------------------------------------------------------
library(rstan)
rstan_options("required" = FALSE)
library(bayesplot)
library(loo)
#library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = "-march=corei7 -mtune=corei7")

## Online: for diagnosis and visualisations
library(shinystan)

#-------------------------------------------------------------------------------
## Load data
#-------------------------------------------------------------------------------
load("m14.1Data.RData")
d <- m14.1Data
rm(m14.1Data)

#-------------------------------------------------------------------------------
## Design matrix
#-------------------------------------------------------------------------------
X <- unname(
  model.matrix(
    object = wait ~ 1 + afternoon ,
    data = d
  )
) # "stats"
attr( X , "assign" ) <- NULL
head( X , n = 10 )

#-------------------------------------------------------------------------------
## Specify data list for Stan simulation
#-------------------------------------------------------------------------------
dataList <- list(
  N = nrow(X) ,
  M = ncol(X) ,
  K = length( unique(d$cafe) ) ,
  L = 1L ,
  u = matrix(
    data = rep(1.0, length(unique(d$cafe))),
    nrow = length(unique(d$cafe)),
    ncol = 1L
    ) ,
  afternoon = d$afternoon ,
  X = X ,
  wait = d$wait ,
  cafe = d$cafe #,
  #x_tilde = x_tilde ,
  #N_tilde = length(x_tilde)
)
#rm( d , x_tilde )
rm(X)

#-------------------------------------------------------------------------------
## Stan model: centered parametrisation of Gauß linear mode
#   McElreath's (2020) original model m14.1
#    (compiles alright: So, 07.06.2020)
#-------------------------------------------------------------------------------
set.seed(867530)

# Execution of Stan simulation: saving simulated posterior samples
mod0.stan <- stan(
  file = "m14.1.stan" ,
  data = dataList ,
  chains = 4 ,
  iter = 5000 ,
  warmup = 1000 ,
  thin = 1 ,
  init = "random" ,
  algorithm = "NUTS" ,
  control = list( adapt_delta = 0.95 ,
                  max_treedepth = 15 ) ,
  cores = 3
)

class(mod0.stan)
dim(mod0.stan)

check_hmc_diagnostics(mod0.stan)

print(
  x = mod0.stan ,
  pars = c("a", "a_cafe") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
stan_plot(
  object = mod0.stan ,
  pars = c("a", "a_cafe") ,
  ci_level = 0.89 ,
  outer_level = 0.97
)  # "rstan"
print(
  x = mod0.stan ,
  pars = c("b", "b_cafe") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
stan_plot(
  object = mod0.stan ,
  pars = c("b", "b_cafe") ,
  ci_level = 0.89 ,
  outer_level = 0.97
)  # "rstan"
print(
  x = mod0.stan ,
  pars = c("Rho") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
stan_plot(
  object = mod0.stan ,
  pars = c("Rho") ,
  ci_level = 0.89 ,
  outer_level = 0.97
)  # "rstan"

#-------------------------------------------------------------------------------
# Pareto-smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO)
#-------------------------------------------------------------------------------
loo(
  x = mod0.stan ,
  pars = "log_lik"
)  # "rstan"
looMod0 <- loo(
  x = mod0.stan ,
  pars = "log_lik"
)
plot(looMod2)

#-------------------------------------------------------------------------------
## Stan model: non-centered parametrisation of Gauß linear model
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Stan model code 1 (non-centered parametrisation, looped)
#    (compiles alright: Mi, 03.06.2020)
#     PSIS-LOO CV: ...
#-------------------------------------------------------------------------------
set.seed(867530)

# Execution of Stan simulation: saving simulated posterior samples
mod1.stan <- stan(
  file = "m14.1ncLoop.stan" ,
  data = dataList ,
  chains = 4 ,
  iter = 5000 ,
  warmup = 1000 ,
  thin = 1 ,
  init = "random" ,
  algorithm = "NUTS" ,
  control = list( adapt_delta = 0.95 ,
                  max_treedepth = 15 ) ,
  cores = 3
)

class(mod1.stan)
dim(mod1.stan)

check_hmc_diagnostics(mod1.stan)

stan_trace(
  object = mod1.stan ,
  pars = c("beta", "lp__") ,
  inc_warmup = TRUE
)  # "rstan"
print(
  x = mod1.stan ,
  pars = c("beta") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
stan_plot(
  object = mod1.stan ,
  pars = c("beta") ,
  ci_level = 0.89 ,
  outer_level = 0.97
)  # "rstan"
stan_dens(
  object = mod1.stan ,
  pars = c("beta", "lp__")
)  # "rstan"

stan_trace(
  object = mod1.stan ,
  pars = c("a", "b") ,
  inc_warmup = TRUE
)  # "rstan"
print(
  x = mod1.stan ,
  pars = c("a", "b") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
stan_plot(
  object = mod1.stan ,
  pars = c("a", "b") ,
  ci_level = 0.89 ,
  outer_level = 0.97
)  # "rstan"
stan_dens(
  object = mod1.stan ,
  pars = c("a", "b")
)  # "rstan"

print(
  x = mod1.stan ,
  pars = c("a", "a_cafe") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
print(
  x = mod1.stan ,
  pars = c("b", "b_cafe") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
stan_plot(
  object = mod1.stan ,
  pars = c("b", "b_cafe") ,
  ci_level = 0.89 ,
  outer_level = 0.97
)  # "rstan"
print(
  x = mod1.stan ,
  pars = c("L_Rho") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
stan_plot(
  object = mod1.stan ,
  pars = c("L_Rho") ,
  ci_level = 0.89 ,
  outer_level = 0.97
)  # "rstan"

stan_trace(
  object = mod1.stan ,
  pars = c("Rho") ,
  inc_warmup = TRUE
)  # "rstan"
print(
  x = mod1.stan ,
  pars = c("Rho") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
stan_plot(
  object = mod1.stan ,
  pars = c("Rho") ,
  ci_level = 0.89 ,
  outer_level = 0.97
)  # "rstan"
stan_dens(
  object = mod1.stan ,
  pars = c("Rho")
)  # "rstan"

#-------------------------------------------------------------------------------
# Stan model code 2 (non-centered parametrisation, fully vectorised)
#    (compiles alright: Do, 18.06.2020)
#     PSIS-LOO CV: ...
#-------------------------------------------------------------------------------
set.seed(867530)

# Execution of Stan simulation: saving simulated posterior samples
mod2.stan <- stan(
  file = "m14.1ncVec.stan" ,
  data = dataList ,
  chains = 4 ,
  iter = 5000 ,
  warmup = 1000 ,
  thin = 1 ,
  init = "random" ,
  algorithm = "NUTS" ,
  control = list( adapt_delta = 0.95 ,
                  max_treedepth = 15 ) ,
  cores = 3
)

class(mod2.stan)
dim(mod2.stan)

check_hmc_diagnostics(mod2.stan)
print(
  x = mod2.stan ,
  pars = c("beta") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
stan_trace(
  object = mod2.stan ,
  pars = c("beta") ,
  inc_warmup = TRUE
)  # "rstan"
stan_plot(
  object = mod2.stan ,
  pars = c("beta") ,
  ci_level = 0.89 ,
  outer_level = 0.97
)  # "rstan"

print(
  x = mod2.stan ,
  pars = c("L_Rho") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
stan_plot(
  object = mod2.stan ,
  pars = c("L_Rho") ,
  ci_level = 0.89 ,
  outer_level = 0.97
)  # "rstan"
print(
  x = mod2.stan ,
  pars = c("Rho") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
stan_plot(
  object = mod2.stan ,
  pars = c("Rho") ,
  ci_level = 0.89 ,
  outer_level = 0.97
)  # "rstan"

print(
  x = mod2.stan ,
  pars = c("a_cafe") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
stan_plot(
  object = mod2.stan ,
  pars = c("a_cafe") ,
  ci_level = 0.89 ,
  outer_level = 0.97
)  # "rstan"
print(
  x = mod2.stan ,
  pars = c("b_cafe") ,
  probs = c(0.015, 0.25, 0.50, 0.75, 0.985)
)
stan_plot(
  object = mod2.stan ,
  pars = c("b_cafe") ,
  ci_level = 0.89 ,
  outer_level = 0.97
)  # "rstan"

#-------------------------------------------------------------------------------
# Visualisation of fitted model
#-------------------------------------------------------------------------------
post <- extract(object = mod2.stan)
View(post)

a1 <- sapply( 1:length(unique(d$cafe)) ,
              function(i) { mean(d$wait[d$cafe==i & d$afternoon==0]) } )
b1 <- sapply( 1:length(unique(d$cafe)) ,
              function(i) { mean(d$wait[d$cafe==i & d$afternoon==1]) } ) - a1

a2 <- apply( post$a_cafe , 2 , mean )
b2 <- apply( post$b_cafe , 2 , mean )

plot( a1 ,b1 , xlab = "intercept" , ylab = "slope" ,
      pch = 16 , col = "blue" , ylim = c( min(b1)-0.1 , max(b1)+0.1 ) ,
      xlim=c( min(a1)-0.1 , max(a1)+0.1 ) )
points( a2 , b2 , pch = 1 )
for ( i in 1:length(unique(d$cafe)) ) {
  lines( c(a1[i], a2[i]) , c(b1[i], b2[i]) )
}

Mu_est <- c( mean(post$gamma[,1,1]) , mean(post$gamma[,1,2]) )
rho_est <- mean( post$Rho[,1,2] )
sa_est <- mean( post$sigma_cafe[,1] )
sb_est <- mean( post$sigma_cafe[,2] )
cov_ab <- sa_est*sb_est*rho_est
Sigma_est <- matrix( c(sa_est^2,cov_ab,cov_ab,sb_est^2) , ncol = 2 )

library(ellipse)
for ( l in c(0.1,0.3,0.5,0.8,0.99) ) {
  lines(ellipse(Sigma_est, centre = Mu_est, level = l) ,
        col = "green")
}

wait_morning_1 <- (a1)
wait_afternoon_1 <- (a1 + b1)
wait_morning_2 <- (a2)
wait_afternoon_2 <- (a2 + b2)

plot( wait_morning_1 , wait_afternoon_1 , xlab = "morning wait" ,
      ylab = "afternoon wait" , pch = 16 , col="blue" ,
      ylim = c( min(wait_afternoon_1)-0.1 , max(wait_afternoon_1)+0.1 ) ,
      xlim = c( min(wait_morning_1)-0.1 , max(wait_morning_1)+0.1 ) )
points( wait_morning_2 , wait_afternoon_2 , pch = 1 )
for ( i in 1:length(unique(d$cafe)) ) {
  lines( c(wait_morning_1[i], wait_morning_2[i]) ,
         c(wait_afternoon_1[i], wait_afternoon_2[i]) )
}
abline( a = 0 , b = 1 , lty = 2 )

library(MASS)
v <- mvrnorm( 1e4 , Mu_est , Sigma_est )
v[,2] <- v[,1] + v[,2] # calculate afternoon wait
Sigma_est2 <- cov(v)
Mu_est2 <- Mu_est
Mu_est2[2] <- Mu_est[1] + Mu_est[2]

for ( l in c(0.1,0.3,0.5,0.8,0.99) ) {
  lines(ellipse(Sigma_est2, centre = Mu_est2, level = l) ,
        col = "green")
}

#-------------------------------------------------------------------------------
# Pareto-smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO)
#-------------------------------------------------------------------------------
loo(
  x = mod2.stan ,
  pars = "log_lik"
)  # "rstan"
looMod2 <- loo(
  x = mod2.stan ,
  pars = "log_lik"
)
plot(looMod2)

#-------------------------------------------------------------------------------
## Session info
#-------------------------------------------------------------------------------
print( sessionInfo() )

################################################################################
################################################################################