# Chimp Nurturing Paper Analysis ####

# 1 Script intro ####
# aim of paper is to look at how maternal behaviour (nurturing/IBI) varies 
# in different mothers based on their:
# age, condition and current ecological circumstances; as well as offspring:
# age and sex

# This script includes:
# (a) Data exploration for IBI and behavioural models
# (b) IBI models and validation
# (i) Behavioural models and validation

# 2 Data, packages, functions ####

library(tidyverse)
library(reshape)
library(brms)
library(cmdstanr)
library(GGally)
library(beepr)
library(loo)
library(car)
library(paletteer)
library(rethinking)
library(patchwork)
library(bayestestR)
library(posterior)
library(tidybayes)
library(modelr)
library(MetBrewer)
library(bayesplot)
library(coda)

# function to extract mode from each row (via liz)
mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

set.seed(1977)

ls()

# mld = dataset for model 1 (IBI with rank data)
# m2d = dataset for model 2 (IBI with imputed rank data)
# m3d = dataset for model 3 (IBI with age data as proxy of rank)
# m4d = dataset for model 4 (Nurturing youngest offspring)
# m5d = dataset for model 5 (Nurturing any offspring)
# m6d = dataset for model 6 (IBI as predictor for nurturing)

# 3 IBI data exploration ####

str(m1d)

'$ community    : chr [1:68] Community of chimp
 $ name_code    : chr [1:68] Offspring name
 $ mother_name  : chr [1:68] Mother name
 $ offspring_sex: chr [1:68] Sex of offspring (male or female)
 $ parity       : num [1:68] Birth order of offspring
 $ off_sex      : chr [1:68] Duplicate of above for some reason
 $ IBI          : num [1:68] Inter-birth interval in years
 $ off_name     : chr [1:68] Another duplicate
 $ mother_age   : num [1:68] Age of mother in years
 $ rank         : num [1:68] Maternal rank (0 low; 1 high)
 $ rank2        : num [1:68] Imputed rank if rank unavailable (0 low; 1 high)
 $ zage         : num [1:68] -1.1495 -0.638 -0.5366 0.0783 0.285 ...
 $ zrank        : num [1:68] 0.402 0.55 -0.771 -0.328 0.965 ...'


names(m1d)
covees <- m1d[,c('mother_age', 'rank', 'off_sex')]
ggpairs(covees)
# weak correlation between rank 
# and mother's age: r = 0.215

corr_table <- as.data.frame(matrix(data = NA, nrow = 25, ncol = 3))

# 4 IBI model data prep ####

# one dataset with ranks, one without

# distribution of IBIs
hist(m1d$IBI)
hist(m2d$IBI)
hist(m3d$IBI)

# all should be ok gaussian

#need to scale rank and age
names(ibi1)

#age
m1d$zage = (m1d$mother_age - mean(m1d$mother_age))/sd(m1d$mother_age)
range(m1d$zage)
m2d$zage = (m2d$mother_age - mean(m2d$mother_age))/sd(m2d$mother_age)
range(m2d$zage)
m3d$zage = (m3d$mother_age - mean(m3d$mother_age))/sd(m3d$mother_age)
range(m3d$zage)

#rank
m1d$zrank = (m1d$rank - mean(m1d$rank))/sd(m1d$rank)
range(m1d$zrank)
m2d$zrank = (m2d$rank2 - mean(m2d$rank2))/sd(m2d$rank2)
range(m2d$zrank)
# not needed for model 3

length(unique(m1d$mother_name)) # 41 mothers for 68 data points
length(unique(m2d$mother_name)) # 47 mothers for 98 data points
length(unique(m3d$mother_name)) # 69 mothers for 135 data points

length(unique(m1d$name_code)) # 68 offspring for 68 data points
length(unique(m2d$name_code)) # 98 offspring for 98 data points
length(unique(m3d$name_code)) # 135 offspring for 135 data points

vif1 <- lm(IBI ~ zage + zrank + off_sex, data = m1d)
vif(vif1)
'    zage    zrank  off_sex 
1.083164 1.049149 1.033286 '

vif2 <- lm(IBI ~ zage + zrank + off_sex, data = m2d)
vif(vif2)
'zage    zrank  off_sex 
1.066856 1.029579 1.044141' 

vif3 <- lm(IBI ~ zage + off_sex, data = m3d)
vif(vif3)
'    zage  off_sex 
1.003919 1.003919 '

# 5 M1: IBIs with ranks ####

# first model is using only individuals we have information for
# let's try and do some prior predictive plots

set.seed(1977) 
N <- 100 # 100 lines
a <- rnorm( N , 0 , 1 )
b <- rnorm( N , 0 , 1 )

plot( NULL , xlim=range(m1d$mother_age) , ylim=c(0,20) , 
      xlab="age" , ylab="IBI" )
abline( h=0 , lty=2 )
abline( h=10 , lty=1 , lwd=0.5 )
xbar <- mean(m1d$mother_age)
for ( i in 1:N ) curve( a[i] + b[i]*(x - xbar) ,
                        from=min(m1d$mother_age) , to=max(m1d$mother_age) , add=TRUE)
# wow, super bad
# perhaps try the log normal as richard suggests in his book

a <- rlnorm( 1e4 , 1 , 0.5 )
dens( a )

b <- rnorm( 1e4 , 0 , 1 )
dens( b )

b <- rnorm( 1e4 , 0 , 0.075 )
dens( b )

set.seed(1977) 
N <- 100 # 100 lines
a <- rlnorm(N, 1 , 0.5 )
b <- rnorm(N, 0 , 1 )

plot( NULL , xlim=range(m1d$zage) , ylim=c(0,20) , 
      xlab="age" , ylab="IBI" )
abline( h=0 , lty=2 )
abline( h=10 , lty=1 , lwd=0.5 )
xbar <- mean(m1d$zage)
for ( i in 1:N ) curve( a[i] + b[i]*(x - xbar) ,
                        from=min(m1d$zage) , to=max(m1d$zage) , add=TRUE,
                        col=col.alpha("black",0.2))

# bit better but some negative IBIs... and some very low ones too

a <- rlnorm( 1e4 , 1.5 , 0.4 )
dens( a )

set.seed(1977) 
N <- 100 # 100 lines
a <- rlnorm(N, 1.5 , 0.4 )
b <- rnorm(N, 0 , 1 )

plot( NULL , xlim=range(m1d$zage) , ylim=c(0,20) , 
      xlab="age" , ylab="IBI" )
abline( h=0 , lty=2 )
abline( h=10 , lty=1 , lwd=0.5 )
xbar <- mean(m1d$zage)
for ( i in 1:N ) curve( a[i] + b[i]*(x - xbar) ,
                        from=min(m1d$zage) , to=max(m1d$zage) , add=TRUE,
                        col=col.alpha("black",0.2))

# looks a lot better; should also try for rank

set.seed(1977) 
N <- 100 # 100 lines
a <- rlnorm(N, 1.5 , 0.4 )
b <- rnorm(N, 0 , 1 )

plot( NULL , xlim=range(m1d$zrank) , ylim=c(0,20) , 
      xlab="rank" , ylab="IBI" )
abline( h=0 , lty=2 )
abline( h=10 , lty=1 , lwd=0.5 )
xbar <- mean(m1d$zrank)
for ( i in 1:N ) curve( a[i] + b[i]*(x - xbar) ,
                        from=min(m1d$zrank) , to=max(m1d$zrank) , add=TRUE,
                        col=col.alpha("black",0.2))

m1prior = c(prior(lognormal(1.5,0.4), class=Intercept),
            prior(normal(0,1), class=b),
            prior(exponential(1), class=sd))

make_stancode(
  IBI ~ zage + 
    zrank*off_sex +
    (1|mother_name),
  data = m1d,
  family = gaussian(link="identity"),
  prior = m1prior)


start_time <- Sys.time()
prior_check <- brm(
  IBI ~ zage + 
    zrank*off_sex +
    (1|mother_name),
  data = m1d, 
  family = gaussian(link="identity"), 
  prior = m1prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975), 
  sample_prior = "only")
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 21.59565 secs
mcmc_plot(prior_check)
rm(prior_check)

start_time <- Sys.time()
m1 <- brm(
  IBI ~ zage + 
    zrank*off_sex +
    (1|mother_name),
  data = m1d, 
  family = gaussian(link="identity"), 
  prior = m1prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 21.59565 secs

windows()
plot(m1)
summary(m1)

m1_ppc <- pp_check(m1, ndraws=100)
m1_ppc <- m1_ppc+labs(title="M1: Weakly regularising priors")
windows()
m1_ppc # looks ok

# 6 M1: uniform priors ####

m1_uniprior <- get_prior(IBI ~ zage + zrank*off_sex +
                           (1|mother_name),
                         data = m1d, family = gaussian(link="identity"))

make_stancode(
  IBI ~ zage + 
    zrank*off_sex +
    (1|mother_name),
  data = m1d,
  family = gaussian(link="identity"),
  prior = m1_uniprior)


start_time <- Sys.time()
m1uni <- brm(
  IBI ~ zage + 
    zrank*off_sex +
    (1|mother_name),
  data = m1d, 
  family = gaussian(link="identity"), 
  prior = m1_uniprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 26.59565 secs

plot(m1uni)
summary(m1uni)

m1uni_ppc <- pp_check(m1uni, ndraws=100)
m1uni_ppc <- m1uni_ppc+labs(title="M1: Uniform priors")
windows()
m1uni_ppc # looks ok

# 7 M1: strong priors ####

a <- rlnorm( 1e4 , 1 , 0.5 )
dens( a )

b <- rnorm( 1e4 , 0 , 1 )
dens( b )

b <- rnorm( 1e4 , 0 , 0.075 )
dens( b )

set.seed(1977) 
N <- 100 # 100 lines
a <- rlnorm(N, 1 , 0.25 )
b <- rnorm(N, 0 , 0.5 )

plot( NULL , xlim=range(m1d$zage) , ylim=c(0,20) , 
      xlab="age" , ylab="IBI" )
abline( h=0 , lty=2 )
abline( h=10 , lty=1 , lwd=0.5 )
xbar <- mean(m1d$zage)
for ( i in 1:N ) curve( a[i] + b[i]*(x - xbar) ,
                        from=min(m1d$zage) , to=max(m1d$zage) , add=TRUE,
                        col=col.alpha("black",0.2))

m1strprior = c(prior(lognormal(1,0.25), class=Intercept),
               prior(normal(0,0.5), class=b),
               prior(exponential(1), class=sd))

make_stancode(
  IBI ~ zage + 
    zrank*off_sex +
    (1|mother_name),
  data = m1d,
  family = gaussian(link="identity"),
  prior = m1strprior)


start_time <- Sys.time()
m1str <- brm(
  IBI ~ zage + 
    zrank*off_sex +
    (1|mother_name),
  data = m1d, 
  family = gaussian(link="identity"), 
  prior = m1strprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 38.25217 secs

plot(m1str)
summary(m1str)

m1str_ppc <- pp_check(m1str, ndraws=100)
m1str_ppc <- m1str_ppc+labs(title="M1: Strongly regularising priors")
windows()
m1str_ppc # looks ok

# 8 M2: Inter-birth interval model with imputed ranks ####

m2prior = c(prior(lognormal(1.5,0.4), class=Intercept),
            prior(normal(0,1), class=b),
            prior(exponential(1), class=sd))

make_stancode(
  IBI ~ zage + 
    zrank*off_sex +
    (1|mother_name),
  data = m2d,
  family = gaussian(link="identity"),
  prior = m2prior)

start_time <- Sys.time()
m2 <- brm(
  IBI ~ zage + 
    zrank*off_sex +
    (1|mother_name),
  data = m2d, 
  family = gaussian(link="identity"), 
  prior = m2prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 21.59565 secs

windows()
plot(m2)
summary(m2)

m2_ppc <- pp_check(m2, ndraws=100)
m2_ppc <- m2_ppc+labs(title="M2: Weakly regularising priors")
windows()

# 9 M1: uniform priors ####

m2_uniprior <- get_prior(IBI ~ zage + zrank*off_sex +
                           (1|mother_name),
                         data = m2d, family = gaussian(link="identity"))

make_stancode(
  IBI ~ zage + 
    zrank*off_sex +
    (1|mother_name),
  data = m2d,
  family = gaussian(link="identity"),
  prior = m2_uniprior)


start_time <- Sys.time()
m2uni <- brm(
  IBI ~ zage + 
    zrank*off_sex +
    (1|mother_name),
  data = m2d, 
  family = gaussian(link="identity"), 
  prior = m2_uniprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 26.59565 secs

plot(m2uni)
summary(m2uni)

m2uni_ppc <- pp_check(m2uni, ndraws=100)
m2uni_ppc <- m2uni_ppc+labs(title="M2: Uniform priors")
windows()

# 10 M2: strong priors ####

m2strprior = c(prior(lognormal(1,0.25), class=Intercept),
               prior(normal(0,0.5), class=b),
               prior(exponential(1), class=sd))

make_stancode(
  IBI ~ zage + 
    zrank*off_sex +
    (1|mother_name),
  data = m2d,
  family = gaussian(link="identity"),
  prior = m2strprior)

start_time <- Sys.time()
m2str <- brm(
  IBI ~ zage + 
    zrank*off_sex +
    (1|mother_name),
  data = m2d, 
  family = gaussian(link="identity"), 
  prior = m2strprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 38.25217 secs

plot(m2str)
summary(m2str)

m2str_ppc <- pp_check(m2str, ndraws=100)
m2str_ppc <- m2str_ppc+labs(title="M2: Strongly regularising priors")
windows()
m2str_ppc # looks ok

# 11 M3: Age rank proxy model ####

# let's try and do some prior predictive plots

set.seed(1977) 

# let's try the same prior as last time
a <- rlnorm( 1e4 , 1.5 , 0.4 )
dens( a )

set.seed(1977) 
N <- 100 # 100 lines
a <- rlnorm(N, 1.5 , 0.4 )
b <- rnorm(N, 0 , 1 )

plot( NULL , xlim=range(m3d$zage) , ylim=c(0,20) , 
      xlab="age" , ylab="IBI" )
abline( h=0 , lty=2 )
abline( h=10 , lty=1 , lwd=0.5 )
xbar <- mean(m3d$zage)
for ( i in 1:N ) curve( a[i] + b[i]*(x - xbar) ,
                        from=min(m3d$zage) , to=max(m3d$zage) , add=TRUE,
                        col=col.alpha("black",0.2))

# looks ok

m3prior = c(prior(lognormal(1.5,0.4), class=Intercept),
            prior(normal(0,1), class=b),
            prior(exponential(1), class=sd))

start_time <- Sys.time()
prior_check <- brm(
  IBI ~ zage*off_sex +
    (1|mother_name),
  data = m3d, 
  family = gaussian(link="identity"), 
  prior = m3prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975), 
  sample_prior = "only")
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 21.59565 secs
mcmc_plot(prior_check)
rm(prior_check)

start_time <- Sys.time()
m3 <- brm(
  IBI ~ zage*off_sex +
    (1|mother_name),
  data = m3d, 
  family = gaussian(link="identity"), 
  prior = m3prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 21.59565 secs

plot(m3)
summary(m3)

m3_ppc <- pp_check(m3, ndraws=100)
m3_ppc <- m3_ppc+labs(title="M3: Weakly regularising priors")
windows()
m3_ppc # looks ok

# 12 M3: uniform priors ####

m3_uniprior <- get_prior(IBI ~ zage*off_sex +
                           (1|mother_name),
                         data = m3d, family = gaussian(link="identity"))

make_stancode(
  IBI ~ zage*off_sex +
    (1|mother_name),
  data = m3d,
  family = gaussian(link="identity"),
  prior = m3_uniprior)


start_time <- Sys.time()
m3uni <- brm(
  IBI ~ zage*off_sex +
    (1|mother_name),
  data = m3d, 
  family = gaussian(link="identity"), 
  prior = m3_uniprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 26.59565 secs

plot(m3uni)
summary(m3uni)


m3uni_ppc <- pp_check(m3uni, ndraws=100)
m3uni_ppc <- m3uni_ppc+labs(title="M3: Uniform priors")
windows()
m3uni_ppc # looks ok

# 13 M3: strong priors ####

set.seed(1977) 
N <- 100 # 100 lines
a <- rlnorm(N, 1 , 0.25 )
b <- rnorm(N, 0 , 0.5 )

plot( NULL , xlim=range(m3d$zage) , ylim=c(0,20) , 
      xlab="age" , ylab="IBI" )
abline( h=0 , lty=2 )
abline( h=10 , lty=1 , lwd=0.5 )
xbar <- mean(m3d$zage)
for ( i in 1:N ) curve( a[i] + b[i]*(x - xbar) ,
                        from=min(m3d$zage) , to=max(m3d$zage) , add=TRUE,
                        col=col.alpha("black",0.2))

m3strprior = c(prior(lognormal(1,0.25), class=Intercept),
               prior(normal(0,0.5), class=b),
               prior(exponential(1), class=sd))

make_stancode(
  IBI ~ zage*off_sex +
    (1|mother_name),
  data = m3d,
  family = gaussian(link="identity"),
  prior = m3strprior)


start_time <- Sys.time()
m3str <- brm(
  IBI ~ zage*off_sex +
    (1|mother_name),
  data = m3d, 
  family = gaussian(link="identity"), 
  prior = m3strprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 38.25217 secs

plot(m3str)
summary(m3str)

m3str_ppc <- pp_check(m3str, ndraws=100)
m3str_ppc <- m3str_ppc+labs(title="M3: Strongly regularising priors")
windows()
m3str_ppc # looks ok

# 14 M4 data exploration ####

names(m4d)
str(m4d)

'$ subject             : chr  name of mother
$ rrdate              : chr  date (rdate format)
$ time.n              : num  start of behaviour: minutes since midnight
$ end                 : num  end of behaviour: minutes since midnight
$ observationID       : chr  ID of focal observation
$ day                 : num  30 30 30 30 30 30 30 30 30 30 ...
$ month               : num  12 12 12 12 12 12 12 12 12 12 ...
$ year                : num  2016 2016 2016 2016 2016 ...
$ sample_day          : num  ID of day of observation
$ youngest_offspring  : chr  ID of youngest offspring for subject
$ eldest_offspring_age: num  Age of eldest offspring of subject
$ scan                : chr  Social choice category
$ pregnant            : num  Refers to subject/mother; 1 = yes; 0 = no
$ youngest_off_age    : num  Age of youngest offspring of subject
$ youngest_off_sex    : chr  Sex of youngest offspring of subject
$ community           : chr  Community of subject
$ food_availability   : num  Index of food availability (0-3)
$ male.in.party       : int  Adult males in sub-party
$ female.in.party     : int  Adult females in sub-party
$ party_adults        : int  All adults sub-party
$ rank                : num  Rank of subject/mother (0 low; 1 high)
$ imm_kids            : num  Number of immature offspring of mother/subject
$ age                 : num  Age of subject/mother
$ database            : chr  Project where data came from
$ oestrous            : num  Mother in full sexual swelling (0 no; 1 yes)
$ frequency           : num  Variable just to get counts of different rows
$ zoffage             : num  0.652 0.652 0.652 0.652 0.652 ...
$ zage                : num  -0.969 -0.969 -0.969 -0.969 -0.969 ...
$ zrank               : num  0.267 0.267 0.267 0.267 0.267 ...
$ zfai                : num  0.00656 0.00656 0.00656 0.00656 0.00656 ...
$ ztime               : num  -1.25 -1.19 -1.14 -1.08 -1.02 ...
$ zparty              : num  -0.783 -1.267 -1.025 -0.783 -0.783 ...
$ zkids               : num  -0.725 -0.725 -0.725 -0.725 -0.725 ...
$ y                   : num  Numeric version of "scan"'

covees <- m4d[,c('youngest_off_age', 'age', 'rank',
                 'youngest_off_sex', 'food_availability',
                 'time.n', 'party_adults','pregnant',
                 'oestrous', 'imm_kids')]
windows()
ggpairs(covees)
# highest correlation between rank and mother's age: r = 0.685
# bit borderline so justifies model comparison with and without rank
rm(covees)

# 15 M4 data prep ####

# need to scale our continuous variables
names(m4d)

#ages
m4d$zoffage = (m4d$youngest_off_age - mean(m4d$youngest_off_age))/sd(m4d$youngest_off_age)
range(m4d$zoffage)
m4d$zage = (m4d$age - mean(m4d$age))/sd(m4d$age)
range(m4d$zage)

#rank
m4d$zrank = (m4d$rank - mean(m4d$rank))/sd(m4d$rank)
range(m4d$zrank)

#fai
m4d <- filter(m4d, !is.na(food_availability))
m4d$zfai = (m4d$food_availability - mean(m4d$food_availability))/sd(m4d$food_availability)
range(m4d$zfai)

#time
m4d <- filter(m4d, !is.na(time.n))
m4d$ztime = (m4d$time.n - mean(m4d$time.n))/sd(m4d$time.n)
range(m4d$ztime)

#party size
m4d <- filter(m4d, !is.na(party_adults))
m4d$zparty = (m4d$party_adults - mean(m4d$party_adults))/sd(m4d$party_adults)
range(m4d$zparty)

#other immature offspring
m4d <- filter(m4d, !is.na(imm_kids))
m4d$zkids = (m4d$imm_kids - mean(m4d$imm_kids))/sd(m4d$imm_kids)
range(m4d$zkids)

#check other variables
summary(m4d$pregnant)
m4d$pregnant[is.na(m4d$pregnant)]=0
summary(m4d$pregnant)
summary(m4d$oestrous)

# random effects
length(unique(m4d$youngest_offspring))# 42 levels
sort(unique(m4d$subject))# 31 levels
length(unique(m4d$observationID))
sort(unique(m4d$observationID))# 560 levels

# check vifs

m4d$y <- as.numeric(as.factor(m4d$scan))
sort(unique(m4d$scan))
range(m4d$y)

vif1 <- lm(y ~ 
             zoffage + 
             zage + 
             zrank + 
             youngest_off_sex +
             zfai + 
             ztime + 
             zparty + 
             oestrous +
             pregnant+
             zkids, data = m4d)
vif2 <- as.data.frame(vif(vif1))
#2.55 all good

# 16 M4: Nurturing youngest offspring model ####

# sorting priors for this might be tricky
# let's see what brms comes up with and start without correlations among all the slopes

m4prior <- get_prior(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
    (1|observationID),
  data = m4d,
  family = categorical(link="logit"))

m4prior$prior[7:20] <- "normal(0,1)"
m4prior$prior[44:55] <- "normal(0,1)"
m4prior$prior[77:90] <- "normal(0,1)"
m4prior$prior[112:125] <- "normal(0,1)"

m4prior2 = c(prior(normal(0,1), class=Intercept, dpar=mubfoodshareoffspring),
             prior(normal(0,1), class=Intercept, dpar=mucnurseoffspring),
             prior(normal(0,1), class=Intercept, dpar=mudsociooffspring),
             prior(normal(0,1), class=Intercept, dpar=muesocioothers),
             prior(normal(0,1), class=b, dpar=mubfoodshareoffspring),
             prior(normal(0,1), class=b, dpar=mucnurseoffspring),
             prior(normal(0,1), class=b, dpar=mudsociooffspring),
             prior(normal(0,1), class=b, dpar=muesocioothers),
             prior(exponential(1), class=sd, dpar=mubfoodshareoffspring),
             prior(exponential(1), class=sd, dpar=mucnurseoffspring),
             prior(exponential(1), class=sd, dpar=mudsociooffspring),
             prior(exponential(1), class=sd, dpar=muesocioothers))


make_stancode(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
    (1|observationID),
  data = m4d,
  family = categorical(link="logit"),
  prior = m4prior)

make_stancode(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
    (1|observationID),
  data = m4d,
  family = categorical(link="logit"),
  prior = m4prior2)

start_time <- Sys.time()
prior_check1 <- brm(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
    (1|observationID),
  data = m4d,
  family = categorical(link="logit"),
  prior = m4prior,
  chains = 1, 
  #threads = threading(3),
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 2500, 
  warmup = 500, 
  thin = 2,
  control = list(adapt_delta = 0.975), 
  sample_prior = "only"
)

prior_check2 <- brm(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
    (1|observationID),
  data = m4d,
  family = categorical(link="logit"),
  prior = m4prior2,
  chains = 1, 
  #threads = threading(3),
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 2500, 
  warmup = 500, 
  thin = 2,
  control = list(adapt_delta = 0.975), 
  sample_prior = "only"
)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time

mcmc_plot(prior_check1, variable = "^b_", regex=TRUE)
mcmc_plot(prior_check2, variable = "^b_", regex=TRUE)
# so prior one seems to be more regularising

m4prior3 = c(prior(normal(0,5), class=Intercept, dpar=mubfoodshareoffspring),
             prior(normal(0,5), class=Intercept, dpar=mucnurseoffspring),
             prior(normal(0,5), class=Intercept, dpar=mudsociooffspring),
             prior(normal(0,5), class=Intercept, dpar=muesocioothers),
             prior(normal(0,1), class=b, dpar=mubfoodshareoffspring),
             prior(normal(0,1), class=b, dpar=mucnurseoffspring),
             prior(normal(0,1), class=b, dpar=mudsociooffspring),
             prior(normal(0,1), class=b, dpar=muesocioothers),
             prior(exponential(1), class=sd, dpar=mubfoodshareoffspring),
             prior(exponential(1), class=sd, dpar=mucnurseoffspring),
             prior(exponential(1), class=sd, dpar=mudsociooffspring),
             prior(exponential(1), class=sd, dpar=muesocioothers))

prior_check3 <- brm(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
    (1|observationID),
  data = m4d,
  family = categorical(link="logit"),
  prior = m4prior3,
  chains = 1, 
  #threads = threading(3),
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 2500, 
  warmup = 500, 
  thin = 2,
  control = list(adapt_delta = 0.975), 
  sample_prior = "only"
)

mcmc_plot(prior_check3, variable = "^b_", regex=TRUE)
# ok, so this at least make sense as being less regularising

m4prior4 = c(prior(normal(0,1), class=b, dpar=mubfoodshareoffspring),
             prior(normal(0,1), class=b, dpar=mucnurseoffspring),
             prior(normal(0,1), class=b, dpar=mudsociooffspring),
             prior(normal(0,1), class=b, dpar=muesocioothers),
             prior(exponential(1), class=sd, dpar=mubfoodshareoffspring),
             prior(exponential(1), class=sd, dpar=mucnurseoffspring),
             prior(exponential(1), class=sd, dpar=mudsociooffspring),
             prior(exponential(1), class=sd, dpar=muesocioothers))

prior_check4 <- brm(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
    (1|observationID),
  data = m4d,
  family = categorical(link="logit"),
  prior = m4prior4,
  chains = 1, 
  #threads = threading(3),
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 2500, 
  warmup = 500, 
  thin = 2,
  control = list(adapt_delta = 0.975), 
  sample_prior = "only"
)

mcmc_plot(prior_check4, variable = "^b_", regex=TRUE)
mcmc_plot(prior_check3, variable = "^b_", regex=TRUE)
mcmc_plot(prior_check2, variable = "^b_", regex=TRUE)

# ok, so I don't know how to make a generative model of something this complicated
# so will just have to work with these different priors and see how we get on

m4prior = c(prior(normal(0,1), class=Intercept, dpar=mubfoodshareoffspring),
            prior(normal(0,1), class=Intercept, dpar=mucnurseoffspring),
            prior(normal(0,1), class=Intercept, dpar=mudsociooffspring),
            prior(normal(0,1), class=Intercept, dpar=muesocioothers),
            prior(normal(0,1), class=b, dpar=mubfoodshareoffspring),
            prior(normal(0,1), class=b, dpar=mucnurseoffspring),
            prior(normal(0,1), class=b, dpar=mudsociooffspring),
            prior(normal(0,1), class=b, dpar=muesocioothers),
            prior(exponential(1), class=sd, dpar=mubfoodshareoffspring),
            prior(exponential(1), class=sd, dpar=mucnurseoffspring),
            prior(exponential(1), class=sd, dpar=mudsociooffspring),
            prior(exponential(1), class=sd, dpar=muesocioothers))

make_stancode(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
    (1|observationID),
  data = m4d,
  family = categorical(link="logit"),
  prior = m4prior)


start_time <- Sys.time()
m4 <- brm(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
    (1|observationID),
  data = m4d,
  family = categorical(link="logit"),
  prior = m4prior,
  chains = 3, 
  cores = 9,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
# Time difference of 2.061436 days

plot(m4)
summary(m4)

# now the ppc plot
# we need a plot with four columns per scan type
# two columns will be a 0 and 1 for the observed data for that scan
# other column will be the mode of the posterior draws for a given observation, with 0/1 scan type

pred_vals = predict(m4, summary = FALSE)
predictions = as.data.frame(t(pred_vals))
pred_vals2 <- as.data.frame(pred_vals)

smode = mode(predictions)
range(smode)
colnames(smode) <- c("predictions")
sort(unique(m4d$scan))

ppc1d <- m4d
ppc1d$predictions <- smode$predictions
head(ppc1d)

xx <- aggregate(ppc1d$subject, by=list(ppc1d$scan),length)
yy <- aggregate(ppc1d$subject, by=list(ppc1d$predictions),length)

xx$mode_preds <- yy$x
colnames(xx) <- c("scan", "observed", "predicted")

ppc_bard <- xx
rm(xx,yy)

pred_vals2$non_social = rowSums(pred_vals2 == "1")
pred_vals2$foodshare = rowSums(pred_vals2 == "2")
pred_vals2$nurse = rowSums(pred_vals2 == "3")
pred_vals2$socioff = rowSums(pred_vals2 == "4")
pred_vals2$socioth = rowSums(pred_vals2 == "5")

quantile(pred_vals2$non_social,prob=c(0.0,1.0))
# 8502 8985
quantile(pred_vals2$foodshare,prob=c(0.0,1.0))
# 220 339
quantile(pred_vals2$nurse,prob=c(0.0,1.0))
# 708 972
quantile(pred_vals2$socioff,prob=c(0.0,1.0))
# 946 1268
quantile(pred_vals2$socioth,prob=c(0.0,1.0))
# 654 970

ppc_bard$lower <- c(8502,220,708,946,654)
ppc_bard$upper <- c(8985,339,972,1268,970)

# reformat to make easier for plotting

xx <- ppc_bard$observed
yy <- ppc_bard$predicted
zz <- c(xx,yy)
xx <- rep("observed",5)
yy <- rep("posterior",5)
zzz <- c(xx,yy)
xxx <- c("Non social", "Foodshare with offspring",
         "Nurse", "Sociopositive with offspring",
         "Sociopositive with others")
xxx <- c(xxx,xxx)

test <- as.data.frame(zzz)
test$scan <- xxx
test$count <- zz
lower <- ppc_bard$lower
lower <- c(NA,NA,NA,NA,NA,lower)
test$lower <- lower
upper <- ppc_bard$upper
upper <- c(NA,NA,NA,NA,NA,upper)
test$upper <- upper

ppc_bard <- test
rm(test,lower,upper,xx,xxx,yy,zz,zzz)
ppc_bard$y <- ppc_bard$zzz

ppc_bard$y <- factor(ppc_bard$y, levels=c("posterior","observed"))
ppc_bard$scan <- factor(ppc_bard$scan, 
                        levels=c("Non social","Foodshare with offspring",
                                 "Nurse","Sociopositive with offspring",
                                 "Sociopositive with others"))

m4_ppc1 <- ggplot(ppc_bard, 
                  aes(x=scan, y=count, fill=y)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 13))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))

m4_ppc1 <- m4_ppc1 + scale_fill_brewer(palette="Paired") + 
  theme_classic()

m4_ppc1 <- m4_ppc1+labs(title="M4:  Weakly regularising priors",x="Social choice",
                        y="Count")
m4_ppc1 <- m4_ppc1+xlab("Social choice")
m4_ppc1 <- m4_ppc1+theme(plot.title = element_text(face="bold", hjust = 0.5))
windows()
m4_ppc1

# 17 M4: Youngest offspring model - uniform prior ####

m4uniprior = get_prior(scan ~
                         (zoffage+I(zoffage^2))+
                         zage+
                         zrank*youngest_off_sex+
                         zfai+
                         community+
                         ztime+
                         zparty+
                         zkids+
                         pregnant+
                         oestrous+
                         (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
                         (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
                         (1|observationID),
                       data = m4d,
                       family = categorical(link="logit"))

m4uniprior$prior[c(21,56,91,126)] <- ""
m4uniprior$prior[c(22,57,92,127)] <- "exponential(1)"

make_stancode(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
    (1|observationID),
  data = m4d,
  family = categorical(link="logit"),
  prior = m4uniprior)


start_time <- Sys.time()
m4_uni <- brm(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
    (1|observationID),
  data = m4d,
  family = categorical(link="logit"),
  prior = m4uniprior,
  chains = 3, 
  cores = 9,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
# Time difference of 2.146646 days

plot(m4_uni)
summary(m4_uni)

# now the ppc plot
# we need a plot with four columns per scan type
# two columns will be a 0 and 1 for the observed data for that scan
# other column will be the mode of the posterior draws for a given observation, with 0/1 scan type

pred_vals = predict(m4_uni, summary = FALSE)
predictions = as.data.frame(t(pred_vals))
pred_vals2 <- as.data.frame(pred_vals)

smode = mode(predictions)
range(smode)
colnames(smode) <- c("predictions")
sort(unique(m4d$scan))

ppc1d <- m4d
ppc1d$predictions <- smode$predictions
head(ppc1d)
beep(8)
xx <- aggregate(ppc1d$subject, by=list(ppc1d$scan),length)
yy <- aggregate(ppc1d$subject, by=list(ppc1d$predictions),length)

xx$mode_preds <- yy$x
colnames(xx) <- c("scan", "observed", "predicted")

ppc_bard <- xx
rm(xx,yy)

pred_vals2$non_social = rowSums(pred_vals2 == "1")
pred_vals2$foodshare = rowSums(pred_vals2 == "2")
pred_vals2$nurse = rowSums(pred_vals2 == "3")
pred_vals2$socioff = rowSums(pred_vals2 == "4")
pred_vals2$socioth = rowSums(pred_vals2 == "5")

quantile(pred_vals2$non_social,prob=c(0.0,1.0))
# 8535 9002
quantile(pred_vals2$foodshare,prob=c(0.0,1.0))
# 210 391
quantile(pred_vals2$nurse,prob=c(0.0,1.0))
# 704 977
quantile(pred_vals2$socioff,prob=c(0.0,1.0))
# 939 1264
quantile(pred_vals2$socioth,prob=c(0.0,1.0))
# 665 958

ppc_bard$lower <- c(8535,210,704,939,665)
ppc_bard$upper <- c(9002,391,977,1264,958)

# reformat to make easier for plotting

xx <- ppc_bard$observed
yy <- ppc_bard$predicted
zz <- c(xx,yy)
xx <- rep("observed",5)
yy <- rep("posterior",5)
zzz <- c(xx,yy)
xxx <- c("Non social", "Foodshare with offspring",
         "Nurse", "Sociopositive with offspring",
         "Sociopositive with others")
xxx <- c(xxx,xxx)

test <- as.data.frame(zzz)
test$scan <- xxx
test$count <- zz
lower <- ppc_bard$lower
lower <- c(NA,NA,NA,NA,NA,lower)
test$lower <- lower
upper <- ppc_bard$upper
upper <- c(NA,NA,NA,NA,NA,upper)
test$upper <- upper

ppc_bard <- test
rm(test,lower,upper,xx,xxx,yy,zz,zzz)
ppc_bard$y <- ppc_bard$zzz

ppc_bard$y <- factor(ppc_bard$y, levels=c("posterior","observed"))
ppc_bard$scan <- factor(ppc_bard$scan, 
                        levels=c("Non social","Foodshare with offspring",
                                 "Nurse","Sociopositive with offspring",
                                 "Sociopositive with others"))

m4uni_ppc1 <- ggplot(ppc_bard, 
                     aes(x=scan, y=count, fill=y)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 13))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))

m4uni_ppc1 <- m4uni_ppc1 + scale_fill_brewer(palette="Paired") + 
  theme_classic()

m4uni_ppc1 <- m4uni_ppc1+labs(title="M4:  Uniform priors",x="Social choice",
                              y="Count")
m4uni_ppc1 <- m4uni_ppc1+theme(plot.title = element_text(face="bold", hjust = 0.5))
windows()
m4uni_ppc1

# 18 M4: Youngest offspring model - strongly regularising priors ####

m4strprior = c(prior(normal(0,0.5), class=Intercept, dpar=mubfoodshareoffspring),
               prior(normal(0,0.5), class=Intercept, dpar=mucnurseoffspring),
               prior(normal(0,0.5), class=Intercept, dpar=mudsociooffspring),
               prior(normal(0,0.5), class=Intercept, dpar=muesocioothers),
               prior(normal(0,0.5), class=b, dpar=mubfoodshareoffspring),
               prior(normal(0,0.5), class=b, dpar=mucnurseoffspring),
               prior(normal(0,0.5), class=b, dpar=mudsociooffspring),
               prior(normal(0,0.5), class=b, dpar=muesocioothers),
               prior(exponential(2), class=sd, dpar=mubfoodshareoffspring),
               prior(exponential(2), class=sd, dpar=mucnurseoffspring),
               prior(exponential(2), class=sd, dpar=mudsociooffspring),
               prior(exponential(2), class=sd, dpar=muesocioothers))

make_stancode(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
    (1|observationID),
  data = m4d,
  family = categorical(link="logit"),
  prior = m4strprior)

start_time <- Sys.time()
m4_str <- brm(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|youngest_offspring)+
    (1|observationID),
  data = m4d,
  family = categorical(link="logit"),
  prior = m4strprior,
  chains = 3, 
  cores = 9,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
# Time difference of 1.921808 days

plot(m4_str)
summary(m4_str)

# now the ppc plot
# we need a plot with four columns per scan type
# two columns will be a 0 and 1 for the observed data for that scan
# other column will be the mode of the posterior draws for a given observation, with 0/1 scan type

pred_vals = predict(m4_str, summary = FALSE)
predictions = as.data.frame(t(pred_vals))
pred_vals2 <- as.data.frame(pred_vals)

smode = mode(predictions)
range(smode)
colnames(smode) <- c("predictions")
sort(unique(m4d$scan))

ppc1d <- m4d
ppc1d$predictions <- smode$predictions
head(ppc1d)
beep(8)

xx <- aggregate(ppc1d$subject, by=list(ppc1d$scan),length)
yy <- aggregate(ppc1d$subject, by=list(ppc1d$predictions),length)

xx$mode_preds <- yy$x
colnames(xx) <- c("scan", "observed", "predicted")

ppc_bard <- xx
rm(xx,yy)

pred_vals2$non_social = rowSums(pred_vals2 == "1")
pred_vals2$foodshare = rowSums(pred_vals2 == "2")
pred_vals2$nurse = rowSums(pred_vals2 == "3")
pred_vals2$socioff = rowSums(pred_vals2 == "4")
pred_vals2$socioth = rowSums(pred_vals2 == "5")

quantile(pred_vals2$non_social,prob=c(0.0,1.0))
# 8471 8940
quantile(pred_vals2$foodshare,prob=c(0.0,1.0))
# 235 416
quantile(pred_vals2$nurse,prob=c(0.0,1.0))
# 700 995
quantile(pred_vals2$socioff,prob=c(0.0,1.0))
# 963 1281
quantile(pred_vals2$socioth,prob=c(0.0,1.0))
# 681 941

ppc_bard$lower <- c(8471,235,700,963,681)
ppc_bard$upper <- c(8940,416,995,1281,941)

# reformat to make easier for plotting

xx <- ppc_bard$observed
yy <- ppc_bard$predicted
zz <- c(xx,yy)
xx <- rep("observed",5)
yy <- rep("posterior",5)
zzz <- c(xx,yy)
xxx <- c("Non social", "Foodshare with offspring",
         "Nurse", "Sociopositive with offspring",
         "Sociopositive with others")
xxx <- c(xxx,xxx)

test <- as.data.frame(zzz)
test$scan <- xxx
test$count <- zz
lower <- ppc_bard$lower
lower <- c(NA,NA,NA,NA,NA,lower)
test$lower <- lower
upper <- ppc_bard$upper
upper <- c(NA,NA,NA,NA,NA,upper)
test$upper <- upper

ppc_bard <- test
rm(test,lower,upper,xx,xxx,yy,zz,zzz)
ppc_bard$y <- ppc_bard$zzz

ppc_bard$y <- factor(ppc_bard$y, levels=c("posterior","observed"))
ppc_bard$scan <- factor(ppc_bard$scan, 
                        levels=c("Non social","Foodshare with offspring",
                                 "Nurse","Sociopositive with offspring",
                                 "Sociopositive with others"))

m4str_ppc1 <- ggplot(ppc_bard, 
                     aes(x=scan, y=count, fill=y)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 13))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))

m4str_ppc1 <- m4str_ppc1 + scale_fill_brewer(palette="Paired") + 
  theme_classic()

m4str_ppc1 <- m4str_ppc1+labs(title="M4:  Strongly regularising priors",x="Social choice",
                              y="Count")
m4str_ppc1 <- m4str_ppc1+theme(plot.title = element_text(face="bold", hjust = 0.5))

# 19 M5 data exploration ####

names(m5d)
str(m5d)
# same as m4d except for "new_scan" variable
# this is the social choice of the mother also taking into 
# account her choices in relation to her older immature offspring

names(m5d)

covees <- m5d[,c('eldest_offspring_age','youngest_off_age', 'age', 'rank',
                 'youngest_off_sex','eldest_off_sex', 'food_availability',
                 'time.n', 'party_adults',
                 'imm_kids')]
windows()
ggpairs(covees)
# highest correlation between rank and mother's age: r = 0.653
# bit borderline so justifies model comparison with and without rank
rm(covees)

# 20 M5 data prep ####

# need to scale our continuous variables
names(m5d)

#ages
m5d$zeoffage = (m5d$eldest_offspring_age - mean(m5d$eldest_offspring_age))/sd(m5d$eldest_offspring_age)
range(m5d$zeoffage)
m5d$zyoffage = (m5d$youngest_off_age - mean(m5d$youngest_off_age))/sd(m5d$youngest_off_age)
range(m5d$zyoffage)
m5d$zage = (m5d$age - mean(m5d$age))/sd(m5d$age)
range(m5d$zage)

#rank
m5d$zrank = (m5d$rank - mean(m5d$rank))/sd(m5d$rank)
range(m5d$zrank)

#fai
m5d <- filter(m5d, !is.na(food_availability))
m5d$zfai = (m5d$food_availability - mean(m5d$food_availability))/sd(m5d$food_availability)
range(m5d$zfai)

#time
m5d <- filter(m5d, !is.na(time.n))
m5d$ztime = (m5d$time.n - mean(m5d$time.n))/sd(m5d$time.n)
range(m5d$ztime)

#party size
m5d <- filter(m5d, !is.na(party_adults))
m5d$zparty = (m5d$party_adults - mean(m5d$party_adults))/sd(m5d$party_adults)
range(m5d$zparty)

#other immature offspring
m5d <- filter(m5d, !is.na(imm_kids))
m5d$zkids = (m5d$imm_kids - mean(m5d$imm_kids))/sd(m5d$imm_kids)
range(m5d$zkids)

#check other variables
summary(m5d$pregnant)
hist(m5d$pregnant)
summary(m5d$oestrous)
hist(m5d$oestrous)
length(m5d$oestrous[m5d$oestrous==1])

# random effects
names(m5d)
sort(unique(m5d$youngest_offspring))#19 levels
sort(unique(m5d$eldest_offspring)) # 20 levels - so no different combinations of offspring
sort(unique(m5d$subject))# 17 levels
length(unique(m5d$subject))
length(unique(m5d$observationID))
sort(unique(m5d$observationID))# 224 levels

table(m5d$subject, m5d$youngest_offspring)
table(m5d$subject, m5d$offspring_list)

# only Rewenzori appears with more than one youngest offspring, so no need for that
# random effect in this model
# Rwenzori and Toumai only mothers with different combinations of offspring 
# so we only need mother ID and observation ID as random effects.

# check vifs

m5d$y <- as.numeric(as.factor(m5d$new_scan))
range(m5d$y)

vif1 <- lm(y ~ 
             zeoffage + 
             zyoffage + 
             zage + 
             zrank + 
             eldest_off_sex +
             youngest_off_sex +
             zfai + 
             ztime + 
             zparty + 
             zkids, data = m5d)
vif2 <- as.data.frame(vif(vif1))

# 21 M5: Dependent offspring model ####

# going to run a model with default priors to see if it runs
# also with simpler random effect structure

m5prior <- get_prior(
  new_scan ~
    (zyoffage+I(zyoffage^2))+
    zeoffage+
    zage+
    zrank*eldest_off_sex+
    zrank+
    eldest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    (1|subject)+
    (1|observationID),
  data = m5d,
  family = categorical(link="logit")
)

make_stancode(
  new_scan ~
    (zyoffage+I(zyoffage^2))+
    zeoffage+
    zage+
    zrank*youngest_off_sex+
    zrank+
    eldest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    (1|subject)+
    (1|observationID),
  data = m5d,
  family = categorical(link="logit"),
  prior = m5prior)

start_time <- Sys.time()
m5 <- brm(
  new_scan ~
    (zyoffage+I(zyoffage^2))+
    zeoffage+
    zage+
    zrank*youngest_off_sex+
    zrank+
    eldest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    (1|subject)+
    (1|observationID),
  data = m5d,
  family = categorical(link="logit"),
  prior = m5prior,
  chains = 3, 
  cores = 9,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
# Time difference of 4.044248 hours

# ok, so that runs with good diagnostics... let's run it with our weakly regularising priors
sort(unique(m5d$new_scan))
m5prior = c(prior(normal(0,1), class=Intercept, dpar=mubnutrinonyoungest),
            prior(normal(0,1), class=Intercept, dpar=mucsociononyoungest),
            prior(normal(0,1), class=Intercept, dpar=mudnurtureboth),
            prior(normal(0,1), class=Intercept, dpar=muenurtureyoungest),
            prior(normal(0,1), class=Intercept, dpar=mufsocialnonoffspring),
            prior(normal(0,1), class=b, dpar=mubnutrinonyoungest),
            prior(normal(0,1), class=b, dpar=mucsociononyoungest),
            prior(normal(0,1), class=b, dpar=mudnurtureboth),
            prior(normal(0,1), class=b, dpar=muenurtureyoungest),
            prior(normal(0,1), class=b, dpar=mufsocialnonoffspring),
            prior(exponential(1), class=sd, dpar=mubnutrinonyoungest),
            prior(exponential(1), class=sd, dpar=mucsociononyoungest),
            prior(exponential(1), class=sd, dpar=mudnurtureboth),
            prior(exponential(1), class=sd, dpar=muenurtureyoungest),
            prior(exponential(1), class=sd, dpar=mufsocialnonoffspring)
)

make_stancode(
  new_scan ~
    (zyoffage+I(zyoffage^2))+
    zeoffage+
    zage+
    zrank*youngest_off_sex+
    zrank*eldest_off_sex+
    zrank+
    eldest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    ((zyoffage+I(zyoffage^2))+
       zeoffage+
       zage+ztime|subject)+
    (1|observationID),
  data = m5d,
  family = categorical(link="logit"),
  prior = m5prior)

start_time <- Sys.time()
m5 <- brm(
  new_scan ~
    (zyoffage+I(zyoffage^2))+
    zeoffage+
    zage+
    zrank*youngest_off_sex+
    zrank*eldest_off_sex+
    zrank+
    eldest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    ((zyoffage+I(zyoffage^2))+
       zeoffage+
       zage+ztime|subject)+
    (1|observationID),
  data = m5d,
  family = categorical(link="logit"),
  prior = m5prior,
  chains = 3, 
  cores = 9,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
# Time difference of 10.83071 hours

plot(m5)
summary(m5)

# now the ppc plot
# we need a plot with four columns per scan type
# two columns will be a 0 and 1 for the observed data for that scan
# other column will be the mode of the posterior draws for a given observation, with 0/1 scan type

pred_vals = predict(m5, summary = FALSE)
predictions = as.data.frame(t(pred_vals))
pred_vals2 <- as.data.frame(pred_vals)

smode = mode(predictions)
range(smode)
colnames(smode) <- c("predictions")

ppc5d <- m5d
ppc5d$predictions <- smode$predictions
head(ppc5d)

xx <- aggregate(ppc5d$subject, by=list(ppc5d$new_scan),length)
yy <- aggregate(ppc5d$subject, by=list(ppc5d$predictions),length)

xx$mode_preds <- yy$x
colnames(xx) <- c("scan", "observed", "predicted")

ppc_bard <- xx
rm(xx,yy)

pred_vals2$non_social = rowSums(pred_vals2 == "1")
pred_vals2$nutri_eldest = rowSums(pred_vals2 == "2")
pred_vals2$socio_eldest = rowSums(pred_vals2 == "3")
pred_vals2$nurture_both = rowSums(pred_vals2 == "4")
pred_vals2$nurture_youngest = rowSums(pred_vals2 == "5")
pred_vals2$social_others = rowSums(pred_vals2 == "6")

quantile(pred_vals2$non_social,prob=c(0.0,1.0))
# 2883 3166
quantile(pred_vals2$nutri_eldest,prob=c(0.0,1.0))
# 15 77
quantile(pred_vals2$socio_eldest,prob=c(0.0,1.0))
# 95 215
quantile(pred_vals2$nurture_both,prob=c(0.0,1.0))
# 51 153
quantile(pred_vals2$nurture_youngest,prob=c(0.0,1.0))
# 572 802
quantile(pred_vals2$social_others,prob=c(0.0,1.0))
# 169 315

ppc_bard$lower <- c(2883,15,95,51,572,169)
ppc_bard$upper <- c(3166,77,215,153,802,315)

# reformat to make easier for plotting

xx <- ppc_bard$observed
yy <- ppc_bard$predicted
zz <- c(xx,yy)
xx <- rep("observed",6)
yy <- rep("posterior",6)
zzz <- c(xx,yy)
xxx <- c("Non social", "Nutritional with eldest", "Sociopositive with eldest",
         "Nurture both", "Nurture youngest",
         "Sociopositive with others")
xxx <- c(xxx,xxx)

test <- as.data.frame(zzz)
test$scan <- xxx
test$count <- zz
lower <- ppc_bard$lower
lower <- c(NA,NA,NA,NA,NA,NA,lower)
test$lower <- lower
upper <- ppc_bard$upper
upper <- c(NA,NA,NA,NA,NA,NA,upper)
test$upper <- upper

ppc_bard <- test
rm(test,lower,upper,xx,xxx,yy,zz,zzz)
ppc_bard$y <- ppc_bard$zzz

ppc_bard$y <- factor(ppc_bard$y, levels=c("posterior","observed"))
ppc_bard$scan <- factor(ppc_bard$scan, 
                        levels=c("Non social","Nutritional with eldest","Sociopositive with eldest",
                                 "Nurture both","Nurture youngest",
                                 "Sociopositive with others"))

m5ppc1 <- ggplot(ppc_bard, 
                 aes(x=scan, y=count, fill=y)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 13))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))

m5ppc1 <- m5ppc1 + scale_fill_brewer(palette="Paired") + 
  theme_classic()

m5ppc1 <- m5ppc1+labs(title="M5:  Weakly regularising priors",x="Social choice",
                      y="Count")
m5ppc1 <- m5ppc1+xlab("Social choice")
m5ppc1 <- m5ppc1+theme(plot.title = element_text(face="bold", hjust = 0.5))

# 22 M5: Uniform prior ####

m5uniprior = get_prior(new_scan ~
                         (zyoffage+I(zyoffage^2))+
                         zeoffage+
                         zage+
                         zrank*youngest_off_sex+
                         zrank*eldest_off_sex+
                         zrank+
                         eldest_off_sex+
                         zfai+
                         community+
                         ztime+
                         zparty+
                         zkids+
                         ((zyoffage+I(zyoffage^2))+
                            zeoffage+
                            zage+ztime|subject)+
                         (1|observationID),
                       data = m5d,
                       family = categorical(link="logit"))


m5uniprior = c(prior(exponential(1), class=sd, dpar=mubnutrinonyoungest),
               prior(exponential(1), class=sd, dpar=mucsociononyoungest),
               prior(exponential(1), class=sd, dpar=mudnurtureboth),
               prior(exponential(1), class=sd, dpar=muenurtureyoungest),
               prior(exponential(1), class=sd, dpar=mufsocialnonoffspring)
)


#m5uniprior$prior[c(21,48,75,102,129)] <- ""
#m5uniprior$prior[c(22,49,76,103,130)] <- "exponential(1)"

make_stancode(
  new_scan ~
    (zyoffage+I(zyoffage^2))+
    zeoffage+
    zage+
    zrank*youngest_off_sex+
    zrank*eldest_off_sex+
    zrank+
    eldest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    ((zyoffage+I(zyoffage^2))+
       zeoffage+
       zage+ztime|subject)+
    (1|observationID),
  data = m5d,
  family = categorical(link="logit"),
  prior = m5uniprior)

start_time <- Sys.time()
m5uni <- brm(
  new_scan ~
    (zyoffage+I(zyoffage^2))+
    zeoffage+
    zage+
    zrank*youngest_off_sex+
    zrank*eldest_off_sex+
    zrank+
    eldest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    ((zyoffage+I(zyoffage^2))+
       zeoffage+
       zage+ztime|subject)+
    (1|observationID),
  data = m5d,
  family = categorical(link="logit"),
  prior = m5uniprior,
  chains = 3, 
  cores = 9,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
mod_time
# Time difference of 1.150235 days
beep(8)

plot(m5uni)
summary(m5uni)

# now the ppc plot
# we need a plot with four columns per scan type
# two columns will be a 0 and 1 for the observed data for that scan
# other column will be the mode of the posterior draws for a given observation, with 0/1 scan type

pred_vals = predict(m5uni, summary = FALSE)
predictions = as.data.frame(t(pred_vals))
pred_vals2 <- as.data.frame(pred_vals)

smode = mode(predictions)
range(smode)
colnames(smode) <- c("predictions")

ppc5d <- m5d
ppc5d$predictions <- smode$predictions
head(ppc5d)

xx <- aggregate(ppc5d$subject, by=list(ppc5d$new_scan),length)
yy <- aggregate(ppc5d$subject, by=list(ppc5d$predictions),length)

xx$mode_preds <- yy$x
colnames(xx) <- c("scan", "observed", "predicted")

ppc_bard <- xx
rm(xx,yy)

pred_vals2$non_social = rowSums(pred_vals2 == "1")
pred_vals2$nutri_eldest = rowSums(pred_vals2 == "2")
pred_vals2$socio_eldest = rowSums(pred_vals2 == "3")
pred_vals2$nurture_both = rowSums(pred_vals2 == "4")
pred_vals2$nurture_youngest = rowSums(pred_vals2 == "5")
pred_vals2$social_others = rowSums(pred_vals2 == "6")

quantile(pred_vals2$non_social,prob=c(0.0,1.0))
# 2887 3212
quantile(pred_vals2$nutri_eldest,prob=c(0.0,1.0))
# 13 66
quantile(pred_vals2$socio_eldest,prob=c(0.0,1.0))
# 97 221
quantile(pred_vals2$nurture_both,prob=c(0.0,1.0))
# 40 148
quantile(pred_vals2$nurture_youngest,prob=c(0.0,1.0))
# 560 800
quantile(pred_vals2$social_others,prob=c(0.0,1.0))
# 167 321

ppc_bard$lower <- c(2887,13,97,40,560,167)
ppc_bard$upper <- c(3212,66,221,148,800,321)

# reformat to make easier for plotting

xx <- ppc_bard$observed
yy <- ppc_bard$predicted
zz <- c(xx,yy)
xx <- rep("observed",6)
yy <- rep("posterior",6)
zzz <- c(xx,yy)
xxx <- c("Non social", "Nutritional with eldest", "Sociopositive with eldest",
         "Nurture both", "Nurture youngest",
         "Sociopositive with others")
xxx <- c(xxx,xxx)

test <- as.data.frame(zzz)
test$scan <- xxx
test$count <- zz
lower <- ppc_bard$lower
lower <- c(NA,NA,NA,NA,NA,NA,lower)
test$lower <- lower
upper <- ppc_bard$upper
upper <- c(NA,NA,NA,NA,NA,NA,upper)
test$upper <- upper

ppc_bard <- test
rm(test,lower,upper,xx,xxx,yy,zz,zzz)
ppc_bard$y <- ppc_bard$zzz

ppc_bard$y <- factor(ppc_bard$y, levels=c("posterior","observed"))
ppc_bard$scan <- factor(ppc_bard$scan, 
                        levels=c("Non social","Nutritional with eldest","Sociopositive with eldest",
                                 "Nurture both","Nurture youngest",
                                 "Sociopositive with others"))

m5unippc1 <- ggplot(ppc_bard, 
                    aes(x=scan, y=count, fill=y)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 13))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))

m5unippc1 <- m5unippc1 + scale_fill_brewer(palette="Paired") + 
  theme_classic()

m5unippc1 <- m5unippc1+labs(title="M5:  Uniform priors",x="Social choice",
                            y="Count")
m5unippc1 <- m5unippc1+xlab("Social choice")
m5unippc1 <- m5unippc1+theme(plot.title = element_text(face="bold", hjust = 0.5))

# 23 M5: Strong prior ####

m5strprior = c(prior(normal(0,0.5), class=Intercept, dpar=mubnutrinonyoungest),
               prior(normal(0,0.5), class=Intercept, dpar=mucsociononyoungest),
               prior(normal(0,0.5), class=Intercept, dpar=mudnurtureboth),
               prior(normal(0,0.5), class=Intercept, dpar=muenurtureyoungest),
               prior(normal(0,0.5), class=Intercept, dpar=mufsocialnonoffspring),
               prior(normal(0,0.5), class=b, dpar=mubnutrinonyoungest),
               prior(normal(0,0.5), class=b, dpar=mucsociononyoungest),
               prior(normal(0,0.5), class=b, dpar=mudnurtureboth),
               prior(normal(0,0.5), class=b, dpar=muenurtureyoungest),
               prior(normal(0,0.5), class=b, dpar=mufsocialnonoffspring),
               prior(exponential(2), class=sd, dpar=mubnutrinonyoungest),
               prior(exponential(2), class=sd, dpar=mucsociononyoungest),
               prior(exponential(2), class=sd, dpar=mudnurtureboth),
               prior(exponential(2), class=sd, dpar=muenurtureyoungest),
               prior(exponential(2), class=sd, dpar=mufsocialnonoffspring)
)

make_stancode(
  new_scan ~
    (zyoffage+I(zyoffage^2))+
    zeoffage+
    zage+
    zrank*youngest_off_sex+
    zrank*eldest_off_sex+
    zrank+
    eldest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    ((zyoffage+I(zyoffage^2))+
       zeoffage+
       zage+ztime|subject)+
    (1|observationID),
  data = m5d,
  family = categorical(link="logit"),
  prior = m5strprior)

start_time <- Sys.time()
m5str <- brm(
  new_scan ~
    (zyoffage+I(zyoffage^2))+
    zeoffage+
    zage+
    zrank*youngest_off_sex+
    zrank*eldest_off_sex+
    zrank+
    eldest_off_sex+
    zfai+
    community+
    ztime+
    zparty+
    zkids+
    ((zyoffage+I(zyoffage^2))+
       zeoffage+
       zage+ztime|subject)+
    (1|observationID),
  data = m5d,
  family = categorical(link="logit"),
  prior = m5strprior,
  chains = 3, 
  cores = 9,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
# Time difference of 6.939731 hours
saveRDS(m5str,"E:/OneDrive/Project_data/chimp_nurturing/m5str.RDS")

plot(m5str)
summary(m5str)

# now the ppc plot
# we need a plot with four columns per scan type
# two columns will be a 0 and 1 for the observed data for that scan
# other column will be the mode of the posterior draws for a given observation, with 0/1 scan type

pred_vals = predict(m5str, summary = FALSE)
predictions = as.data.frame(t(pred_vals))
pred_vals2 <- as.data.frame(pred_vals)

smode = mode(predictions)
range(smode)
colnames(smode) <- c("predictions")

ppc5d <- m5d
ppc5d$predictions <- smode$predictions
head(ppc5d)

xx <- aggregate(ppc5d$subject, by=list(ppc5d$new_scan),length)
yy <- aggregate(ppc5d$subject, by=list(ppc5d$predictions),length)

xx$mode_preds <- yy$x
colnames(xx) <- c("scan", "observed", "predicted")

ppc_bard <- xx
rm(xx,yy)

pred_vals2$non_social = rowSums(pred_vals2 == "1")
pred_vals2$nutri_eldest = rowSums(pred_vals2 == "2")
pred_vals2$socio_eldest = rowSums(pred_vals2 == "3")
pred_vals2$nurture_both = rowSums(pred_vals2 == "4")
pred_vals2$nurture_youngest = rowSums(pred_vals2 == "5")
pred_vals2$social_others = rowSums(pred_vals2 == "6")

quantile(pred_vals2$non_social,prob=c(0.0,1.0))
# 2854 3136
quantile(pred_vals2$nutri_eldest,prob=c(0.0,1.0))
# 16 75
quantile(pred_vals2$socio_eldest,prob=c(0.0,1.0))
# 108 233
quantile(pred_vals2$nurture_both,prob=c(0.0,1.0))
# 45 150
quantile(pred_vals2$nurture_youngest,prob=c(0.0,1.0))
# 572 799
quantile(pred_vals2$social_others,prob=c(0.0,1.0))
# 182 332

ppc_bard$lower <- c(2854,16,108,45,572,182)
ppc_bard$upper <- c(3136,75,233,150,799,332)

# reformat to make easier for plotting

xx <- ppc_bard$observed
yy <- ppc_bard$predicted
zz <- c(xx,yy)
xx <- rep("observed",6)
yy <- rep("posterior",6)
zzz <- c(xx,yy)
xxx <- c("Non social", "Nutritional with eldest", "Sociopositive with eldest",
         "Nurture both", "Nurture youngest",
         "Sociopositive with others")
xxx <- c(xxx,xxx)

test <- as.data.frame(zzz)
test$scan <- xxx
test$count <- zz
lower <- ppc_bard$lower
lower <- c(NA,NA,NA,NA,NA,NA,lower)
test$lower <- lower
upper <- ppc_bard$upper
upper <- c(NA,NA,NA,NA,NA,NA,upper)
test$upper <- upper

ppc_bard <- test
rm(test,lower,upper,xx,xxx,yy,zz,zzz)
ppc_bard$y <- ppc_bard$zzz

ppc_bard$y <- factor(ppc_bard$y, levels=c("posterior","observed"))
ppc_bard$scan <- factor(ppc_bard$scan, 
                        levels=c("Non social","Nutritional with eldest","Sociopositive with eldest",
                                 "Nurture both","Nurture youngest",
                                 "Sociopositive with others"))

m5strppc1 <- ggplot(ppc_bard, 
                    aes(x=scan, y=count, fill=y)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 13))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))

m5strppc1 <- m5unippc1 + scale_fill_brewer(palette="Paired") + 
  theme_classic()

m5strppc1 <- m5strppc1+labs(title="M5:  Strongly regularising priors",x="Social choice",
                            y="Count")
m5strppc1 <- m5strppc1+xlab("Social choice")
m5strppc1 <- m5strppc1+theme(plot.title = element_text(face="bold", hjust = 0.5))

# 24 M6 data exploration ####

# and now we add in their IBIs from ibi
names(m6d)
str(m6d)

# same as m4d but now includes "ibi"; 
# this is the inter-birth interval for the youngest offspring

m6d <- filter(m6d, !is.na(ibi))

sort(unique(m6d$subject)) # 15 mothers
length(unique(m6d$subject))


names(m6d)
covees <- m6d[,c('ibi','youngest_off_age', 'age', 'rank',
                 'youngest_off_sex','food_availability',
                 'time.n', 'party_adults','pregnant',
                 'oestrous', 'imm_kids')]
windows()
ggpairs(covees)
# highest correlation between rank and mother's age: r = 0.685
# bit borderline so justifies model comparison with and without rank
rm(covees)

m6d$zibi = (m6d$ibi - mean(m6d$ibi))/sd(m6d$ibi)
range(m6d$zibi)

#check other variables
summary(m6d$pregnant)
hist(m6d$pregnant)
summary(m6d$oestrous)
hist(m6d$oestrous)

# random effects
names(m6d)
length(unique(m6d$youngest_offspring))# 16 levels
sort(unique(m6d$youngest_offspring))
sort(unique(m6d$subject))# 15 levels
length(unique(m6d$subject))
length(unique(m6d$observationID))
sort(unique(m6d$observationID))# 228 levels
length(unique(m6d$subject[m6d$pregnant==1])) # 9 different pregnant females
length(unique(m6d$subject[m6d$oestrous==1])) # 10 different oestrous females

table(m6d$subject, m6d$youngest_offspring)
table(m5d$subject, m5d$offspring_list)

# only Asanti appears with more than one youngest offspring here 

# check vifs

m6d$y <- as.numeric(as.factor(m6d$scan))
range(m6d$y)

vif1 <- lm(y ~ 
             zoffage + 
             zage + 
             zrank + 
             zibi +
             youngest_off_sex +
             zfai + 
             ztime + 
             zparty +
             pregnant + 
             oestrous + 
             zkids, data = m6d)
vif2 <- as.data.frame(vif(vif1))

# 25 M6: IBI nurturing predictor model ####

m6prior = c(prior(normal(0,1), class=Intercept, dpar=mubfoodshareoffspring),
            prior(normal(0,1), class=Intercept, dpar=mucnurseoffspring),
            prior(normal(0,1), class=Intercept, dpar=mudsociooffspring),
            prior(normal(0,1), class=Intercept, dpar=muesocioothers),
            prior(normal(0,1), class=b, dpar=mubfoodshareoffspring),
            prior(normal(0,1), class=b, dpar=mucnurseoffspring),
            prior(normal(0,1), class=b, dpar=mudsociooffspring),
            prior(normal(0,1), class=b, dpar=muesocioothers),
            prior(exponential(1), class=sd, dpar=mubfoodshareoffspring),
            prior(exponential(1), class=sd, dpar=mucnurseoffspring),
            prior(exponential(1), class=sd, dpar=mudsociooffspring),
            prior(exponential(1), class=sd, dpar=muesocioothers))

make_stancode(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    zibi+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (1|observationID),
  data = m6d,
  family = categorical(link="logit"),
  prior = m6prior)


start_time <- Sys.time()
m6 <- brm(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    zibi+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (1|observationID),
  data = m6d,
  family = categorical(link="logit"),
  prior = m6prior,
  chains = 3, 
  cores = 9,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
# Time difference of 5.855842 hours

summary(m6)
plot(m6)

# now the ppc plot
# we need a plot with four columns per scan type
# two columns will be a 0 and 1 for the observed data for that scan
# other column will be the mode of the posterior draws for a given observation, with 0/1 scan type

pred_vals = predict(m6, summary = FALSE)
predictions = as.data.frame(t(pred_vals))
pred_vals2 <- as.data.frame(pred_vals)

smode = mode(predictions)
range(smode)
colnames(smode) <- c("predictions")
sort(unique(m6d$scan))

ppc1d <- m6d
ppc1d$predictions <- smode$predictions
head(ppc1d)

xx <- aggregate(ppc1d$subject, by=list(ppc1d$scan),length)
yy <- aggregate(ppc1d$subject, by=list(ppc1d$predictions),length)

xx$mode_preds <- yy$x
colnames(xx) <- c("scan", "observed", "predicted")

ppc_bard <- xx
rm(xx,yy)

pred_vals2$non_social = rowSums(pred_vals2 == "1")
pred_vals2$foodshare = rowSums(pred_vals2 == "2")
pred_vals2$nurse = rowSums(pred_vals2 == "3")
pred_vals2$socioff = rowSums(pred_vals2 == "4")
pred_vals2$socioth = rowSums(pred_vals2 == "5")

quantile(pred_vals2$non_social,prob=c(0.0,1.0))
# 2749 3076
quantile(pred_vals2$foodshare,prob=c(0.0,1.0))
# 116 236
quantile(pred_vals2$nurse,prob=c(0.0,1.0))
# 341 546
quantile(pred_vals2$socioff,prob=c(0.0,1.0))
# 452 671
quantile(pred_vals2$socioth,prob=c(0.0,1.0))
# 220 402

ppc_bard$lower <- c(2749,116,341,452,220)
ppc_bard$upper <- c(3076,236,546,671,402)

# reformat to make easier for plotting

xx <- ppc_bard$observed
yy <- ppc_bard$predicted
zz <- c(xx,yy)
xx <- rep("observed",5)
yy <- rep("posterior",5)
zzz <- c(xx,yy)
xxx <- c("Non social", "Foodshare with offspring",
         "Nurse", "Sociopositive with offspring",
         "Sociopositive with others")
xxx <- c(xxx,xxx)

test <- as.data.frame(zzz)
test$scan <- xxx
test$count <- zz
lower <- ppc_bard$lower
lower <- c(NA,NA,NA,NA,NA,lower)
test$lower <- lower
upper <- ppc_bard$upper
upper <- c(NA,NA,NA,NA,NA,upper)
test$upper <- upper

ppc_bard <- test
rm(test,lower,upper,xx,xxx,yy,zz,zzz)
ppc_bard$y <- ppc_bard$zzz

ppc_bard$y <- factor(ppc_bard$y, levels=c("posterior","observed"))
ppc_bard$scan <- factor(ppc_bard$scan, 
                        levels=c("Non social","Foodshare with offspring",
                                 "Nurse","Sociopositive with offspring",
                                 "Sociopositive with others"))

m6_ppc1 <- ggplot(ppc_bard, 
                  aes(x=scan, y=count, fill=y)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 13))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))

m6_ppc1 <- m6_ppc1 + scale_fill_brewer(palette="Paired") + 
  theme_classic()

m6_ppc1 <- m6_ppc1+labs(title="M6:  Weakly regularising priors",x="Social choice",
                        y="Count")
m6_ppc1 <- m6_ppc1+xlab("Social choice")
m6_ppc1 <- m6_ppc1+theme(plot.title = element_text(face="bold", hjust = 0.5))

# 26 M6 Uniform prior ####

m6uniprior = get_prior(scan ~
                         (zoffage+I(zoffage^2))+
                         zage+
                         zrank*youngest_off_sex+
                         zfai+
                         zibi+
                         community+
                         ztime+
                         zparty+
                         zkids+
                         pregnant+
                         oestrous+
                         (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
                         (1|observationID),
                       data = m6d,
                       family = categorical(link="logit"))

m4uniprior$prior[c(21,49,77,105)] <- ""
m4uniprior$prior[c(22,50,78,106)] <- "exponential(1)"

make_stancode(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    zibi+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (1|observationID),
  data = m6d,
  family = categorical(link="logit"),
  prior = m6uniprior)


start_time <- Sys.time()
m6_uni <- brm(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    zibi+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (1|observationID),
  data = m6d,
  family = categorical(link="logit"),
  prior = m6uniprior,
  chains = 3, 
  cores = 9,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
# Time difference of 10.36316 hours

plot(m6_uni)
summary(m6_uni)

# now the ppc plot
# we need a plot with four columns per scan type
# two columns will be a 0 and 1 for the observed data for that scan
# other column will be the mode of the posterior draws for a given observation, with 0/1 scan type

pred_vals = predict(m6_uni, summary = FALSE)
predictions = as.data.frame(t(pred_vals))
pred_vals2 <- as.data.frame(pred_vals)

smode = mode(predictions)
range(smode)
colnames(smode) <- c("predictions")
sort(unique(m6d$scan))

ppc1d <- m6d
ppc1d$predictions <- smode$predictions
head(ppc1d)

xx <- aggregate(ppc1d$subject, by=list(ppc1d$scan),length)
yy <- aggregate(ppc1d$subject, by=list(ppc1d$predictions),length)

xx$mode_preds <- yy$x
colnames(xx) <- c("scan", "observed", "predicted")

ppc_bard <- xx
rm(xx,yy)

pred_vals2$non_social = rowSums(pred_vals2 == "1")
pred_vals2$foodshare = rowSums(pred_vals2 == "2")
pred_vals2$nurse = rowSums(pred_vals2 == "3")
pred_vals2$socioff = rowSums(pred_vals2 == "4")
pred_vals2$socioth = rowSums(pred_vals2 == "5")

quantile(pred_vals2$non_social,prob=c(0.0,1.0))
# 2796 3109
quantile(pred_vals2$foodshare,prob=c(0.0,1.0))
# 111 232
quantile(pred_vals2$nurse,prob=c(0.0,1.0))
# 348 543
quantile(pred_vals2$socioff,prob=c(0.0,1.0))
# 456 657
quantile(pred_vals2$socioth,prob=c(0.0,1.0))
# 219 423

ppc_bard$lower <- c(2796,111,348,456,219)
ppc_bard$upper <- c(3109,232,543,657,423)

# reformat to make easier for plotting

xx <- ppc_bard$observed
yy <- ppc_bard$predicted
zz <- c(xx,yy)
xx <- rep("observed",5)
yy <- rep("posterior",5)
zzz <- c(xx,yy)
xxx <- c("Non social", "Foodshare with offspring",
         "Nurse", "Sociopositive with offspring",
         "Sociopositive with others")
xxx <- c(xxx,xxx)

test <- as.data.frame(zzz)
test$scan <- xxx
test$count <- zz
lower <- ppc_bard$lower
lower <- c(NA,NA,NA,NA,NA,lower)
test$lower <- lower
upper <- ppc_bard$upper
upper <- c(NA,NA,NA,NA,NA,upper)
test$upper <- upper

ppc_bard <- test
rm(test,lower,upper,xx,xxx,yy,zz,zzz)
ppc_bard$y <- ppc_bard$zzz

ppc_bard$y <- factor(ppc_bard$y, levels=c("posterior","observed"))
ppc_bard$scan <- factor(ppc_bard$scan, 
                        levels=c("Non social","Foodshare with offspring",
                                 "Nurse","Sociopositive with offspring",
                                 "Sociopositive with others"))

m6uni_ppc1 <- ggplot(ppc_bard, 
                     aes(x=scan, y=count, fill=y)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 13))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))

m6uni_ppc1 <- m6uni_ppc1 + scale_fill_brewer(palette="Paired") + 
  theme_classic()

m6uni_ppc1 <- m6uni_ppc1+labs(title="M6:  Uniform priors",x="Social choice",
                              y="Count")
m6uni_ppc1 <- m6uni_ppc1+theme(plot.title = element_text(face="bold", hjust = 0.5))

# 27 M6: strongly regularising priors ####

m6strprior = c(prior(normal(0,0.5), class=Intercept, dpar=mubfoodshareoffspring),
               prior(normal(0,0.5), class=Intercept, dpar=mucnurseoffspring),
               prior(normal(0,0.5), class=Intercept, dpar=mudsociooffspring),
               prior(normal(0,0.5), class=Intercept, dpar=muesocioothers),
               prior(normal(0,0.5), class=b, dpar=mubfoodshareoffspring),
               prior(normal(0,0.5), class=b, dpar=mucnurseoffspring),
               prior(normal(0,0.5), class=b, dpar=mudsociooffspring),
               prior(normal(0,0.5), class=b, dpar=muesocioothers),
               prior(exponential(2), class=sd, dpar=mubfoodshareoffspring),
               prior(exponential(2), class=sd, dpar=mucnurseoffspring),
               prior(exponential(2), class=sd, dpar=mudsociooffspring),
               prior(exponential(2), class=sd, dpar=muesocioothers))

make_stancode(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    zibi+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (1|observationID),
  data = m6d,
  family = categorical(link="logit"),
  prior = m6strprior)


start_time <- Sys.time()
m6_str <- brm(
  scan ~
    (zoffage+I(zoffage^2))+
    zage+
    zrank*youngest_off_sex+
    zfai+
    zibi+
    community+
    ztime+
    zparty+
    zkids+
    pregnant+
    oestrous+
    (zoffage + I(zoffage^2) + zage + zfai + zparty + ztime|subject)+
    (1|observationID),
  data = m6d,
  family = categorical(link="logit"),
  prior = m6strprior,
  chains = 3, 
  cores = 9,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 10000, 
  warmup = 5000, 
  thin = 2,
  control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
# Time difference of 6.247682 hours

plot(m6_str)
summary(m6_str)

# now the ppc plot
# we need a plot with four columns per scan type
# two columns will be a 0 and 1 for the observed data for that scan
# other column will be the mode of the posterior draws for a given observation, with 0/1 scan type

pred_vals = predict(m6_str, summary = FALSE)
predictions = as.data.frame(t(pred_vals))
pred_vals2 <- as.data.frame(pred_vals)

smode = mode(predictions)
range(smode)
colnames(smode) <- c("predictions")
sort(unique(m6d$scan))

ppc1d <- m6d
ppc1d$predictions <- smode$predictions
head(ppc1d)

xx <- aggregate(ppc1d$subject, by=list(ppc1d$scan),length)
yy <- aggregate(ppc1d$subject, by=list(ppc1d$predictions),length)

xx$mode_preds <- yy$x
colnames(xx) <- c("scan", "observed", "predicted")

ppc_bard <- xx
rm(xx,yy)

pred_vals2$non_social = rowSums(pred_vals2 == "1")
pred_vals2$foodshare = rowSums(pred_vals2 == "2")
pred_vals2$nurse = rowSums(pred_vals2 == "3")
pred_vals2$socioff = rowSums(pred_vals2 == "4")
pred_vals2$socioth = rowSums(pred_vals2 == "5")

quantile(pred_vals2$non_social,prob=c(0.0,1.0))
# 2757 3081
quantile(pred_vals2$foodshare,prob=c(0.0,1.0))
# 112 249
quantile(pred_vals2$nurse,prob=c(0.0,1.0))
# 359 562
quantile(pred_vals2$socioff,prob=c(0.0,1.0))
# 460 667
quantile(pred_vals2$socioth,prob=c(0.0,1.0))
# 241 432

ppc_bard$lower <- c(2757,112,359,460,241)
ppc_bard$upper <- c(3081,249,562,667,432)

# reformat to make easier for plotting

xx <- ppc_bard$observed
yy <- ppc_bard$predicted
zz <- c(xx,yy)
xx <- rep("observed",5)
yy <- rep("posterior",5)
zzz <- c(xx,yy)
xxx <- c("Non social", "Foodshare with offspring",
         "Nurse", "Sociopositive with offspring",
         "Sociopositive with others")
xxx <- c(xxx,xxx)

test <- as.data.frame(zzz)
test$scan <- xxx
test$count <- zz
lower <- ppc_bard$lower
lower <- c(NA,NA,NA,NA,NA,lower)
test$lower <- lower
upper <- ppc_bard$upper
upper <- c(NA,NA,NA,NA,NA,upper)
test$upper <- upper

ppc_bard <- test
rm(test,lower,upper,xx,xxx,yy,zz,zzz)
ppc_bard$y <- ppc_bard$zzz

ppc_bard$y <- factor(ppc_bard$y, levels=c("posterior","observed"))
ppc_bard$scan <- factor(ppc_bard$scan, 
                        levels=c("Non social","Foodshare with offspring",
                                 "Nurse","Sociopositive with offspring",
                                 "Sociopositive with others"))

m6str_ppc1 <- ggplot(ppc_bard, 
                     aes(x=scan, y=count, fill=y)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 13))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))

m6str_ppc1 <- m4str_ppc1 + scale_fill_brewer(palette="Paired") + 
  theme_classic()

m6str_ppc1 <- m6str_ppc1+labs(title="M6:  Strongly regularising priors",x="Social choice",
                              y="Count")
m6str_ppc1 <- m6str_ppc1+theme(plot.title = element_text(face="bold", hjust = 0.5))
