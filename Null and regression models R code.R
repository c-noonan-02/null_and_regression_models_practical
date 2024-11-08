# Null and regression models in ecology

rm(list=ls())
getwd()
setwd("C:/Users/charl/Documents/University/Year 5 Masters/Semester 1/Biodiversity Under Pressure/8. Week Eight/Week Eight Computer Practical 1/Data")
getwd()

library(vegan)
library(gam)
library(scam)
library(gdm)
library(car)
library(betapart)
library(tidyverse)
library(zetadiv)


# 1. NULL MODELS: PERMUTATION ALGORITHMS

# load the lists for society
load("islands_Soc_Haw_null_models_practical.RData")
society_data_pa <- dat.soc.pa

# compute the beta diversity values for the two archipelagos
# apply function beta.pair()
# this will generate a list with the Sorensen, Simpson and nestedness beta diversity values
society_beta_obs <- c((beta.pair(society_data_pa)$beta.sor), mean(beta.pair(society_data_pa)$beta.sim))
# compute the average beta diversity value across islands for each index
society_beta_obs <- c(mean(beta.pair(society_data_pa)$beta.sor), mean(beta.pair(society_data_pa)$beta.sim))

# now need to assess if these values are different from what we would expect by chance
# we will use a permutation algorithm to generate beta diversity values under randomness
?permatfull

# impose no constraints
society_data_pa_perm_none <- permatfull(m=society_data_pa, fixedmar = "none", times = 99, mtype = "prab")
# fix the row sums - keeps alpha diversity the same
society_data_pa_perm_row <- permatfull(m=society_data_pa,fixedmar = "rows",times=99,mtype="prab")
# fix the column sums - keeps the occupancy (number sites) the same
society_data_pa_perm_col <- permatfull(m=society_data_pa,fixedmar = "columns",times=99,mtype="prab")
# fix both
society_data_pa_perm_both <- permatfull(m=society_data_pa,fixedmar = "both",times=99,mtype="prab")

# we now need to compute beta diversity values for each permutated matrix
# then compute the average for the beta diversity values corresponding to all pairs of islands

?lapply # returns a list of the same length as X, each element of which is the result of applying FUN to the corresponding element of X
?unlist # simplify a list to produce a vector containing all atomic components that occur in x
?matrix # creates a matrix from the given set of values

# create the function
beta_mean <- function(dat){
  beta.dat <- beta.pair(dat)
  return(c(mean(beta.dat$beta.sor),mean(beta.dat$beta.sim)))
}

# use the function
society_beta_rand_none <- data.frame(matrix(unlist(lapply(society_data_pa_perm_none$perm,beta_mean)),99,2,byrow = TRUE))
names(society_beta_rand_none) <- c("Sorensen","Simpson")

society_beta_rand_row <- data.frame(matrix(unlist(lapply(society_data_pa_perm_row$perm,beta_mean)),99,2,byrow = TRUE))
names(society_beta_rand_row) <- c("Sorensen","Simpson")

society_beta_rand_col <- data.frame(matrix(unlist(lapply(society_data_pa_perm_col$perm,beta_mean)),99,2,byrow = TRUE))
names(society_beta_rand_col) <- c("Sorensen","Simpson")

society_beta_rand_both <- data.frame(matrix(unlist(lapply(society_data_pa_perm_both$perm,beta_mean)),99,2,byrow = TRUE))
names(society_beta_rand_both) <- c("Sorensen","Simpson")

# now we can plot the distributions under the different permutation algorithms
# along with the observed beta diversity values for the Sorensen and Simpson dissimilarity index
# and the 5% and 95% quantiles

par(mfrow=c(2,4))
hist(society_beta_rand_none$Sorensen,breaks=seq(0.5,1,0.01),main="None fixed",xlab="Sorensen")
abline(v=quantile(society_beta_rand_none$Sorensen,0.025),col="blue")
abline(v=quantile(society_beta_rand_none$Sorensen,0.975),col="blue")
abline(v=society_beta_obs[1],col="red")
hist(society_beta_rand_row$Sorensen,breaks=seq(0.5,1,0.01),main="Rows fixed",xlab="Sorensen")
abline(v=quantile(society_beta_rand_row$Sorensen,0.025),col="blue")
abline(v=quantile(society_beta_rand_row$Sorensen,0.975),col="blue")
abline(v=society_beta_obs[1],col="red")
hist(society_beta_rand_col$Sorensen,breaks=seq(0.5,1,0.01),main="Columns fixed",xlab="Sorensen")
abline(v=quantile(society_beta_rand_col$Sorensen,0.025),col="blue")
abline(v=quantile(society_beta_rand_col$Sorensen,0.975),col="blue")
abline(v=society_beta_obs[1],col="red")
hist(society_beta_rand_both$Sorensen,breaks=seq(0.5,1,0.01),main="Both fixed",xlab="Sorensen")
abline(v=quantile(society_beta_rand_both$Sorensen,0.025),col="blue")
abline(v=quantile(society_beta_rand_both$Sorensen,0.975),col="blue")
abline(v=society_beta_obs[1],col="red")

hist(society_beta_rand_none$Simpson,breaks=seq(0.2,1,0.01),main="",xlab="Simpson")
abline(v=quantile(society_beta_rand_none$Simpson,0.025),col="blue")
abline(v=quantile(society_beta_rand_none$Simpson,0.975),col="blue")
abline(v=society_beta_obs[2],col="red")
hist(society_beta_rand_row$Simpson,breaks=seq(0.2,1,0.01),main="",xlab="Simpson")
abline(v=quantile(society_beta_rand_row$Simpson,0.025),col="blue")
abline(v=quantile(society_beta_rand_row$Simpson,0.975),col="blue")
abline(v=society_beta_obs[2],col="red")
hist(society_beta_rand_col$Simpson,breaks=seq(0.2,1,0.01),main="",xlab="Simpson")
abline(v=quantile(society_beta_rand_col$Simpson,0.025),col="blue")
abline(v=quantile(society_beta_rand_col$Simpson,0.975),col="blue")
abline(v=society_beta_obs[2],col="red")
hist(society_beta_rand_both$Simpson,breaks=seq(0.2,1,0.01),main="",xlab="Simpson")
abline(v=quantile(society_beta_rand_both$Simpson,0.025),col="blue")
abline(v=quantile(society_beta_rand_both$Simpson,0.975),col="blue")
abline(v=society_beta_obs[2],col="red")

# None fixed and Rows fixed overestimated the Sorensen and Simpson
# Columns fixed underestimated for Sorensen and overestimated for Simpson
# Both fixed gives the closest estimate to the observed values for both

# now do this all for the hawaii data set as well
# load the lists for Hawaii
load("permutations_hawaii.RData")
hawaii_data_pa <- dat.haw.pa

hawaii_data_pa_perm_none <- permatfull(m=hawaii_data_pa,fixedmar = "none",times=99,mtype="prab")
hawaii_data_pa_perm_row <- permatfull(m=hawaii_data_pa,fixedmar = "rows",times=99,mtype="prab")
hawaii_data_pa.perm_col <- permatfull(m=hawaii_data_pa,fixedmar = "columns",times=99,mtype="prab")
hawaii_data_pa_perm_both <- permatfull(m=hawaii_data_pa,fixedmar = "both",times=99,mtype="prab")

# use the function
hawaii_beta_rand_none <- data.frame(matrix(unlist(lapply(hawaii_data_pa_perm_none$perm,beta_mean)),99,2,byrow = TRUE))
names(hawaii_beta_rand_none) <- c("Sorensen","Simpson")
hawaii_beta_rand_row <- data.frame(matrix(unlist(lapply(hawaii_data_pa_perm_row$perm,beta_mean)),99,2,byrow = TRUE))
names(hawaii_beta_rand_row) <- c("Sorensen","Simpson")
hawaii_beta_rand_col <- data.frame(matrix(unlist(lapply(dat.haw.pa.perm_col$perm,beta_mean)),99,2,byrow = TRUE))
names(hawaii_beta_rand_col) <- c("Sorensen","Simpson")
hawaii_beta_rand_both <- data.frame(matrix(unlist(lapply(hawaii_data_pa_perm_both$perm,beta_mean)),99,2,byrow = TRUE))
names(hawaii_beta_rand_both) <- c("Sorensen","Simpson")

# plot
par(mfrow=c(2,4))
hist(hawaii_beta_rand_none$Sorensen,breaks=seq(0,1,0.025),main="None fixed",xlab="Sorensen")
abline(v=quantile(hawaii_beta_rand_none$Sorensen,0.025),col="blue")
abline(v=quantile(hawaii_beta_rand_none$Sorensen,0.975),col="blue")
abline(v=hawaii_beta_obs[1],col="red")
hist(hawaii_beta_rand_row$Sorensen,breaks=seq(0,1,0.025),main="Rows fixed",xlab="Sorensen")
abline(v=quantile(hawaii_beta_rand_row$Sorensen,0.025),col="blue")
abline(v=quantile(hawaii_beta_rand_row$Sorensen,0.975),col="blue")
abline(v=hawaii_beta_obs[1],col="red")
hist(hawaii_beta_rand_col$Sorensen,breaks=seq(0,1,0.025),main="Columns fixed",xlab="Sorensen")
abline(v=quantile(hawaii_beta_rand_col$Sorensen,0.025),col="blue")
abline(v=quantile(hawaii_beta_rand_col$Sorensen,0.975),col="blue")
abline(v=hawaii_beta_obs[1],col="red")
hist(hawaii_beta_rand_both$Sorensen,breaks=seq(0,1,0.025),main="Both fixed",xlab="Sorensen")
abline(v=quantile(hawaii_beta_rand_both$Sorensen,0.025),col="blue")
abline(v=quantile(hawaii_beta_rand_both$Sorensen,0.975),col="blue")
abline(v=hawaii_beta_obs[1],col="red")

hist(hawaii_beta_rand_none$Simpson,breaks=seq(0,1,0.025),main="",xlab="Simpson")
abline(v=quantile(hawaii_beta_rand_none$Simpson,0.025),col="blue")
abline(v=quantile(hawaii_beta_rand_none$Simpson,0.975),col="blue")
abline(v=hawaii_beta_obs[2],col="red")
hist(hawaii_beta_rand_row$Simpson,breaks=seq(0,1,0.025),main="",xlab="Simpson")
abline(v=quantile(hawaii_beta_rand_row$Simpson,0.025),col="blue")
abline(v=quantile(hawaii_beta_rand_row$Simpson,0.975),col="blue")
abline(v=hawaii_beta_obs[2],col="red")
hist(hawaii_beta_rand_col$Simpson,breaks=seq(0,1,0.025),main="",xlab="Simpson")
abline(v=quantile(hawaii_beta_rand_col$Simpson,0.025),col="blue")
abline(v=quantile(hawaii_beta_rand_col$Simpson,0.975),col="blue")
abline(v=hawaii_beta_obs[2],col="red")
hist(hawaii_beta_rand_both$Simpson,breaks=seq(0,1,0.025),main="",xlab="Simpson")
abline(v=quantile(hawaii_beta_rand_both$Simpson,0.025),col="blue")
abline(v=quantile(hawaii_beta_rand_both$Simpson,0.975),col="blue")
abline(v=hawaii_beta_obs[2],col="red")


# 2. REGRESSION MODELS ON SPECIES RICHNESS
