install.packages("magrittr")
install.packages("dplyr")
install.packages("RMTstat")
library(RMTstat)
library(magrittr)
library(dplyr)
library(tibble)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(car)

## LOAD DATA ##

dir <- "/Users/brydonbrancart/Desktop/stats/final_project/welders_full.csv"
dat <- read.csv(dir, row.names=1, header=TRUE)

## RUN T-TEST ##

# Define subject/control split
subject = 1:9; control = 10:16

#Run t-test
t.stat = apply(X=dat, MARGIN=1, FUN=function(X){t.test(X[subject],X[control])$statistic})
df = 14
m = length(t.stat)
P = 1 - pt(t.stat, df=14)

## PRODUCE HISTOGRAM AND QQPLOT ##

par(mfrow = c(2, 2)) # Create a 2 x 2 plotting matrix
t.hist <- hist(t.stat, 100, probability = T)
tt <- t.hist$mids
lines(tt, dt(tt, df = df), col='red', lwd=2)
# Check normality of original data all together
X = t(scale(t(dat)))
#X <- X[1:2,]
qqPlot(as.vector(X)[sample(1000)]) # left skew

## CHECK FOR CORRELATION ##

#Sample 1000 to check correlation
n <- ncol(dat)
m <- 1000
J <- sample(1:nrow(dat), m)
X <- t(dat[J,])

R.hat <- cor(X)
#heatmap(R.hat, Rowv = NA, Colv = NA, scale = 'none')
offdiag <- function(R) {
  diag(R) <- NA
  return(R)
}

R.hat.offdiag <- offdiag(R.hat)
hist(R.hat.offdiag, 100, freq = F)
hist(sqrt(n-2) * R.hat.offdiag, 100, freq = F)
curve(dnorm(x), from = -8, to = 8, col = 'red', lwd = 2, add = T)
mean(R.hat.offdiag^2, na.rm = T)


## EXAMINE EIGENVALUES ##

"From this we see that there are eigenvalues that diverge from the null."

#Compare eigen values to null
R.hat.eig <- eigen(R.hat)
plot(R.hat.eig$values); abline(h = 1, col='black', lwd=2)
hist(R.hat.eig$values, breaks = 20, freq = F)
curve(dmp(x, n, m), from = (1 - sqrt(m/n))^2, to = max(R.hat.eig$values), col = 'red', lwd = 2, add = T)
plot(ecdf(R.hat.eig$values))
curve(pmp(x, n, m), from = (1 - sqrt(m/n))^2, to = max(R.hat.eig$values), col = 'red', lwd = 2, add = T)

#Confidence intervals
(1 - sqrt(m/n))^2
(1 + sqrt(m/n))^2

## EXPLORE CORRELATION ##

# Attempt regressing out the first principle component

#Define first principle component via defining svd
R.hat.eig <- eigen(R.hat)
plot(R.hat.eig$values)
abline(h = 1)
X.svd <- svd(X)
lines(X.svd$d^2/n, col='red', lwd=2)

# Regress out first principal component
X.res <- apply(X, 2, function(x) lm(x ~ X.svd$u[,1] - 1)$residuals)
R.res <- cor(X.res)
#heatmap(R.res, Rowv = NA, Colv = NA, scale = 'none')
R.res.offdiag <- offdiag(R.res)
hist(sqrt(n-2) * R.res.offdiag, 100, freq = F, ylim= c(0,.5))
curve(dnorm(x), from = -8, to = 8, col = 'red', lwd = 2, add = T)

# Check eigenvalues again
R.res.eig <- eigen(R.res)
plot(R.res.eig$values); abline(h=1)
hist(R.res.eig$values, breaks = 20, freq = F)
curve(dmp(x, n, m), from = (1 - sqrt(m/n))^2, to = (1 + sqrt(m/n))^2, col = 'red', lwd = 2, add = T)
(1 - sqrt(m/n))^2
(1 + sqrt(m/n))^2
# There are still some hidden covariates left

# Regress out two principal components
X.res <- apply(X, 2, function(x) lm(x ~ X.svd$u[,1] - 1)$residuals)
R.res <- cor(X.res)
#heatmap(R.res, Rowv = NA, Colv = NA, scale = 'none')
R.res.offdiag <- offdiag(R.res)
hist(sqrt(n-2) * R.res.offdiag, 100, freq = F, ylim=c(0,0.5))
curve(dnorm(x), from = -8, to = 8, col = 'red', lwd = 2, add = T)

# Check eigenvalues again
R.res.eig <- eigen(R.res)
plot(R.res.eig$values); abline(h=1)
hist(R.res.eig$values, breaks = 20, freq = F)
curve(dmp(x, n, m), from = (1 - sqrt(m/n))^2, to = (1 + sqrt(m/n))^2, col = 'red', lwd = 2, add = T)
(1 - sqrt(m/n))^2
(1 + sqrt(m/n))^2

#### REGRESS OUT CORRELATION ####

# Do again for decorrelated genes
dat.res <- apply(t(dat), 2, function(x) lm(x ~ X.svd$u[,1] - 1)$residuals)
x <- c(rep(1, length(subject)), rep(0, length(control)))
t.stat.res = apply(X=t(dat.res), MARGIN=1, FUN=function(X){summary(lm(unlist(X) ~ x))$coef[2,3]})
plot(t.stat, t.stat.res)
hist(t.stat.res, 100, freq = F, ylim = c(0,0.5))
curve(dnorm(x), from = -5, to = 5, col = 'red', lwd = 2, add = T)

#### PERFORM ERROR CORRECTION ####
# FDR analysis
p.values.res <- 2*(1 - pt(abs(t.stat.res), df = df))
plot(ecdf(p.values.res))
abline(0, 1, col = 'blue')
p.res <- sort(p.values.res)
alpha = 0.05
plot(p.res[1:10]); abline(0, alpha/m, col = 'red')

FDRhat.res = p.res / pmax(ecdf(p.res)(p.res),1/m)
plot(p.res, FDRhat.res, type='l', lwd=2)
abline(h = alpha, lwd=2, col='red')
p.res[2]
