############################################
## R Code to Perform Sample Size Examples from Hulley et al
##
## See Chapter 6 - Sample Size Estimation - Hulley et al (download from following URL)
## https://edisciplinas.usp.br/pluginfile.php/5486505/mod_resource/content/1/Stephen%20B.%20Hulley%2C%20Steven%20R.%20Cummings%2C%20Warren%20S.%20Browner%2C%20Deborah%20G.%20Grady%2C%20Thomas%20B.%20Newm.pdf
##
## Author: Christopher Meaney
## Date: July 2023
############################################

#########
## Dependency packages
#########

## MASS --- for simulating MVN random vectors
library(MASS)


###############################
## Example 1 - Sample size for study design comparing two continuous means 
##
## Parameters/Inputs
## - alpha: desired type-1 error rate (typically 0.05)
## - beta:  desired power (typically >= 0.80)
## - q1/q0: allocation ratios (typically assume 50:50)
## - mu0:   mean of outcome in group 0
## - mu1:   mean of outcome in group 1
## - std:   pooled standard deviation
##
## Formula from URL: https://sample-size.net/sample-size-means/
##
###############################

## UDF for sample size for study design comparing two continuous means
sampsize_2means <- function(alpha, beta, q1, q0, mu0, mu1, std) {
  ## Allocation ratio quantity
  A <- 1/q1 + 1/q0
  ## Numerator quantity (power & type-1 errors)
  B <- (qnorm(alpha/2) + qnorm(beta))^2
  ## (Absolute) Mean difference
  E <- abs(mu0 - mu1)
  ## Sample size
  n <- (A*B)/((E/std)^2)
  ## Return Sample Size to User
  return(ceiling(n/2))
}

## Note: Above UDF computes sample size PER GROUP

## Compute sample size from Example 6.1 Hulley et al
sampsize_ex_6_1 <- sampsize_2means(alpha=0.05, beta=0.2, q1=0.5, q0=0.5, mu0=2.0, mu1=1.80, std=1)
sampsize_ex_6_1 



##
## Verify above analytic calculation using simulation
## We use t-test for comparing two groups --- assuming homogeneous variances
## 
sampsize_2means_sim <- function(n, m0, m1, s, alpha) {
	## Simulate data from random normal data generating mechanism
  y0 <- rnorm(n, mean=m0, sd=s)
	y1 <- rnorm(n, mean=m1, sd=s)
	## Analyze independent random normal samples (y0, y1) using t-test
	t_test <- t.test(x=y0, y=y1, var.equal=TRUE)
	t_test_pval <- t_test$p.value
	t_test_tag <- ifelse(t_test_pval < alpha, 1, 0)
	return(t_test_tag)
}

power_ex_6_1 <- mean(replicate(10000, sampsize_2means_sim(n=sampsize_ex_6_1, m0=1.8, m1=2.0, s=1, alpha=0.05)))
power_ex_6_1

##
## Simulated power matches analytic power almost exactly 
## 








###############################
## Example 2 - Sample size for study design comparing two binomial proportions
##
## Parameters/Inputs
## - alpha: desired type-1 error rate (typically 0.05)
## - beta:  desired power (typically >= 0.80)
## - q1/q0: allocation ratios (typically assume 50:50)
## - p0:    proportion of outcome in group 0
## - p1:    proportion of outcome in group 1
##
## Formula from URL: https://sample-size.net/sample-size-proportions/
##
###############################

## UDF for sample size for study design comparing two binomial proportions
##
## Note: this is calculated WITHOUT a continuity correct included
sampsize_2props <- function(alpha, beta, q1, q0, p0, p1) {
  ## Standard normal variates
  za <- abs(qnorm(alpha/2))
  zb <- abs(qnorm(beta))
  ## Pooled proportion
  P <- (q0*p0) + (q1*p1)
  ## Numerator quantity A
  A <- za * sqrt(P * (1-P) * ( (1/q0) + (1/q1) ) )
  ## Numerator quantity B
  B <- zb * sqrt( p1*(1-p1)*(1/q1) + p0*(1-p0)*(1/q0) )
  ## Denominator quantity C
  C <- (p1 - p0)^2
  ## Sample size
  n <- ((A + B)^2)/C
  ## Return Sample Size to User
  return(ceiling(n/2))
}

## Note: Above computes sample size PER GROUP

## Compute sample size from Example 6.2 Hulley et al
sampsize_ex_6_2 <- sampsize_2props(alpha=0.05, beta=0.2, q1=0.5, q0=0.5, p0=0.2, p1=0.3)
sampsize_ex_6_2 



##
## Verify above analytic calculation using simulation
## We use binomial z-test test for comparing two groups --- assuming no continuity correction
## 
sampsize_2props_sim <- function(n, p0, p1, alpha) {
  ## Simulate data from random binomial data generating mechanism
  y0 <- rbinom(n, size=1, prob=p0)
  y1 <- rbinom(n, size=1, prob=p1)
  ## Analyze independent random binomial samples (y0, y1) using binomial Z-test for independent proportions (no continuity correction)
  n_y0 <- length(y0)
  n_y1 <- length(y1)
  x_y0 <- sum(y0)
  x_y1 <- sum(y1)
  n_vec <- c(n_y0, n_y1)
  x_vec <- c(x_y0, x_y1)
  prop_test <- prop.test(x_vec, n_vec, correct=FALSE)
  prop_test_pval <- prop_test$p.value
  prop_test_tag <- ifelse(prop_test_pval < alpha, 1, 0)
  return(prop_test_tag)
}

power_ex_6_2 <- mean(replicate(10000, sampsize_2props_sim(n=sampsize_ex_6_2, p0=0.3, p1=0.2, alpha=0.05)))
power_ex_6_2

##
## Simulated power matches analytic power almost exactly 
## 








###############################
## Example 3 - Sample size for 1-sided test about correlation coefficient equaling zero or not
##
## Parameters/Inputs
## - alpha: desired type-1 error rate (typically 0.05)
## - beta:  desired power (typically >= 0.80)
## - r:     population level correlation coefficient
##
## Formula from URL: https://sample-size.net/correlation-sample-size/
##
###############################


## UDF for sample size for study design testing whether correlation equal zero (NULL) vs. correlation is non-zero (Alternative)
sampsize_corr <- function(alpha, beta, r) {
  ## Standard normal variates
  za <- abs(qnorm(alpha/2))
  zb <- abs(qnorm(beta))
  ## Numerator quantity C
  C <- 0.5 * log( (1+r)/(1-r) )
  ## Sample size
  n <- ((za + zb)/C)^2 + 3
  ## Return Sample Size to User
  return(ceiling(n))
}

## Note: Above computes overall sample size (since this is one-group design)

## Compute sample size from Example 6.3 Hulley et al
sampsize_ex_6_3 <- sampsize_corr(alpha=0.05, beta=0.1, r=0.3)
sampsize_ex_6_3 


##
## Verify above analytic calculation using simulation
## We use cor.test() function for testing whether correlation is zero/non-zero
## 
sampsize_corr_sim <- function(n, r, alpha) {
  ## Simulate data from MVN data generating mechanism
  mu_vec <- c(0,0)
  sigma_mat <- matrix(c(1,r,r,1), ncol=2, nrow=2, byrow=TRUE)
  y_mat <- mvrnorm(n, mu=mu_vec, Sigma=sigma_mat)
  y0 <- y_mat[,1]
  y1 <- y_mat[,2]
  ## cor.test() to test whether H0: rho=0; vs. HA: rho!=0
  cor_test <- cor.test(y0, y1)
  cor_test_pval <- cor_test$p.value
  cor_test_tag <- ifelse(cor_test_pval < alpha, 1, 0)
  return(cor_test_tag)
}

power_ex_6_3 <- mean(replicate(10000, sampsize_corr_sim(n=sampsize_ex_6_3, r=0.3, alpha=0.05)))
power_ex_6_3

##
## Simulated power matches analytic power almost exactly 
## 









###############################
## Example 4 - Sample size for descriptive (1-sample) design about a single continuous mean
##
## Parameters/Inputs
## - alpha: desired width of (1-aplha)-level confidence interval (typically 0.05)
## - mu:    population mean
## - std:   population standard deviation
## - w:     total width of confidence interval
##
## Formula from URL: https://sample-size.net/sample-size-conf-interval-mean/
##
###############################


## UDF for sample size for study design for precision about a single continuous mean
sampsize_1mean <- function(alpha, std, w) {
  ## Standard normal quantity
  za <- abs(qnorm(alpha/2))
  ## Sample size
  n <- (4 * (za^2) * (std^2))/(w^2)
  ## Return Sample Size to User
  return(ceiling(n))
}

## Compute sample size from Example 6.4 Hulley et al
sampsize_ex_6_4 <- sampsize_1mean(alpha=0.05, std=1, w=0.6)
sampsize_ex_6_4 

##
## Verify above analytic calculation using simulation
## We use t.test() function for computing confidence interval for continuous mean
## 
sampsize_1mean_sim <- function(n, std, alpha) {
  ## Simulate data from random normal data generating mechanism
  y <- rnorm(n, mean=0, sd=std)
  ## t.test() to compute CI
  t_test <- t.test(y, conf.level=(1-alpha))
  t_test_ci <- t_test$conf.int
  return(t_test_ci)
}

precision_ex_6_4 <- mean(apply(replicate(10000, sampsize_1mean_sim(n=sampsize_ex_6_4, std=1, alpha=0.05)), 2, function(x) x[2] - x[1]))
precision_ex_6_4

##
## Simulated coverage (CI width) matches that expected analytic width very closely 
## 









###############################
## Example 5 - Sample size for descriptive (1-sample) design about a single binomial proportion
##
## Parameters/Inputs
## - alpha: desired width of (1-aplha)-level confidence interval (typically 0.05)
## - p:     population proportion
## - w:     total width of confidence interval
##
## Formula from URL: https://sample-size.net/sample-size-conf-interval-proportion/
##
###############################


## UDF for sample size for study design for precision about a single binomial proportion
sampsize_1prop <- function(alpha, p, w) {
  ## Standard normal quantity
  za <- abs(qnorm(alpha/2))
  ## Sample size
  n <- (4 * (za^2) * (p*(1-p)))/(w^2)
  ## Return Sample Size to User
  return(ceiling(n))
}


## Compute sample size from Example 6.5 Hulley et al
sampsize_ex_6_5 <- sampsize_1prop(alpha=0.05, p=0.8, w=0.1)
sampsize_ex_6_5 

##
## Verify above analytic calculation using simulation
## We use prop.test() function for computing confidence interval about a single binomial proportion
## 
sampsize_1prop_sim <- function(n, p, alpha) {
  ## Simulate data from random binomial data generating mechanism
  y <- rbinom(n, size=1, prob=p)
  ## t.test() to compute CI
  binom_test <- binom.test(sum(y), length(y), conf.level=(1-alpha))
  binom_test_ci <- binom_test$conf.int
  return(binom_test_ci)
}

precision_ex_6_5 <- mean(apply(replicate(10000, sampsize_1prop_sim(n=sampsize_ex_6_5, p=0.8, alpha=0.05)), 2, function(x) x[2] - x[1]))
precision_ex_6_5

##
## Simulated coverage (CI width) matches that expected analytic width very closely 
## 



##
## Solve for expected width instead of n (since n is fixed by design)
##
prop1_width_vs_np <- function(n, p, alpha) {
  ## Standard normal quantity
  za <- abs(qnorm(alpha/2))
  ## Width given n, p, za
  w <- sqrt( (4 * (za^2) * (p*(1-p)))/n)
  ## Return width to user
  return(w)
}


n <- c(10, 20)
p <- seq(0.05, 0.95, 0.05)
alpha <- 0.05

df_grid <- data.frame(expand.grid(n,p,alpha))
names(df_grid) <- c("n","p","alpha")
dim(df_grid)

w_list <- list()
for (i in 1:nrow(df_grid)) {
  w_list[i] <- prop1_width_vs_np(n=df_grid$n[i], p=df_grid$p[i], alpha=df_grid$alpha[i])
}

df_grid$w <- unlist(w_list)
df_grid <- with(df_grid, df_grid[order(n,p,alpha),])
df_grid









