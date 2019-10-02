########################################################################################################################################################## 
##############################          Code used in "Is the R coefficient of interest in cluster randomized trials         ############################## 
##############################                               with a binary outcome?" Mbekwe et al.                          ############################## 
########################################################################################################################################################## 


##############################        Simulation to evaluate Crespi et al.'s R coefficient dependence on prevalence         ############################## 


##########     est_FC: a function to estimate the intraclass correlation coefficient using the Fleiss-Cuzick estimator     ##########

#####   Arguments   #####
# g: a group variable
# out: the outcome
# data: a dataframe containing g and out

#####   Value   #####
# res: the estimated intraclass correlation coefficient

est_FC <- function(g, out, data)
{
  N <- nrow(data)
  
  l <- list(g=substitute(g), out=substitute(out))
  d <- data.frame(g=eval(expr = l$g, envir = data), out=eval(expr = l$out, envir = data))
  
  k <- length(levels(d$g))
  
  size <- rowSums(table(d$g,d$out))

  zi <- table(d$g,d$out)[, 2]
  
  pchap <- sum(zi)/N
  
  res <- 1 - sum((zi*(size-zi))/size)/((N - k)*pchap*(1 - pchap))
  
  return(res)
}



##########    sum_f_int: a function to calculate the binary intraclass correlation coefficient     ##########
##########               knowing the continuous intraclass correlation coefficient                 ########## 

f_int <- function(x, h)
{
  (1/sqrt(1-x^2))*exp(-h^2/(1+x))
}

#####   Arguments   #####
# r_tet: the continuous intraclass correlation coefficient
# n: the number of intervals
# h: the threshold of dichotomization
# p: the prevalence

#####   Value   #####
# res: the estimated binary intraclass correlation coefficient

sum_f_int <- function(r_tet, n, h, p)
{ 
  x0 <- 0; s <- 0
  delta <- (r_tet-0)/n
  
  for (i in 1:(n-1))
  {
    x0 <- x0 + delta
    s <- s + f_int(x=x0, h)
  }
  
  s <- delta*((f_int(x=0, h) + f_int(x=r_tet, h))/2 + s)
  
  res <- (1/(2*pi*p*(1-p)))*s
  
  return(res)
}



##########     simul: a function to simulate correlated binary data and estimate the intraclass correlation and the R coefficients     ##########

#####   Arguments   #####
# k: the number of clusters
# m: the average size of clusters
# v: the variance of cluster sizes
# rho: the theoretical intraclass correlation coefficient
# p: the prevalence

#####   Values   #####
# Prop_success: the estimated success proportion
# Rho: the estimated intraclass correlation coefficient
# Coef_R: the estimated R coefficient

simul <-  function(k, m, v, rho, p)
{
  if (v == 0)  {clust_siz <- rep(m,k)}
  else {
    clust_siz <- rnbinom(k, size=m^2/(v-m), mu=m)
    i <- 0
    while (is.element(0, clust_siz) & i < 1000)
    {
      clust_siz <- rnbinom(k, size=m^2/(v-m), mu=m)
      i <- i + 1
    }
  }
  n <- sum(clust_siz)
  
  # simulation of correlated binary outcome using the method of Lunn and Davies
  # model: Xij = (1-Uij)Yij + UijZi
  # Zi ~ Binom(1,pi); Yij ~ Binom(1,pi); Uij ~ Binom(1,sqrt(rho))
  
  Z <- rep(rbinom(k, 1, p), clust_siz)
  Y <- rbinom(n, 1, p)
  U <- rbinom(n, 1, sqrt(rho))
  
  X <- (1-U)*Y + U*Z
  
  while(length(unique(X)) == 1)
  {
    Z <- rep(rbinom(k, 1, p), clust_siz)
    Y <- rbinom(n, 1, p)
    U <- rbinom(n, 1, sqrt(rho))
    
    X <- (1-U)*Y + U*Z
  }
  
  est_p <- as.numeric(prop.table(table(X))[2])
  
  clust <- factor(rep(1:k, clust_siz))
  
  don <- data.frame(Outcome=X, Cluster=clust)
  
  # estimation of rho using the Fleiss-Cuzick estimator
  est_rho <- est_FC(g=Cluster, out=Outcome, data=don)
  
  if (est_rho < 0) {est_rho <- 0}
  
  # calculation of R
  r <- 1 + est_rho*(1-est_p)/est_p
  
  return(list(Prop_success=est_p, Rho=est_rho, Coef_R=r))
}



##########     vecplot_simul: a function to replicate simul     ##########

vecplot_simul <- function(nsimul, k0, m0, v0, rho0, p0)
{
  res <- replicate(nsimul, simul(k=k0, m=m0, v=v0, rho=rho0, p=p0), simplify = F)
  
  prop <- mean(sapply(res, function(a){a$Prop_success}))
  rho <- mean(sapply(res, function(a){a$Rho}))
  r <- mean(sapply(res, function(a){a$Coef_R}))
  
  return(list(prop=prop, rho=rho, r=r))
}


vecplot_simul(nsimul=100, k0=20, m0=25, v0=225, rho0=sum_f_int(r_tet = 0.01, n=100, h=qnorm(1-0.01), p=0.01), p0 = 0.01)
vecplot_simul(nsimul=100, k0=20, m0=25, v0=225, rho0=sum_f_int(r_tet = 0.01, n=100, h=qnorm(1-0.5), p=0.5), p0 = 0.5)
vecplot_simul(nsimul=100, k0=20, m0=25, v0=225, rho0=sum_f_int(r_tet = 0.01, n=100, h=qnorm(1-0.99), p=0.99), p0 = 0.99)





##############################          Power comparison using the three approaches proposed by Crespi et al.          ############################## 



##########    sum_f_int: a function to calculate the continuous intraclass correlation coefficient     ##########
##########                    knowing the binary intraclass correlation coefficient                    ########## 

####   Arguments   ####
# rhobin: the binary intraclass correlation coefficient
# n: the number of intervals
# h: the threshold of dichotomization
# p: the prevalence

#####   Value   #####
# res: the estimated continuous intraclass correlation coefficient

rec_sum_f_int <- function(rhobin, n, h, p)
{
  f <- function(x, n, h, p) sum_f_int(x, n, h, p)  - rhobin
  res <- uniroot(f = f, interval = c(0,1), n = n, h = h, p = p)$root
  return(res)
}

##########    sim: a function to simulate binary correlated data    ########## 

#####   Arguments   #####
# k: the number of clusters
# m: the average size of clusters
# v: the variance of cluster sizes
# rho: the theoretical intraclass correlation coefficient
# p: the prevalence

#####   Value   #####
# don: the simulated dataframe

sim  <-  function(k, m, v, rho, p)
{
  if (v == 0)  {clust_siz <- rep(m,k)}
  
  else {
    
    clust_siz <- rnbinom(k, size=m^2/(v-m), mu=m)
    
    i <- 0
    while (is.element(0, clust_siz) & i < 1000)
    {
      clust_siz <- rnbinom(k, size=m^2/(v-m), mu=m)
      i <- i + 1
    }
  }
  
  n <- sum(clust_siz)
  
  Z <- rep(rbinom(k, 1, p), clust_siz)
  Y <- rbinom(n, 1, p)
  U <- rbinom(n, 1, sqrt(rho))
  
  X <- (1-U)*Y + U*Z
  
  while(length(unique(X)) == 1)
  {
    Z <- rep(rbinom(k, 1, p), clust_siz)
    Y <- rbinom(n, 1, p)
    U <- rbinom(n, 1, sqrt(rho))
    
    X <- (1-U)*Y + U*Z
  }
  
  p <- as.numeric(prop.table(table(X))[2])
  
  clust <- factor(rep(1:k, clust_siz))
  
  don <- data.frame(Outcome=X, Cluster=clust)
  
  return(don)
}



##########    ssr, ssicca, ssiccab: functions to calculate sample size following Crespi et al.'s formulas     ##########

#####   Arguments   #####
# pi1: the expected prevalence in arm 1 of the future study 
# pi2: the expected prevalence in arm 2 of the future study 
# m: the fixed cluster sizes
# r1: the observed R coefficient in arm 1 of the previous study
# r2: the observed R coefficient in arm 2 of the previous study
# rho1: the observed intraclass correlation coefficient in arm 1 of the previous study
# rho2: the observed intraclass correlation coefficient in arm 2 of the previous study
# rho: the common intraclass correlation coefficient 

#####   Value   #####
# k: the number of clusters per arm

ssr <- function(pi1, pi2, alpha, beta, m, r1, r2)
{
  k <- ((qnorm(1 - alpha/2) + qnorm(1-beta))^2*(pi1*(1 - pi1 + (m-1)*(r1 - 1)*pi1) + pi2*(1 - pi2 + (m-1)*(r2 - 1)*pi2)))/(m*(pi1-pi2)^2)
  return(k)
}

ssicca <- function(pi1, pi2, alpha, beta, m, rho1, rho2)
{
  k <- ((qnorm(1 - alpha/2) + qnorm(1-beta))^2*(pi1*(1-pi1)*(1 + (m-1)*rho1) + pi2*(1-pi2)*(1 + (m-1)*rho2)))/(m*(pi1-pi2)^2)
  return(k)
}

ssiccb <- function(pi1, pi2, alpha, beta, m, rho)
{
  k <- ((qnorm(1 - alpha/2) + qnorm(1-beta))^2*(pi1*(1-pi1) + pi2*(1-pi2))*(1 + (m-1)*rho))/(m*(pi1-pi2)^2)
  return(k)
}



##########    adj_chisquare : a function to compute the adjusted chi-square test     ##########

#####   Arguments   #####
# outcome: outcome variable
# cluster: cluster variable
# group: arm variable
# data: a dataframe containing outcome, cluster and group

#####   Value   #####
# P.value: the p-value of the test

adj_chisquare <- function(outcome, cluster, group, data)
{
  l <- list(outcome=substitute(outcome), cluster=substitute(cluster), group=substitute(group))
  d <- data.frame(outcome=eval(expr = l$outcome, envir = data), cluster=eval(expr = l$cluster, envir = data), 
                  group=eval(expr = l$group, envir = data))
  
  p <- mean(d$outcome)
  m <- nrow(d)
  
  d1 <- d[d$group==unique(d$group)[1], ] # data corresponding to the first arm
  size1 <- rowSums(table(d1$cluster,d1$outcome)) 
  m1 <- nrow(d1)
  p1 <- mean(d1$outcome)
  
  
  d2 <- d[d$group==unique(d$group)[2], ] # data corresponding to the second arm
  size2 <- rowSums(table(d2$cluster,d2$outcome)) 
  m2 <- nrow(d2)
  p2 <- mean(d2$outcome)
  
  k <- length(unique(d1$cluster)) + length(unique(d2$cluster))
  
  msc <- (sum(size1*(table(d1$cluster,d1$outcome)[, 2]/rowSums(table(d1$cluster,d1$outcome)) - p1)^2) +
            sum(size2*(table(d2$cluster,d2$outcome)[, 2]/rowSums(table(d2$cluster,d2$outcome)) - p2)^2))/(k - 2)
  
  msw <- (sum(size1*(table(d1$cluster,d1$outcome)[, 2]/rowSums(table(d1$cluster,d1$outcome)))*(1 - table(d1$cluster,d1$outcome)[, 2]/rowSums(table(d1$cluster,d1$outcome)))) + 
            sum(size2*(table(d2$cluster,d2$outcome)[, 2]/rowSums(table(d2$cluster,d2$outcome)))*(1 - table(d2$cluster,d2$outcome)[, 2]/rowSums(table(d2$cluster,d2$outcome))))) /(m-k)
  
  m0 <- (m -  (sum(size1^2)/m1 + sum(size2^2)/m2))/(k -2)
  
  rho <- (msc - msw)/(msc + (m0 - 1)*msw)
  
  c1 <- sum(size1*(1 + (size1 - 1)*rho))/sum(size1)
  c2 <- sum(size2*(1 + (size2 - 1)*rho))/sum(size2)
  
  stat <- (m1*(p1 - p)^2)/(c1*p*(1-p)) + (m2*(p2 - p)^2)/(c2*p*(1-p))
  
  p_val <- 1 - pchisq(stat, 1)
  
  return(list(P.value = p_val))
}



##########    rep_simul: a function to compute empirical power     ##########

#####   Arguments   #####
# nsimul: number of replication
# m: average size of clusters
# v: variance of cluster sizes
# p_obs.1: the observed prevalence in arm 1 of the previous study 
# p_obs.2: the observed prevalence in arm 2 of the previous study 
# p_exp.1: the expected prevalence in arm 1 of the future study 
# p_exp.2: the expected prevalence in arm 2 of the future study 
# r1: the observed R coefficient in arm 1 of the previous study
# r2: the observed R coefficient in arm 2 of the previous study
# method: method for sample size calculation, must be one of "R-based", "ICC A" or "ICC B'".

#####   Values   #####
# Power: the estimated power
# k: the required number of clusters per arm

rep_simul <- function(nsimul, m, v, p_obs.1, p_obs.2, p_exp.1, p_exp.2, r1, r2, alpha, beta, method)
{
  res <- c()
  
  if (method == "R-based")
  {
    k1 <- ceiling(ssr(pi1 = p_exp.1, pi2 = p_exp.2, alpha = alpha, beta = beta, m = m, r1 = r1, r2 = r2))
    k2 <- k1
    
    for (i in 1:nsimul)
    {
      b1 <- sim(k = k1, m = m, v = v, 
                rho = sum_f_int(r_tet = rec_sum_f_int(rhobin = ((r1 - 1)*p_obs.1)/(1 - p_obs.1), n = 100, h = qnorm(p_obs.1), p  = p_obs.1),
                                n = 100, h = qnorm(p_exp.1), p = p_exp.1), p = p_exp.1)
      
      b2 <- sim(k = k2, m = m, v = v, 
                rho = sum_f_int(r_tet = rec_sum_f_int(rhobin = ((r2 - 1)*p_obs.2)/(1 - p_obs.2), n = 100, h = qnorm(p_obs.2), p  = p_obs.2),
                                n = 100, h = qnorm(p_exp.2), p = p_exp.2), p = p_exp.2)
      
      b1$Intervention <- factor(rep(x = "A", times = nrow(b1)))
      b2$Intervention  <- factor(rep(x = "B", times = nrow(b2)))
      
      b <- rbind(b1, b2)
      
      res[i] <- ifelse(adj_chisquare(outcome = Outcome, cluster = Cluster, group = Intervention, data = b)$P.value < 0.05, 1, 0)
    }
    
    return(list(Power = 100*round(mean(res), 2), k = k1))
  }
  
  else if (method == "ICC A")
  {
    k1 <- ceiling(ssicca(pi1 = p_exp.1, pi2 = p_exp.2, alpha = alpha, beta = beta, m = m, rho1 = ((r1 - 1)*p_obs.1)/(1 - p_obs.1), rho2 = ((r2 - 1)*p_obs.2)/(1 - p_obs.2)))
    k2 <- k1
    
    for (i in 1:nsimul)
    {
      b1 <- sim(k = k1, m = m, v = v, 
                rho = sum_f_int(r_tet = rec_sum_f_int(rhobin = ((r1 - 1)*p_obs.1)/(1 - p_obs.1), n = 100, h = qnorm(p_obs.1), p  = p_obs.1),
                                n = 100, h = qnorm(p_exp.1), p = p_exp.1), p = p_exp.1)
      
      b2 <- sim(k = k2, m = m, v = v, 
                rho = sum_f_int(r_tet = rec_sum_f_int(rhobin = ((r2 - 1)*p_obs.2)/(1 - p_obs.2), n = 100, h = qnorm(p_obs.2), p  = p_obs.2),
                                n = 100, h = qnorm(p_exp.2), p = p_exp.2), p = p_exp.2)
      
      b1$Intervention <- factor(rep(x = "A", times = nrow(b1)))
      b2$Intervention  <- factor(rep(x = "B", times = nrow(b2)))
      
      b <- rbind(b1, b2)
      
      res[i] <- ifelse(adj_chisquare(outcome = Outcome, cluster = Cluster, group = Intervention, data = b)$P.value < 0.05, 1, 0)
    }
    
    return(list(Power = 100*round(mean(res), 2), k = k1))
  }
  
  else if (method == "ICC B")
  {
    p <- (p_obs.1 + p_obs.2)/2
    k1 <- ceiling(ssiccb(pi1 = p_exp.1, pi2 = p_exp.2, alpha = alpha, beta = beta, m = m, rho = (((r1 - 1)*p_obs.1)/(1 - p_obs.1) + ((r2 - 1)*p_obs.2)/(1 - p_obs.2))/2))
    k2 <- k1
    
    for (i in 1:nsimul)
    {
      b1 <- sim(k = k1, m = m, v = v, 
                rho = sum_f_int(r_tet = rec_sum_f_int(rhobin = ((r1 - 1)*p_obs.1)/(1 - p_obs.1), n = 100, h = qnorm(p_obs.1), p  = p_obs.1),
                                n = 100, h = qnorm(p_exp.1), p = p_exp.1), p = p_exp.1)
      
      b2 <- sim(k = k2, m = m, v = v, 
                rho = sum_f_int(r_tet = rec_sum_f_int(rhobin = ((r2 - 1)*p_obs.2)/(1 - p_obs.2), n = 100, h = qnorm(p_obs.2), p  = p_obs.2),
                                n = 100, h = qnorm(p_exp.2), p = p_exp.2), p = p_exp.2)
      
      b1$Intervention <- factor(rep(x = "A", times = nrow(b1)))
      b2$Intervention  <- factor(rep(x = "B", times = nrow(b2)))
      
      b <- rbind(b1, b2)
      
      res[i] <- ifelse(adj_chisquare(outcome = Outcome, cluster = Cluster, group = Intervention, data = b)$P.value < 0.05, 1, 0)
    }
    
    return(list(Power = 100*round(mean(res), 2), k = k1))
  }
  
  else (print("This method is not implemented"))
}

rep_simul(nsimul = 100, m = 20, v = 0, p_obs.1 = 0.7, p_obs.2 = 0.5, p_exp.1 = 0.5, p_exp.2 = 0.3, r1 = 1.02, r2 = 1.2, alpha = 0.05, beta = 0.2,
          method = "R-based")
rep_simul(nsimul = 100, m = 20, v = 0, p_obs.1 = 0.7, p_obs.2 = 0.5, p_exp.1 = 0.5, p_exp.2 = 0.3, r1 = 1.02, r2 = 1.2, alpha = 0.05, beta = 0.2,
          method = "ICC A")
rep_simul(nsimul = 100, m = 20, v = 0, p_obs.1 = 0.7, p_obs.2 = 0.5, p_exp.1 = 0.5, p_exp.2 = 0.3, r1 = 1.02, r2 = 1.2, alpha = 0.05, beta = 0.2, 
          method = "ICC B")