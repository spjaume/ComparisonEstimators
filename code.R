# load packages needed
library(ThresholdROC)
library(tidyverse)

# suppress warnigns from thres2()
options(warn=-1)

# generation of dataset and threshold estimation
sim_sample <- function(n1, n2, par1.1, par1.2, par2.1, par2.2, rho){
  # n1: sample size of healthy subpopulation
  # n2: sample size of diseased subpopulation
  # par1.1: mean for the healthy subpopulation
  # par1.2: std dev for the healthy subpopulation
  # par2.1: mean for the diseased subpopulation
  # par2.2: std dev for the diseased subpopulation
  # rho: disease prevalence
  
  # healthy/diseased continous diagnostic test
  k1 <- rnorm(n1, par1.1, par1.2)
  k2 <- rnorm(n2, par2.1, par2.2)
  
  # threshold estimation
  t.eq <- thres2(k1, k2, rho, method="eq", ci.method="d") # equal estimator
  t.uneq <- try(thres2(k1, k2, rho, method="uneq", ci.method="d")) # unequal estimator
  
  # extract estimates, standard errors and confidence intervals
  th.eq <- t.eq$T$thres
  th.se.eq <- t.eq$CI$se
  th.ll.eq <- t.eq$CI$lower
  th.ul.eq <- t.eq$CI$upper
  
  
  if (class(t.uneq)!="try-error"){
    th.un <- t.uneq$T$thres
    th.se.un <- t.uneq$CI$se
    th.ll.un <- t.uneq$CI$lower
    th.ul.un <- t.uneq$CI$upper
  }else{
    th.un <- NA
    th.se.un <- NA
    th.ll.un <- NA
    th.ul.un <- NA
  }
  # F test p-value
  varp <- var.test(k1,k2)$p.value
  data.frame(th.eq, th.un, th.se.eq, th.se.un, th.ll.eq, th.ul.eq, th.ll.un, th.ul.un, varp, sd(k1), sd(k2))
}

# this function returns 1 if tv is inside the interval (ll, ul) and 0 otherwise
within_int <- function(ll, ul, tv){
  return((ll<tv)*(ul>tv))
}


sim_var <- function(k, dist1, dist2, n0, n1, par1.1, par1.2, par2.1, par2.2, rho){
  
  # real threshold (theoretical)
  th <- thresTH2(dist1, dist2, par1.1, par1.2, par2.1, par2.2, rho)$thres
  # AUC being generated
  md <- par1.1-par2.1
  sd <- sqrt(par1.2^2+par2.2^2)
  AUC <- pnorm(0, md, sd)
  
  
  out <- 1:k %>% map_dfr(function(i) sim_sample(n0, n1, par1.1, par1.2, par2.1, par2.2, rho))
  out <- out %>% mutate(rvar=ifelse(sd.k1.>sd.k2., (sd.k1./sd.k2.)^2, (sd.k2./sd.k1.)^2)) # greater variance in the numerator
  
  out <- out %>% mutate(th11=ifelse(rvar<1.1, th.eq, th.un),
                        th13=ifelse(rvar<1.3, th.eq, th.un),
                        sd11=ifelse(rvar<1.1, th.se.eq, th.se.un),
                        sd13=ifelse(rvar<1.3, th.se.eq, th.se.un),
                        thp=ifelse(varp>0.05, th.eq, th.un),
                        sdp=ifelse(varp>0.05, th.se.eq, th.se.un),
                        th105=ifelse(rvar<1.05, th.eq, th.un),
                        sd105=ifelse(rvar<1.05, th.se.eq, th.se.un),
                        ll11=ifelse(rvar<1.1, th.ll.eq, th.ll.un),
                        ul11=ifelse(rvar<1.1, th.ul.eq, th.ul.un),
                        ll13=ifelse(rvar<1.3, th.ll.eq, th.ll.un),
                        ul13=ifelse(rvar<1.3, th.ul.eq, th.ul.un),
                        llp=ifelse(varp>0.05, th.ll.eq, th.ll.un),
                        ulp=ifelse(varp>0.05, th.ul.eq, th.ul.un),
                        ll105=ifelse(rvar<1.05, th.ll.eq, th.ll.un),
                        ul105=ifelse(rvar<1.05, th.ul.eq, th.ul.un))
  
  
  # % equal estimator
  pvar <- sum(out$varp>0.05)/k
  p105 <- sum(out$rvar<1.05)/k
  p11 <- sum(out$rvar<1.1)/k
  p13 <- sum(out$rvar<1.3)/k
  
  # bias
  b_eq <- mean(out$th.eq)-th
  b_un <- mean(out$th.un)-th
  b_p <- mean(out$thp)-th
  b_105 <- mean(out$th105)-th
  b_11 <- mean(out$th11)-th
  b_13 <- mean(out$th13)-th
  
  # MSE
  mse_eq <- mean((out$th.eq-th)^2)
  mse_un <- mean((out$th.un-th)^2)
  mse_p <- mean((out$thp-th)^2)
  mse_105 <- mean((out$th105-th)^2)
  mse_11 <- mean((out$th11-th)^2)
  mse_13 <- mean((out$th13-th)^2)
  
  # EC
  ic_eq <- matrix(c(out$th.ll.eq, out$th.ul.eq), ncol=2)
  ic_un <- matrix(c(out$th.ll.un, out$th.ul.un), ncol=2)
  ic_p <- matrix(c(out$llp, out$ulp), ncol=2)
  ic_105 <- matrix(c(out$ll105, out$ul105), ncol=2)
  ic_11 <- matrix(c(out$ll11, out$ul11), ncol=2)
  ic_13 <- matrix(c(out$ll13, out$ul13), ncol=2)
  
  cov_eq <- 1:nrow(ic_eq) %>% map_dfr(function(i) data.frame(x=within_int(ic_eq[i, 1], ic_eq[i, 2], th))) %>% summarise(x=mean(x))
  cov_un <- 1:nrow(ic_un) %>% map_dfr(function(i) data.frame(x=within_int(ic_un[i, 1], ic_un[i, 2], th))) %>% summarise(x=mean(x))
  cov_p <- 1:nrow(ic_p) %>% map_dfr(function(i) data.frame(x=within_int(ic_p[i, 1], ic_p[i, 2], th))) %>% summarise(x=mean(x))
  cov_105 <- 1:nrow(ic_105) %>% map_dfr(function(i) data.frame(x=within_int(ic_105[i, 1], ic_105[i, 2], th))) %>% summarise(x=mean(x))
  cov_11 <- 1:nrow(ic_11) %>% map_dfr(function(i) data.frame(x=within_int(ic_11[i, 1], ic_11[i, 2], th))) %>% summarise(x=mean(x))
  cov_13 <- 1:nrow(ic_13) %>% map_dfr(function(i) data.frame(x=within_int(ic_13[i, 1], ic_13[i, 2], th))) %>% summarise(x=mean(x))

  # output
  res <- data.frame(m0=par1.1, m1=par2.1, s0=par1.2, s1=par2.2, AUC, n0=n0, n1=n1, rho=rho,
                    pvar=pvar, p105=p105, p11=p11, p13=p13,
                    b_eq=b_eq, b_un=b_un, b_p=b_p, b_105=b_105, b_11=b_11, b_13=b_13,
                    mse_eq=mse_eq, mse_un=mse_un, mse_p=mse_p, mse_105=mse_105, mse_11=mse_11, mse_13=mse_13,
                    cov_eq=cov_eq$x*100, cov_un=cov_un$x*100, cov_p=cov_p$x*100, cov_105=cov_105$x*100, cov_11=cov_11$x*100, cov_13=cov_13$x*100)
  return(res)
}



#### scenarios with ratio of variances equal to 1, 1.1, 1.3
set.seed(200)

## p=0.3

# high AUC 
# m0=0,s0=1,m0=2,s0=1,rho=0.3
aux1 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 2, 1, 0.3)
aux2 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 2, 1, 0.3)
aux3 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 2, 1, 0.3)

# m0=0,s0=1,m0=2,s0=1*sqrt(1.1),rho=0.3
aux4 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 2, 1*sqrt(1.1), 0.3)
aux5 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 2, 1*sqrt(1.1), 0.3)
aux6 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 2, 1*sqrt(1.1), 0.3)

# m0=0,s0=1,m0=2,s0=1*sqrt(1.3),rho=0.3
aux7 <- sim_var(10000,"norm","norm",20,20,0,1,2,1*sqrt(1.3),0.3)
aux8 <- sim_var(10000,"norm","norm",100,100,0,1,2,1*sqrt(1.3),0.3)
aux9 <- sim_var(10000,"norm","norm",500,500,0,1,2,1*sqrt(1.3),0.3)


# moderate AUC
# m0=0,s0=1,m0=1,s0=1,rho=0.3
aux10 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 1, 1, 0.3)
aux11 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 1, 1, 0.3)
aux12 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 1, 1, 0.3)

# m0=0,s0=1,m0=1,s0=1*sqrt(1.1),rho=0.3
aux13 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 1, 1*sqrt(1.1), 0.3)
aux14 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 1, 1*sqrt(1.1), 0.3)
aux15 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 1, 1*sqrt(1.1), 0.3)

# m0=0,s0=1,m0=1,s0=1*sqrt(1.3),rho=0.3
aux16 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 1, 1*sqrt(1.3), 0.3)
aux17 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 1, 1*sqrt(1.3), 0.3)
aux18 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 1, 1*sqrt(1.3), 0.3)


# low AUC
# m0=0,s0=1,m0=0.5,s0=1,rho=0.3
aux19 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 0.5, 1, 0.3)
aux20 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 0.5, 1, 0.3)
aux21 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 0.5, 1, 0.3)

# m0=0,s0=1,m0=0.5,s0=1*sqrt(1.1),rho=0.3
aux22 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 0.5, 1*sqrt(1.1), 0.3)
aux23 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 0.5, 1*sqrt(1.1), 0.3)
aux24 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 0.5, 1*sqrt(1.1), 0.3)

# m0=0,s0=1,m0=0.5,s0=1*sqrt(1.3),rho=0.3
aux25 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 0.5, 1*sqrt(1.3), 0.3)
aux26 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 0.5, 1*sqrt(1.3), 0.3)
aux27 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 0.5, 1*sqrt(1.3), 0.3)


res1 <- rbind(aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8, aux9, aux10,
              aux11, aux12, aux13, aux14, aux15, aux16, aux17, aux18, aux19, aux20,
              aux21, aux22, aux23, aux24, aux25, aux26, aux27)


## p=0.3

# high AUC
# m0=0,s0=1,m0=2,s0=1,rho=0.1
aux28 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 2, 1, 0.1)
aux29 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 2, 1, 0.1)
aux30 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 2, 1, 0.1)

# m0=0,s0=1,m0=2,s0=1*sqrt(1.1),rho=0.1
aux31 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 2, 1*sqrt(1.1), 0.1)
aux32 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 2, 1*sqrt(1.1), 0.1)
aux33 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 2, 1*sqrt(1.1), 0.1)

# m0=0,s0=1,m0=2,s0=1*sqrt(1.3),rho=0.1
aux34 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 2, 1*sqrt(1.3), 0.1)
aux35 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 2, 1*sqrt(1.3), 0.1)
aux36 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 2, 1*sqrt(1.3), 0.1)



# moderate AUC
# m0=0,s0=1,m0=1,s0=1,rho=0.1
aux37 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 1, 1, 0.1)
aux38 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 1, 1, 0.1)
aux39 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 1, 1, 0.1)

# m0=0,s0=1,m0=1,s0=1*sqrt(1.1),rho=0.1
aux40 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 1, 1*sqrt(1.1), 0.1)
aux41 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 1, 1*sqrt(1.1), 0.1)
aux42 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 1, 1*sqrt(1.1), 0.1)

# m0=0,s0=1,m0=1,s0=1*sqrt(1.3),rho=0.1
aux43 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 1, 1*sqrt(1.3), 0.1)
aux44 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 1, 1*sqrt(1.3), 0.1)
aux45 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 1, 1*sqrt(1.3), 0.1)

# low AUC
# m0=0,s0=1,m0=0.5,s0=1,rho=0.1
aux46 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 0.5, 1, 0.1)
aux47 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 0.5, 1, 0.1)
aux48 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 0.5, 1, 0.1)

# m0=0,s0=1,m0=0.5,s0=1*sqrt(1.1),rho=0.1
aux49 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 0.5, 1*sqrt(1.1), 0.1)
aux50 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 0.5, 1*sqrt(1.1), 0.1)
aux51 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 0.5, 1*sqrt(1.1), 0.1)

# m0=0,s0=1,m0=0.5,s0=1*sqrt(1.3),rho=0.1
aux52 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 0.5,1*sqrt(1.3), 0.1)
aux53 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 0.5,1*sqrt(1.3), 0.1)
aux54 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 0.5, 1*sqrt(1.3), 0.1)



res2 <- rbind(aux27, aux28, aux29, aux30,
              aux31, aux32, aux33, aux34, aux35, aux36, aux37, aux38, aux39, aux40,
              aux41, aux42, aux43, aux44, aux45, aux46, aux47, aux48, aux49, aux50,
              aux51, aux52, aux53, aux54)



#### scenarios with ratio of variances equal to 1.05
set.seed(200)

## p=0.1
aux55 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 2, 1*sqrt(1.05), 0.1)
aux56 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 2, 1*sqrt(1.05), 0.1)
aux57 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 2, 1*sqrt(1.05), 0.1)
aux58 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 1, 1*sqrt(1.05), 0.1)
aux59 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 1, 1*sqrt(1.05), 0.1)
aux60 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 1, 1*sqrt(1.05), 0.1)
aux61 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 0.5, 1*sqrt(1.05), 0.1)
aux62 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 0.5, 1*sqrt(1.05), 0.1)
aux63 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 0.5, 1*sqrt(1.05), 0.1)



## p=0.3
aux64 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 2, 1*sqrt(1.05), 0.3)
aux65 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 2, 1*sqrt(1.05), 0.3)
aux66 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 2, 1*sqrt(1.05), 0.3)
aux67 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 1, 1*sqrt(1.05), 0.3)
aux68 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 1, 1*sqrt(1.05), 0.3)
aux69 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 1, 1*sqrt(1.05), 0.3)
aux70 <- sim_var(10000, "norm", "norm", 20, 20, 0, 1, 0.5, 1*sqrt(1.05), 0.3)
aux71 <- sim_var(10000, "norm", "norm", 100, 100, 0, 1, 0.5, 1*sqrt(1.05), 0.3)
aux72 <- sim_var(10000, "norm", "norm", 500, 500, 0, 1, 0.5, 1*sqrt(1.05), 0.3)



res3 <- rbind(aux55, aux56, aux57, aux58, aux59, aux60,
              aux61, aux62, aux63, aux64, aux65, aux66, aux67, aux68, aux69, aux70,
              aux71, aux72)



res_tot <- rbind(res1, res2, res3)

result <- res_tot[c(29:31, 56:58, 32:40, 59:61, 41:49, 62:64, 50:55, 1:3, 65:67, 4:12, 68:70, 13:21, 71:73, 22:27), ]
# 'result' contains the data presented in the tables of the article
