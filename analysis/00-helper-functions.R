# script: 00-helper-functions
# author: Charlotte Mann
# purpose: Implementation of point and variance estimators for pCRTs

#estimate k* for Des Raj difference estimator in Middleton & Aronow (2015)
est_k <- function(dat){

  pids <- unique(dat$P)

  k <- foreach(i = pids, .combine = rbind) %do% {
    mod <- lm(Y ~ . -Tr -P , data = dat %>% filter(P != i))

    coefficients(mod)[-1]
  }

  colnames(k) <- paste0("k",colnames(k),"k")
  rownames(k) <- NULL

  k <- k %>%
    data.frame() %>%
    mutate(P = pids)

  return(k)

}

# Des Raj Difference estimator from Middleton & Aronow (2015)
# Y are average cluster outcomes
# Z are total cluster covariates
des_raj_dif <- function(Y,Tr,Z=NULL,P,n,
                        type = c("adj","ht")){

  # n will be a covariate always following the MA paper for des raj
  if(is.null(Z)){
    Z = n
  }else{
    Z = cbind(n,Z)
  }

  Z <- as.matrix(Z)
  colnames(Z) <- paste0("Z", 1:(ncol(Z)))

  #want totals for this
  Y = Y*n

  ordered = pair(Y,Tr,Z,P,n)$ordered
  dat = data.frame(P,Y,Tr,Z)

  #leave-one-out estimations of k vector
  k_ests = est_k(dat)

  #make sure it all lines up
  ordered = ordered %>%
    left_join(k_ests, by = "P")

  #treatment and control covariates as matrices
  cov1 = ordered %>% select(starts_with("Z")&ends_with("1")) %>% as.matrix()
  cov0 = ordered %>% select(starts_with("Z")&ends_with("2")) %>% as.matrix()

  k <- ordered %>% select(starts_with("k")) %>% as.matrix()

  #difference estimator
  M = nrow(ordered)
  N = sum(n)

  cov_means = rep(1,M)%*%t(apply(Z,2, mean))

  #also allow to calculate the actual horvitz thompson estimator and MA variance
  if(type == "ht"){
    U1 = ordered$Y1
    U0 = ordered$Y2
  }else{
    U1 = ordered$Y1 - apply(k*(cov1 - cov_means),1, sum, na.rm = T)
    U0 = ordered$Y2 - apply(k*(cov0 - cov_means),1, sum, na.rm = T)
  }

  tauhat = 2*(sum(U1) - sum(U0))/N

  varhat = 16*M^2*var(c(U1,U0))/(N^2*(2*M-1))

  return(list(tauhat = tauhat, varhat = varhat))

}


# variance estimators for Hajek in de Chaisemartin and Ramirez-Cuella (2020)
hajek_var_cr <- function(mod.dat, ypred){

  mod.dat$yhat = ypred
  dat <- mod.dat %>%
    select(P, Y, Tr, n, yhat, pair_n) %>%
    mutate(r = Y*n - yhat*n) %>%
    pivot_wider(names_from = pair_n,
                values_from = c(Y, Tr, n, yhat, r),
                names_glue = "{.value}{pair_n}")

  Rt <- dat$Tr1*dat$r1 + dat$Tr2*dat$r2
  Rc <- (1-dat$Tr1)*dat$r1 + (1-dat$Tr2)*dat$r2
  Nt <- sum(dat$Tr1*dat$n1 + dat$Tr2*dat$n2)
  Nc <- sum((1-dat$Tr1)*dat$n1 + (1-dat$Tr2)*dat$n2)

  var_pair = sum((Rt/Nt - Rc/Nc)^2)
  var_unit = sum((Rt/Nt)^2 + (Rc/Nc)^2)

  return(list(vhat_pair = var_pair,
              vhat_unit = var_unit))

}

# variance estimator for fe estimate from Schochet et al (2021)
spmk_var_est <- function(mod.dat, ypred, v=0){

  mod.dat$yhat = ypred
  dat <- mod.dat %>%
    select(P, Y, Tr, n, yhat, pair_n) %>%
    mutate(r = Y - yhat) %>%
    pivot_wider(names_from = pair_n,
                values_from = c(Y, Tr, n, yhat, r),
                names_glue = "{.value}{pair_n}")

  w_i = dat$n1*dat$n2/(dat$n1+dat$n2)
  M = nrow(dat)

  vhat = 2*M/((M-v-1)*sum(w_i)^2)*sum(w_i^2*(dat$r1^2 + dat$r2^2))

  return(vhat)
}

# effect point estimator from Imai, King, Nall (2009), returning variance estimators as well from
# Imai, King, Nall (2009) and de Chaisemartin and Ramirez-Cuella (2020)
weighted_effect_est <- function(ordered, n_ordered,
                                weight = c("arithmetic","harmonic")){

  tau_i = ordered$Y1 - ordered$Y2

  if(weight == "arithmetic")  w_i = n_ordered$n1 + n_ordered$n2
  if(weight == "harmonic")    w_i = n_ordered$n1*n_ordered$n2/(n_ordered$n1 + n_ordered$n2)

  tauhat = (1/sum(w_i))*sum(w_i*tau_i)

  M = length(w_i)
  W = sum(w_i)

  varhat_ikn = M/(M-1)*sum((w_i*tau_i/W - (1/M)*tauhat)^2)
  varhat_cr = NULL

  if(weight == "harmonic") varhat_cr =  sum((w_i/W)^2*(tau_i - tauhat)^2)

  return(list(tau = tauhat,
              vhat_ikn = varhat_ikn,
              vhat_cr = varhat_cr))

}
