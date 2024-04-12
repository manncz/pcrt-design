# script: 00-sim-helpers
# author: Charlotte Mann
# purpose: Helper functions to run simulations under different settings

#####################      FUNCTION FOR GENERATING CLUSTER SIZES    #################################

gen_nki <- function(M, a, b, c, d){
  
  if(a == b){
    ni <- rep(a, M, replace = T)
  }else{
    ni <- sample(seq(a, b), M, replace = T)
  }

  m1i <- sample(seq(c, d), M, replace = T)
  m2i <- sample(seq(c, d), M, replace = T)

  n <- cbind(ni + m1i, ni + m2i)

  return(n)

}

#################  FUNCTION FOR GENERATING CLUSTER-LEVEL TREATMENT EFFECTS  #########################


gen_tauki <- function(base_tau = 1,
                      tau_var = 0,
                      M, clust_size, cov,
                      fslope, fintercept, gslope=0){

  f = apply(clust_size, c(1,2), FUN = function(x){fslope*(x-fintercept)})
  g = gslope*cov
  tau_i1 = rnorm(M, mean = base_tau, sd = sqrt(tau_var))
  tau_i = cbind(tau_i1, tau_i1) + f + g

  return(tau_i)

}


###################################### GENERATE SIMULATION DATA  #########################################

gen_clust_dat <- function(M = 20,                                                        # number of pairs
                          nimin = 15, nimax = 35, mkimin = 0, mkimax = 0,                # parameters for generating n
                          cov_clustv = 0,                                                # cluster-level variance explained by covariate
                          alpha = 5, pairv = 1, clusterv = 1,                            # parameters for generating intercept and cluster / pair effects
                          indv = 1,                                                      # variance of individual error term
                          nbeta = 0,                                                    # parameter for effect of cluster size on control outcome
                          base_tau = 1, tau_var = 0, tau_fslope = 0, tau_fintercept = 0, # parameters for generating tau
                          tau_gslope = 0,
                          loc_shift = 0,                                                 # ability to add a location shift
                          seed = NULL
                          ){
  # set seed
  if(!is.null(seed)){set.seed(seed)}
  print(paste("The seed is set to:", seed))

  #cluster sizes
  n <- gen_nki(M, a = nimin, b = nimax, c = mkimin, d = mkimax)

  #errors with variance indv
  err <- apply(n, c(1,2), FUN = function(x){sqrt(indv)*rnorm(x,0,sd = 1)}, simplify = F)

  #cluster-level element of covariate
  cov <- sqrt(cov_clustv)*matrix(rnorm(M*2), ncol = 2)

  #cluster-level error
  gamma_ki = sqrt(clusterv)*matrix(rnorm(M*2), ncol = 2)

  #pair-level error
  alpha_i <- rnorm(M, 0, sd = sqrt(pairv))
  alpha_i <- cbind(alpha_i, alpha_i)

  # add together errors and covariate contribution with loop function to maintain structure
  # calculate control potential outcomes, either with a covariate or not
  po_c = err

  for(i in 1:M){
    for(k in 1:2){
      po_c[i,k][[1]] = alpha + alpha_i[i,k] + cov[i,k] + gamma_ki[i,k] + nbeta*n[i,k] + err[i,k][[1]] + loc_shift
    }
  }

  #subsume the pair effect into the covariate
  cov = cov + alpha_i
  colnames(cov) <- NULL

  # treatment effect - allow for different cluster-level tau, but not individual
  tau_i <- gen_tauki(base_tau, tau_var = tau_var, clust_size = n, M = M,
                     fslope = tau_fslope, fintercept = tau_fintercept,
                     cov = cov, gslope = tau_gslope)

  colnames(tau_i) <- c("tau1","tau2")
  rownames(tau_i) <- paste0("P", 1:M)

  # create individual data
  ind.dat <- foreach(i = 1:M, .combine = rbind) %do% {

    dat1 <- data.frame(P = paste0("P", i),
                       pair_id = 1,
                       n = n[i,1],
                       Z1 = cov[i,1],
                       po_c = po_c[i,1][[1]],
                       tau = tau_i[i,1])

    dat2 <- data.frame(P = paste0("P", i),
                       pair_id = 2,
                       n = n[i,2],
                       Z1 = cov[i,2],
                       po_c = po_c[i,2][[1]],
                       tau = tau_i[i,2])

   check <- dat1 %>%
      rbind(dat2) %>%
      mutate(po_t = po_c + tau)

  }

  # calculate cluster means and treatment potential outcomes
  po <- ind.dat %>%
    group_by(P, pair_id) %>%
    summarize(po_t = mean(po_t), po_c = mean(po_c)) %>%
    pivot_wider(values_from = po_t:po_c,
                names_from = pair_id,
                names_sep = "") %>%
    ungroup()

  po_mc <- po %>%
    select(starts_with("po_c")) %>%
    as.matrix()

  po_mt <- po %>%
    select(starts_with("po_t")) %>%
    as.matrix()

  rownames(po_mc) <- rownames(po_mt) <- po$P

  # generate treatment vector
  T_i <- sample(c(0,1), M, replace = T)

  # combine into wide dataset
  dat.orig <- ind.dat %>%
    group_by(P, pair_id) %>%
    summarize_all(mean) %>%
    pivot_wider(values_from = n:po_t,
                names_from = pair_id,
                names_sep = "") %>%
    ungroup() %>%
    mutate(Tr1 = T_i, Tr2 = (1-T_i),
           Y1 = Tr1*po_t1 + (1-Tr1)*po_c1,
           Y2 = Tr2*po_t2 + (1-Tr2)*po_c2) %>%
    select(starts_with("Tr"),starts_with("Y"), starts_with("Z"), starts_with("n"), P)

  #ensure tau is in the right order, after other code changes
  tau = tau_i[dat.orig$P,]

  return(list(tau = tau,
              ind_dat = ind.dat,
              po_c = po_mc,
              po_t = po_mt,
              dat = dat.orig))

}


###################################### RUN SIMULATIONS #########################################


################################# SIM FOR ONE DATA GEN  ########################################


######################   WITH ALL ESTIMATORS  ###########################

all_est_sims <- function(n.sim = 1000, seed = NULL,
                           dat.orig, po_mt, po_mc, tau){

  # set seed
  if(!is.null(seed)){set.seed(seed)}
  print(paste("The seed is set to:", seed))

  M = nrow(dat.orig)
  N = sum(dat.orig$n1, dat.orig$n2)

  z <- qnorm(p = 1-(.05/2))
  treg <- qt(p = 1-(.05/2), df = 2*M-2)
  tfe <- qt(p = 1-(.05/2), df = M-1)

  tregcov <- qt(p = 1-(.05/2), df = 2*M-3)
  tregcovn <- qt(p = 1-(.05/2), df = 2*M-4)
  tfecov <- qt(p = 1-(.05/2), df = M-2)
  tfecovn <- qt(p = 1-(.05/2), df = M-3)

  sim.out <- foreach(i = 1:n.sim, .combine=rbind) %do% {

    it.ob <- c()

    Tr <- sample(c(0,1), M, replace = T)
    # create data with new treatment effect vector and observed outcomes
    dat <- dat.orig %>%
      mutate(Tr1 = Tr, Tr2 = 1-Tr,
             Y1 = Tr1*(po_mt[,1]) + (1-Tr1)*po_mc[,1],
             Y2 = (1-Tr1)*(po_mt[,2]) + (Tr1)*po_mc[,2]
      )

    #reshape longer as input for functions
    dat.long <- dat %>%
      pivot_longer(cols = Tr1:n2,
                   names_pattern = "(.*)(\\d)",
                   names_to = c(".value", "pair_n"))


    #===========================================  Hajek   ===========================================#
    # equivalent to regression estimator with no covariates
    dmlm = lm(Y ~ Tr, data = dat.long, weights = n)

    it.ob["dm_tau"] = coefficients(dmlm)["Tr"]
    it.ob["dm_var.reg"] = vcov(dmlm)["Tr","Tr"]
    it.ob["dm_covg.reg"] = as.numeric(tau >= it.ob["dm_tau"]-(treg*sqrt(it.ob["dm_var.reg"])) &
                                        tau <= it.ob["dm_tau"]+(treg*sqrt(it.ob["dm_var.reg"])))

    #robust SEs
    #this matches what you would get with STATA "robust"
    it.ob["dm_var.robust"] = vcovHC(dmlm, type = "HC1")["Tr","Tr"]
    it.ob["dm_covg.robust"] = as.numeric(tau >= it.ob["dm_tau"]-(treg*sqrt(it.ob["dm_var.robust"])) &
                                           tau <= it.ob["dm_tau"]+(treg*sqrt(it.ob["dm_var.robust"])))

    #c&r SEs "paired"
    yhat = predict(dmlm)
    it.ob["dm_var.cr"] = hajek_var_cr(mod.dat = dat.long, yhat)$vhat_pair
    it.ob["dm_covg.cr"] = as.numeric(tau >= it.ob["dm_tau"]-(treg*sqrt(it.ob["dm_var.cr"])) &
                                       tau <= it.ob["dm_tau"]+(treg*sqrt(it.ob["dm_var.cr"])))

    #===================================     WLS with Covariates    =================================#
    #-------------- WLS + N   -----------------#

    reglm = lm(Y ~ Tr + n, data = dat.long, weights = n)

    it.ob["wlsn_tau"] = coefficients(reglm)["Tr"]
    it.ob["wlsn_var.robust"] = vcovHC(reglm, type = "HC1")["Tr","Tr"]
    it.ob["wlsn_covg.robust"] = as.numeric(tau >= it.ob["wlsn_tau"]-(tregcov*sqrt(it.ob["wlsn_var.robust"])) &
                                            tau <= it.ob["wlsn_tau"]+(tregcov*sqrt(it.ob["wlsn_var.robust"])))

    it.ob["wlsn_r2"] = summary(reglm)$adj.r.squared

    #-------------- WLS + Cov -----------------#
    reglm = lm(Y ~ Tr + Z1, data = dat.long, weights = n)

    it.ob["wlsx_tau"] = coefficients(reglm)["Tr"]
    it.ob["wlsx_var.robust"] = vcovHC(reglm, type = "HC1")["Tr","Tr"]
    it.ob["wlsx_covg.robust"] = as.numeric(tau >= it.ob["wlsx_tau"]-(tregcov*sqrt(it.ob["wlsx_var.robust"])) &
                                            tau <= it.ob["wlsx_tau"]+(tregcov*sqrt(it.ob["wlsx_var.robust"])))

    it.ob["wlsx_r2"] = summary(reglm)$adj.r.squared

    #-------------- WLS + Cov & N -----------------#

    reglm = lm(Y ~ Tr + Z1 + n, data = dat.long, weights = n)

    it.ob["wlsxn_tau"] = coefficients(reglm)["Tr"]
    it.ob["wlsxn_var.robust"] = vcovHC(reglm, type = "HC1")["Tr","Tr"]
    it.ob["wlsxn_covg.robust"] = as.numeric(tau >= it.ob["wlsxn_tau"]-(tregcovn*sqrt(it.ob["wlsxn_var.robust"])) &
                                            tau <= it.ob["wlsxn_tau"]+(tregcovn*sqrt(it.ob["wlsxn_var.robust"])))

    it.ob["wlsxn_r2"] = summary(reglm)$adj.r.squared


    #=====================  Imai, King, Nall estimator (arithmetic mean weights) ====================#

    #Imai weighted estimate - this is also the same as WLS with cluster weight as the total inds in the pair
    #function pair is from the dRCT package
    paired <- pair(Y=dat.long$Y, Tr=dat.long$Tr, Z=dat.long$n, P = dat.long$P, n =dat.long$n)

    #imai weighted estimate with arithmethic mean weights
    waest = weighted_effect_est(ordered = paired$ordered, n_ordered = paired$n_ordered,
                                weight = "arithmetic")

    it.ob["amw_tau"] = waest$tau

    #ikn variance
    it.ob["amw_var.ikn"] = waest$vhat_ikn
    it.ob["amw_covg.ikn"] = as.numeric(tau >= it.ob["amw_tau"]-(tfe*sqrt(it.ob["amw_var.ikn"])) &
                                        tau <= it.ob["amw_tau"]+(tfe*sqrt(it.ob["amw_var.ikn"])))

    #WLS equivalent
    dat.ws <- dat.long %>%
      group_by(P) %>%
      mutate(nP = sum(n))

    wslm = lm(Y ~ Tr, data = dat.ws, weights = nP)

    #regression variance
    it.ob["amw_var.reg"] = vcov(wslm)["Tr","Tr"]
    it.ob["amw_covg.reg"] = as.numeric(tau >= it.ob["amw_tau"]-(treg*sqrt(it.ob["amw_var.reg"])) &
                                        tau <= it.ob["amw_tau"]+(treg*sqrt(it.ob["amw_var.reg"])))

    #robust variance
    it.ob["amw_var.robust"] = vcovHC(wslm, type = "HC1")["Tr","Tr"]
    it.ob["amw_covg.robust"] = as.numeric(tau >= it.ob["amw_tau"]-(treg*sqrt(it.ob["amw_var.robust"])) &
                                           tau <= it.ob["amw_tau"]+(treg*sqrt(it.ob["amw_var.robust"])))

    #========================  Fixed Effect Estimator (harmonic mean weights) =======================#
    #regression estimator with fixed effects (Chaisemartin & Ramirez-Cuellar working paper)

    #imai weighted estimate but with harmonic mean weights
    feest = weighted_effect_est(ordered = paired$ordered, n_ordered = paired$n_ordered,
                                weight = "harmonic")

    it.ob["fe_tau"] = feest$tau

    #ikn variance estimator
    it.ob["fe_var.ikn"] = feest$vhat_ikn
    it.ob["fe_covg.ikn"] = as.numeric(tau >= it.ob["fe_tau"]-(tfe*sqrt(it.ob["fe_var.ikn"])) &
                                        tau <= it.ob["fe_tau"]+(tfe*sqrt(it.ob["fe_var.ikn"])))

    #c&r variance estimator
    it.ob["fe_var.cr"] = feest$vhat_cr
    it.ob["fe_covg.cr"] = as.numeric(tau >= it.ob["fe_tau"]-(tfe*sqrt(it.ob["fe_var.cr"])) &
                                       tau <= it.ob["fe_tau"]+(tfe*sqrt(it.ob["fe_var.cr"])))

    #the fe_tau above is equivalent to the coefficient from this estimator
    felm = lm(Y ~ Tr + P, data = dat.long, weights = n)

    #spmk variance estimator
    yhat = predict(felm)

    it.ob["fe_var.spmk"] = spmk_var_est(dat.long, ypred=yhat, v = 0)
    it.ob["fe_covg.spmk"] = as.numeric(tau >= it.ob["fe_tau"]-(tfe*sqrt(it.ob["fe_var.spmk"])) &
                                         tau <= it.ob["fe_tau"]+(tfe*sqrt(it.ob["fe_var.spmk"])))

    #standard wlm variance estimator
    it.ob["fe_var.reg"] = vcov(felm)["Tr","Tr"]
    it.ob["fe_covg.reg"] = as.numeric(tau >= it.ob["fe_tau"]-(tfe*sqrt(it.ob["fe_var.reg"])) &
                                        tau <= it.ob["fe_tau"]+(tfe*sqrt(it.ob["fe_var.reg"])))

    #robust wlm variance estimator
    it.ob["fe_var.robust"] = vcovHC(felm, type = "HC1")["Tr","Tr"]
    it.ob["fe_covg.robust"] = as.numeric(tau >= it.ob["fe_tau"]-(tfe*sqrt(it.ob["fe_var.robust"])) &
                                           tau <= it.ob["fe_tau"]+(tfe*sqrt(it.ob["fe_var.robust"])))


    #======================================  Fixed Effect + COV =====================================#
    #-------------- FE + N -----------------#
    felmcov = lm(Y ~ Tr + n + P, data = dat.long, weights = n)
    it.ob["fen_tau"] =coefficients(felmcov)["Tr"]

    #standard var
    it.ob["fen_var.reg"] = vcov(felmcov)["Tr","Tr"]
    it.ob["fen_covg.reg"] = as.numeric(tau >= it.ob["fen_tau"]-(tfecov*sqrt(it.ob["fen_var.reg"])) &
                                           tau <= it.ob["fen_tau"]+(tfecov*sqrt(it.ob["fen_var.reg"])))

    #robust var
    it.ob["fen_var.robust"] = vcovHC(felmcov, type = "HC1")["Tr","Tr"]
    it.ob["fen_covg.robust"] = as.numeric(tau >= it.ob["fen_tau"]-(tfecov*sqrt(it.ob["fen_var.robust"])) &
                                              tau <= it.ob["fen_tau"]+(tfecov*sqrt(it.ob["fen_var.robust"])))

    it.ob["fen_r2"] = summary(felmcov)$adj.r.squared


    #spmk var est
    yhat = predict(felmcov)

    it.ob["fen_var.spmk"] = spmk_var_est(dat.long, ypred=yhat, v = 0)
    it.ob["fen_covg.spmk"] = as.numeric(tau >= it.ob["fen_tau"]-(tfe*sqrt(it.ob["fen_var.spmk"])) &
                                         tau <= it.ob["fen_tau"]+(tfe*sqrt(it.ob["fen_var.spmk"])))


    #-------------- FE + Cov -----------------#
    felmcov = lm(Y ~ Tr + Z1 + P, data = dat.long, weights = n)
    it.ob["feadj_tau"] =coefficients(felmcov)["Tr"]

    #standard var
    it.ob["feadj_var.reg"] = vcov(felmcov)["Tr","Tr"]
    it.ob["feadj_covg.reg"] = as.numeric(tau >= it.ob["feadj_tau"]-(tfecov*sqrt(it.ob["feadj_var.reg"])) &
                                              tau <= it.ob["feadj_tau"]+(tfecov*sqrt(it.ob["feadj_var.reg"])))

    #robust var
    it.ob["feadj_var.robust"] = vcovHC(felmcov, type = "HC1")["Tr","Tr"]
    it.ob["feadj_covg.robust"] = as.numeric(tau >= it.ob["feadj_tau"]-(tfecov*sqrt(it.ob["feadj_var.robust"])) &
                                              tau <= it.ob["feadj_tau"]+(tfecov*sqrt(it.ob["feadj_var.robust"])))

    it.ob["feadj_r2"] = summary(felmcov)$adj.r.squared


    #-------------- FE + Cov & N -----------------#
    felmcov = lm(Y ~ Tr + Z1 + P + n, data = dat.long, weights = n)
    it.ob["feadjn_tau"] =coefficients(felmcov)["Tr"]

    #standard var
    it.ob["feadjn_var.reg"] = vcov(felmcov)["Tr","Tr"]
    it.ob["feadjn_covg.reg"] = as.numeric(tau >= it.ob["feadjn_tau"]-(tfecov*sqrt(it.ob["feadjn_var.reg"])) &
                                           tau <= it.ob["feadjn_tau"]+(tfecov*sqrt(it.ob["feadjn_var.reg"])))

    #robust var
    it.ob["feadjn_var.robust"] = vcovHC(felmcov, type = "HC1")["Tr","Tr"]
    it.ob["feadjn_covg.robust"] = as.numeric(tau >= it.ob["feadjn_tau"]-(tfecov*sqrt(it.ob["feadjn_var.robust"])) &
                                              tau <= it.ob["feadjn_tau"]+(tfecov*sqrt(it.ob["feadjn_var.robust"])))

    yhat = predict(felmcov)

    it.ob["feadjn_var.spmk"] = spmk_var_est(dat.long, ypred=yhat, v = 0)
    it.ob["feadjn_covg.spmk"] = as.numeric(tau >= it.ob["feadjn_tau"]-(tfe*sqrt(it.ob["feadjn_var.spmk"])) &
                                          tau <= it.ob["feadjn_tau"]+(tfe*sqrt(it.ob["feadjn_var.spmk"])))

    #======================================= Horvitz-Thompson =======================================#

    #Horvitz Thompson w Middleton & Aronow variance est
    htout <- des_raj_dif(Y = dat.long$Y, Tr = dat.long$Tr, P = dat.long$P, n = dat.long$n, type = "ht")

    it.ob["ht_tau"] = htout$tauhat


    it.ob["ht_var.ma"] = htout$varhat
    it.ob["ht_covg.ma"] = as.numeric(tau >= it.ob["ht_tau"]-(z*sqrt(it.ob["ht_var.ma"])) &
                                       tau <= it.ob["ht_tau"]+(z*sqrt(it.ob["ht_var.ma"])))


    # loop variance estimator
    it.ob["ht_var.loop"] =  p_loop_var(assigned = paired$agg, v1 = rep(0,M), v2 = rep(0,M), n_assigned = paired$n_assigned)
    it.ob["ht_covg.loop"] = as.numeric(tau >= it.ob["ht_tau"]-(z*sqrt(it.ob["ht_var.loop"])) &
                                       tau <= it.ob["ht_tau"]+(z*sqrt(it.ob["ht_var.loop"])))

    #Su & Ding - OLS equivalent to HT estimator (total outcome and total covariates)
    dat.long <- dat.long %>%
      ungroup() %>%
      mutate(Ytilde = n*Y*(2*M/N),
             ntilde = n - mean(n),
             xtilde = n*(Z1-sum(Z1*n)/N))

    htlm <- lm(Ytilde ~ Tr, data = dat.long)

    #general OLS
    it.ob["ht_var.reg"] = vcov(htlm)["Tr","Tr"]
    it.ob["ht_covg.reg"] = as.numeric(tau >= it.ob["ht_tau"]-(treg*sqrt(it.ob["ht_var.reg"])) &
                                        tau <= it.ob["ht_tau"]+(treg*sqrt(it.ob["ht_var.reg"])))
    #robust OLS
    it.ob["ht_var.robust"] = vcovHC(htlm, type = "HC1")["Tr","Tr"]
    it.ob["ht_covg.robust"] = as.numeric(tau >= it.ob["ht_tau"]-(treg*sqrt(it.ob["ht_var.robust"])) &
                                           tau <= it.ob["ht_tau"]+(treg*sqrt(it.ob["ht_var.robust"])))

    #==================================== Des Raj Difference    ====================================#

    #Middleton & Aronow - des raj difference with only n
    drdout <- des_raj_dif(Y = dat.long$Y, Tr = dat.long$Tr, P = dat.long$P, n = dat.long$n, type = "adj")
    it.ob["drd_tau"] = drdout$tauhat
    it.ob["drd_var.ma"] = drdout$varhat
    it.ob["drd_covg.ma"] = as.numeric(tau >= drdout$tauhat-(z*sqrt(drdout$varhat)) &
                                        tau <= drdout$tauhat+(z*sqrt(drdout$varhat)))

    #Middleton & Aronow - des raj difference with n and cov
    drdout <- des_raj_dif(Y = dat.long$Y, Tr = dat.long$Tr, P = dat.long$P, n = dat.long$n, Z = dat.long$Z1*dat.long$n, type = "adj")
    it.ob["drdcov_tau"] = drdout$tauhat
    it.ob["drdcov_var.ma"] = drdout$varhat
    it.ob["drdcov_covg.ma"] = as.numeric(tau >= drdout$tauhat-(z*sqrt(drdout$varhat)) &
                                        tau <= drdout$tauhat+(z*sqrt(drdout$varhat)))

    #=========================================== Su & Ding ==========================================#
    tsu1 <- qt(p = 1-(.05/2), df = 2*M-4)

    # recommended Su&Ding estimator for cluster randomized trials
    htnlm <- lm(Ytilde ~ Tr*ntilde, data = dat.long)

    it.ob["htadjn_tau"] = coefficients(htnlm)["Tr"]

    #robust OLS as suggested by Su&Ding
    it.ob["htadjn_var.robust"] = vcovHC(htnlm, type = "HC1")["Tr","Tr"]
    it.ob["htadjn_covg.robust"] = as.numeric(tau >= it.ob["htadjn_tau"]-(tsu1*sqrt(it.ob["htadjn_var.robust"])) &
                                              tau <= it.ob["htadjn_tau"]+(tsu1*sqrt(it.ob["htadjn_var.robust"])))

    tsu2 <- qt(p = 1-(.05/2), df = 2*M-6)

    # with covariate
    htnlm <- lm(Ytilde ~ Tr*(ntilde + xtilde), data = dat.long)

    it.ob["htadjcov_tau"] = coefficients(htnlm)["Tr"]

    #robust OLS
    it.ob["htadjcov_var.robust"] = vcovHC(htnlm, type = "HC1")["Tr","Tr"]
    it.ob["htadjcov_covg.robust"] = as.numeric(tau >= it.ob["htadjcov_tau"]-(tsu2*sqrt(it.ob["htadjcov_var.robust"])) &
                                              tau <= it.ob["htadjcov_tau"]+(tsu2*sqrt(it.ob["htadjcov_var.robust"])))

    #======================================  Random Effects =====================================#

    relm <- lmer(Y ~ Tr + (1|P), data = dat.long, weights = n)

    it.ob["re_tau"] = fixef(relm)["Tr"]
    it.ob["re_var.reg"] = vcov(relm)["Tr","Tr"]

    #use the confidence interval calculation built into LMER (does something complicated)
    #confinttemp = suppressMessages(confint(relm))
    # that gave WAY to narrow of intervals... do it by hand instead. Not clear what the degrees of freedom should be, so
    # use the same as ols
    it.ob["re_covg.reg"] = as.numeric(tau >= it.ob["re_tau"]-(treg*sqrt(it.ob["re_var.reg"])) &
                                          tau <= it.ob["re_tau"]+(treg*sqrt(it.ob["re_var.reg"])))


    #with clustersize
    relm <- lmer(Y ~ Tr + n + (1|P), data = dat.long, weights = n)

    it.ob["ren_tau"] = fixef(relm)["Tr"]
    it.ob["ren_var.reg"] = vcov(relm)["Tr","Tr"]
    it.ob["ren_covg.reg"] = as.numeric(tau >= it.ob["ren_tau"]-(tregcov*sqrt(it.ob["ren_var.reg"])) &
                                           tau <= it.ob["ren_tau"]+(tregcov*sqrt(it.ob["ren_var.reg"])))

    #with covariate and cluster size
    relm <- lmer(Y ~ Tr + n + Z1 + (1|P), data = dat.long, weights = n)

    it.ob["readjn_tau"] = fixef(relm)["Tr"]
    it.ob["readjn_var.reg"] = vcov(relm)["Tr","Tr"]
    it.ob["readjn_covg.reg"] = as.numeric(tau >= it.ob["readjn_tau"]-(tregcov*sqrt(it.ob["readjn_var.reg"])) &
                                       tau <= it.ob["readjn_tau"]+(tregcov*sqrt(it.ob["readjn_var.reg"])))

    #============================================ P-LOOP ============================================#

    #LOO MI point and variance estimators
    lpout = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, n =dat.long$n, pred = p_simple)

    it.ob["p_mi_tau"] = lpout["tauhat"]
    it.ob["p_mi_var.loop"] = lpout["varhat"]
    it.ob["p_mi_covg.loop"] = as.numeric(tau >= lpout["tauhat"]-(z*sqrt(lpout["varhat"])) &
                                              tau <= lpout["tauhat"]+(z*sqrt(lpout["varhat"])))

      #see what coverage is like w t-dist
    it.ob["p_mi_covg.loopt"] = as.numeric(tau >= lpout["tauhat"]-(tfe*sqrt(lpout["varhat"])) &
                                              tau <= lpout["tauhat"]+(tfe*sqrt(lpout["varhat"])))

    #our point and variance estimators - using n as a covariate
    lpoutn = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, n =dat.long$n, Z = dat.long$n, pred = p_ols_interp,
                    weighted_imp = T)

    it.ob["p_wlsn_tau"] = lpoutn["tauhat"]
    it.ob["p_wlsn_var.loop"] = lpoutn["varhat"]
    it.ob["p_wlsn_covg.loop"] = as.numeric(tau >= lpoutn["tauhat"]-(z*sqrt(lpoutn["varhat"])) &
                                             tau <= lpoutn["tauhat"]+(z*sqrt(lpoutn["varhat"])))

    it.ob["p_wlsn_covg.loopt"] = as.numeric(tau >= lpoutn["tauhat"]-(tfecov*sqrt(lpoutn["varhat"])) &
                                             tau <= lpoutn["tauhat"]+(tfecov*sqrt(lpoutn["varhat"])))
    #
    # #our point and variance estimators with covariate
    # lpoutcov = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, n =dat.long$n, Z = dat.long$Z1, pred = p_ols_interp)
    #
    # it.ob["p_adj_tau"] = lpoutcov["tauhat"]
    # it.ob["p_adj_var.loop"] = lpoutcov["varhat"]
    # it.ob["p_adj_covg.loop"] = as.numeric(tau >= lpoutcov["tauhat"]-(z*sqrt(lpoutcov["varhat"])) &
    #                                          tau <= lpoutcov["tauhat"]+(z*sqrt(lpoutcov["varhat"])))
    #
    # it.ob["p_adj_covg.loopt"] = as.numeric(tau >= lpoutcov["tauhat"]-(tfecov*sqrt(lpoutcov["varhat"])) &
    #                                         tau <= lpoutcov["tauhat"]+(tfecov*sqrt(lpoutcov["varhat"])))


    #our point and variance estimators with covariate and weighted
    lpoutcovw = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, n =dat.long$n, Z = dat.long$Z1,
                       pred = p_ols_interp, weighted_imp = T)

    it.ob["p_wlsx_tau"] = lpoutcovw["tauhat"]
    it.ob["p_wlsx_var.loop"] = lpoutcovw["varhat"]
    it.ob["p_wlsx_covg.loop"] = as.numeric(tau >= lpoutcovw["tauhat"]-(z*sqrt(lpoutcovw["varhat"])) &
                                            tau <= lpoutcovw["tauhat"]+(z*sqrt(lpoutcovw["varhat"])))

    it.ob["p_wlsx_covg.loopt"] = as.numeric(tau >= lpoutcovw["tauhat"]-(tfecov*sqrt(lpoutcovw["varhat"])) &
                                             tau <= lpoutcovw["tauhat"]+(tfecov*sqrt(lpoutcovw["varhat"])))


    #our point and variance estimators with covariate and cluster size and weighted
    lpoutcovnw = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, n =dat.long$n, Z = cbind(dat.long$Z1,dat.long$n),
                       pred = p_ols_interp, weighted_imp = T)

    it.ob["p_wlsxn_tau"] = lpoutcovnw["tauhat"]
    it.ob["p_wlsxn_var.loop"] = lpoutcovnw["varhat"]
    it.ob["p_wlsxn_covg.loop"] = as.numeric(tau >= lpoutcovnw["tauhat"]-(z*sqrt(lpoutcovnw["varhat"])) &
                                             tau <= lpoutcovnw["tauhat"]+(z*sqrt(lpoutcovnw["varhat"])))

    it.ob["p_wlsxn_covg.loopt"] = as.numeric(tau >= lpoutcovnw["tauhat"]-(tfecov*sqrt(lpoutcovnw["varhat"])) &
                                              tau <= lpoutcovnw["tauhat"]+(tfecov*sqrt(lpoutcovnw["varhat"])))


    it.ob


  }

  return(sim.out)

}



######################   WITH ALL NON COVARIATE ADJUSTED ESTIMATORS  ###########################

unadj_est_sims <- function(n.sim = 1000, seed = NULL,
                           dat.orig, po_mt, po_mc, tau){

  # set seed
  if(!is.null(seed)){set.seed(seed)}
  print(paste("The seed is set to:", seed))

  M = nrow(dat.orig)
  N = sum(dat.orig$n1, dat.orig$n2)

  z <- qnorm(p = 1-(.05/2))
  treg <- qt(p = 1-(.05/2), df = 2*M-2)
  tfe <- qt(p = 1-(.05/2), df = M-1)

  trusums = data.frame(s1 = (po_mt[,1] + po_mc[,2])/2,
                       s2 = (po_mt[,2] + po_mc[,1])/2) %>%
    mutate(P = dat.orig$P) %>%
    arrange(P) %>%
    select(-P) %>%
    as.matrix()

  sim.out <- foreach(i = 1:n.sim, .combine=rbind) %do% {

    it.ob <- c()

    Tr <- sample(c(0,1), M, replace = T)
    # create data with new treatment effect vector and observed outcomes
    dat <- dat.orig %>%
      mutate(Tr1 = Tr, Tr2 = 1-Tr,
             Y1 = Tr1*(po_mt[,1]) + (1-Tr1)*po_mc[,1],
             Y2 = (1-Tr1)*(po_mt[,2]) + (Tr1)*po_mc[,2]
      )

    #reshape longer as input for functions
    dat.long <- dat %>%
      pivot_longer(cols = Tr1:n2,
                   names_pattern = "(.*)(\\d)",
                   names_to = c(".value", "pair_n"))


    #===========================================  Hajek   ===========================================#
    # equivalent to regression estimator with no covariates
    dmlm = lm(Y ~ Tr, data = dat.long, weights = n)

    it.ob["dm_tau"] = coefficients(dmlm)["Tr"]
    it.ob["dm_var.reg"] = vcov(dmlm)["Tr","Tr"]
    it.ob["dm_covg.reg"] = as.numeric(tau >= it.ob["dm_tau"]-(treg*sqrt(it.ob["dm_var.reg"])) &
                                    tau <= it.ob["dm_tau"]+(treg*sqrt(it.ob["dm_var.reg"])))

    #robust SEs
    #this matches what you would get with STATA "robust"
    it.ob["dm_var.robust"] = vcovHC(dmlm, type = "HC1")["Tr","Tr"]
    it.ob["dm_covg.robust"] = as.numeric(tau >= it.ob["dm_tau"]-(treg*sqrt(it.ob["dm_var.robust"])) &
                                           tau <= it.ob["dm_tau"]+(treg*sqrt(it.ob["dm_var.robust"])))

    #c&r SEs "paired"
    yhat = predict(dmlm)
    it.ob["dm_var.cr"] = hajek_var_cr(mod.dat = dat.long, yhat)$vhat_pair
    it.ob["dm_covg.cr"] = as.numeric(tau >= it.ob["dm_tau"]-(tfe*sqrt(it.ob["dm_var.cr"])) &
                                           tau <= it.ob["dm_tau"]+(tfe*sqrt(it.ob["dm_var.cr"])))

    #our variance estimation
    it.ob["dm_var.loop"] = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, n =dat.long$n,
                                  pred = p_simple, loo = F, weighted_imp = T)["varhat"]
    it.ob["dm_covg.loop"] = as.numeric(tau >= it.ob["dm_tau"]-(z*sqrt(it.ob["dm_var.loop"])) &
                                       tau <= it.ob["dm_tau"]+(z*sqrt(it.ob["dm_var.loop"])))

    #=====================  Imai, King, Nall estimator (arithmetic mean weights) ====================#

    #Imai weighted estimate - this is also the same as WLS with cluster weight as the total inds in the pair
    paired <- pair(Y=dat.long$Y, Tr=dat.long$Tr, Z=dat.long$n, P = dat.long$P, n =dat.long$n)

    #imai weighted estimate with arithmethic mean weights
    waest = weighted_effect_est(ordered = paired$ordered, n_ordered = paired$n_ordered,
                                weight = "arithmetic")

    it.ob["ws_tau"] = waest$tau

    #ikn variance
    it.ob["ws_var.ikn"] = waest$vhat_ikn
    it.ob["ws_covg.ikn"] = as.numeric(tau >= it.ob["ws_tau"]-(tfe*sqrt(it.ob["ws_var.ikn"])) &
                                    tau <= it.ob["ws_tau"]+(tfe*sqrt(it.ob["ws_var.ikn"])))

    #WLS equivalent
    dat.ws <- dat.long %>%
      group_by(P) %>%
      mutate(nP = sum(n))

    wslm = lm(Y ~ Tr, data = dat.ws, weights = nP)

    #regression variance
    it.ob["ws_var.reg"] = vcov(wslm)["Tr","Tr"]
    it.ob["ws_covg.reg"] = as.numeric(tau >= it.ob["ws_tau"]-(treg*sqrt(it.ob["ws_var.reg"])) &
                                           tau <= it.ob["ws_tau"]+(treg*sqrt(it.ob["ws_var.reg"])))

    #robust variance
    it.ob["ws_var.robust"] = vcovHC(wslm, type = "HC1")["Tr","Tr"]
    it.ob["ws_covg.robust"] = as.numeric(tau >= it.ob["ws_tau"]-(treg*sqrt(it.ob["ws_var.robust"])) &
                                           tau <= it.ob["ws_tau"]+(treg*sqrt(it.ob["ws_var.robust"])))

    #========================  Fixed Effect Estimator (harmonic mean weights) =======================#
    #regression estimator with fixed effects (Chaisemartin & Ramirez-Cuellar working paper)

    #imai weighted estimate but with harmonic mean weights
    feest = weighted_effect_est(ordered = paired$ordered, n_ordered = paired$n_ordered,
                               weight = "harmonic")

    it.ob["fe_tau"] = feest$tau

    #ikn variance estimator
    it.ob["fe_var.ikn"] = feest$vhat_ikn
    it.ob["fe_covg.ikn"] = as.numeric(tau >= it.ob["fe_tau"]-(tfe*sqrt(it.ob["fe_var.ikn"])) &
                                        tau <= it.ob["fe_tau"]+(tfe*sqrt(it.ob["fe_var.ikn"])))

    #c&r variance estimator
    it.ob["fe_var.cr"] = feest$vhat_cr
    it.ob["fe_covg.cr"] = as.numeric(tau >= it.ob["fe_tau"]-(tfe*sqrt(it.ob["fe_var.cr"])) &
                                           tau <= it.ob["fe_tau"]+(tfe*sqrt(it.ob["fe_var.cr"])))

    #the fe_tau above is equivalent to the coefficient from this estimator
    felm = lm(Y ~ Tr + P, data = dat.long, weights = n)

    #spmk variance estimator
    yhat = predict(felm)
    it.ob["fe_var.spmk"] = spmk_var_est(dat.long, ypred=yhat, v = 0)
    it.ob["fe_covg.spmk"] = as.numeric(tau >= it.ob["fe_tau"]-(tfe*sqrt(it.ob["fe_var.spmk"])) &
                                        tau <= it.ob["fe_tau"]+(tfe*sqrt(it.ob["fe_var.spmk"])))

    #standard wlm variance estimator
    it.ob["fe_var.reg"] = vcov(felm)["Tr","Tr"]
    it.ob["fe_covg.reg"] = as.numeric(tau >= it.ob["fe_tau"]-(tfe*sqrt(it.ob["fe_var.reg"])) &
                                           tau <= it.ob["fe_tau"]+(tfe*sqrt(it.ob["fe_var.reg"])))

    #robust wlm variance estimator
    it.ob["fe_var.robust"] = vcovHC(felm, type = "HC1")["Tr","Tr"]
    it.ob["fe_covg.robust"] = as.numeric(tau >= it.ob["fe_tau"]-(tfe*sqrt(it.ob["fe_var.robust"])) &
                                           tau <= it.ob["fe_tau"]+(tfe*sqrt(it.ob["fe_var.robust"])))


    #======================================= Horvitz-Thompson =======================================#

    #Horvitz Thompson w Middleton & Aronow variance est
    htout <- des_raj_dif(Y = dat.long$Y, Tr = dat.long$Tr, P = dat.long$P, n = dat.long$n, type = "ht")

    it.ob["ht_tau"] = htout$tauhat


    it.ob["ht_var.ma"] = htout$varhat
    it.ob["ht_covg.ma"] = as.numeric(tau >= it.ob["ht_tau"]-(z*sqrt(it.ob["ht_var.ma"])) &
                                    tau <= it.ob["ht_tau"]+(z*sqrt(it.ob["ht_var.ma"])))

    #Su & Ding - OLS equivalent to HT estimator

    dat.long <- dat.long %>%
      ungroup() %>%
      mutate(Ytilde = n*Y*(2*M/N),
             ntilde = n - mean(n))

    htlm <- lm(Ytilde ~ Tr, data = dat.long)

    #general OLS
    it.ob["ht_var.reg"] = vcov(htlm)["Tr","Tr"]
    it.ob["ht_covg.reg"] = as.numeric(tau >= it.ob["ht_tau"]-(treg*sqrt(it.ob["ht_var.reg"])) &
                                           tau <= it.ob["ht_tau"]+(treg*sqrt(it.ob["ht_var.reg"])))
    #robust OLS
    it.ob["ht_var.robust"] = vcovHC(htlm, type = "HC1")["Tr","Tr"]
    it.ob["ht_covg.robust"] = as.numeric(tau >= it.ob["ht_tau"]-(treg*sqrt(it.ob["ht_var.robust"])) &
                                           tau <= it.ob["ht_tau"]+(treg*sqrt(it.ob["ht_var.robust"])))

    #==================================== Des Raj Difference w N ====================================#

    #Middleton & Aronow - des raj difference with only n
    drdout <- des_raj_dif(Y = dat.long$Y, Tr = dat.long$Tr, P = dat.long$P, n = dat.long$n, type = "adj")
    it.ob["drd_tau"] = drdout$tauhat
    it.ob["drd_var.ma"] = drdout$varhat
    it.ob["drd_covg.ma"] = as.numeric(tau >= drdout$tauhat-(z*sqrt(drdout$varhat)) &
                                        tau <= drdout$tauhat+(z*sqrt(drdout$varhat)))

    #=========================================== Su & Ding ==========================================#
    treg2 <- qt(p = 1-(.05/2), df = 2*M-4)

    # recommended Su&Ding estimator for cluster randomized trials
    htnlm <- lm(Ytilde ~ Tr*ntilde, data = dat.long)

    it.ob["htadj_tau"] = coefficients(htnlm)["Tr"]

    #general OLS
    it.ob["htadj_var.reg"] = vcov(htnlm)["Tr","Tr"]
    it.ob["htadj_covg.reg"] = as.numeric(tau >= it.ob["htadj_tau"]-(treg2*sqrt(it.ob["htadj_var.reg"])) &
                                              tau <= it.ob["htadj_tau"]+(treg2*sqrt(it.ob["htadj_var.reg"])))

    #robust OLS
    it.ob["htadj_var.robust"] = vcovHC(htnlm, type = "HC1")["Tr","Tr"]
    it.ob["htadj_covg.robust"] = as.numeric(tau >= it.ob["htadj_tau"]-(treg2*sqrt(it.ob["htadj_var.robust"])) &
                                           tau <= it.ob["htadj_tau"]+(treg2*sqrt(it.ob["htadj_var.robust"])))

    #============================================ P-LOOP ============================================#

    #our point and variance estimators
    lpout = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, n =dat.long$n, pred = p_simple)

    it.ob["p_unadj_tau"] = lpout["tauhat"]
    it.ob["p_unadj_var.loop"] = lpout["varhat"]
    it.ob["p_unadj_covg.loop"] = as.numeric(tau >= lpout["tauhat"]-(z*sqrt(lpout["varhat"])) &
                                         tau <= lpout["tauhat"]+(z*sqrt(lpout["varhat"])))

    #our point and variance estimators using weighted loo imputation
    lpoutw = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, n =dat.long$n, pred = p_simple, weighted_imp = T)

    it.ob["p_unadjw_tau"] = lpoutw["tauhat"]
    it.ob["p_unadjw_var.loop"] = lpoutw["varhat"]
    it.ob["p_unadjw_covg.loop"] = as.numeric(tau >= lpoutw["tauhat"]-(z*sqrt(lpoutw["varhat"])) &
                                              tau <= lpoutw["tauhat"]+(z*sqrt(lpoutw["varhat"])))


    #our point and variance estimators - using full mean imputation, not loo
    lpoutl = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, n =dat.long$n, pred = p_simple, loo = F)

    it.ob["p_unadjmi_tau"] = lpoutl["tauhat"]
    it.ob["p_unadjmi_var.loop"] = lpoutl["varhat"]
    it.ob["p_unadjmi_covg.loop"] = as.numeric(tau >= lpoutl["tauhat"]-(z*sqrt(lpoutl["varhat"])) &
                                                   tau <= lpoutl["tauhat"]+(z*sqrt(lpoutl["varhat"])))

    #our point and variance estimators - using the truth
    lpoutt = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, n =dat.long$n, pred = p_simple, true_sum = trusums)

    it.ob["p_unadjtruth_tau"] = lpoutt["tauhat"]
    it.ob["p_unadjtruth_var.loop"] = lpoutt["varhat"]
    it.ob["p_unadjtruth_covg.loop"] = as.numeric(tau >= lpoutt["tauhat"]-(z*sqrt(lpoutt["varhat"])) &
                                              tau <= lpoutt["tauhat"]+(z*sqrt(lpoutt["varhat"])))

    #our point and variance estimators - using n as a covariate
    lpoutn = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, n =dat.long$n, Z = dat.long$n, pred = p_ols_interp)

    it.ob["p_nadj_tau"] = lpoutn["tauhat"]
    it.ob["p_nadj_var.loop"] = lpoutn["varhat"]
    it.ob["p_nadj_covg.loop"] = as.numeric(tau >= lpoutn["tauhat"]-(z*sqrt(lpoutn["varhat"])) &
                                              tau <= lpoutn["tauhat"]+(z*sqrt(lpoutn["varhat"])))

    it.ob


  }

  return(sim.out)

}


###################################### RUN COV ADJ ESTIMATORS #########################################

adj_est_sims <- function(n.sim = 1000, seed = NULL,
                           dat.orig, po_mt, po_mc, tau){

  # set seed
  if(!is.null(seed)){set.seed(seed)}
  print(paste("The seed is set to:", seed))

  M = nrow(dat.orig)
  N = sum(dat.orig$n1, dat.orig$n2)

  z <- qnorm(p = 1-(.05/2))
  tdm <- qt(p = 1-(.05/2), df = 2*M-2)
  treg <- qt(p = 1-(.05/2), df = 2*M-3)
  tfe <- qt(p = 1-(.05/2), df = M-1)
  tfecov <- qt(p = 1-(.05/2), df = M-2)

  sim.out <- foreach(i = 1:n.sim, .combine=rbind, .errorhandling = "remove") %do% {

    it.ob <- c()

    Tr <- sample(c(0,1), M, replace = T)
    # create data with new treatment effect vector and observed outcomes
    dat <- dat.orig %>%
      mutate(Tr1 = Tr, Tr2 = 1-Tr,
             Y1 = Tr1*(po_mt[,1]) + (1-Tr1)*po_mc[,1],
             Y2 = (1-Tr1)*(po_mt[,2]) + (Tr1)*po_mc[,2]
      )

    #reshape longer as input for functions
    dat.long <- dat %>%
      pivot_longer(cols = Tr1:n2,
                   names_pattern = "(.*)(\\d)",
                   names_to = c(".value", "pair_n"))


    #========================================    Hajek     ========================================#

    dmlm = lm(Y ~ Tr, data = dat.long, weights = n)
    it.ob["dm_tau"] = coefficients(dmlm)["Tr"]
    it.ob["dm_var.robust"] = vcovHC(dmlm, type = "HC1")["Tr","Tr"]
    it.ob["dm_covg.robust"] = as.numeric(tau >= it.ob["dm_tau"]-(tdm*sqrt(it.ob["dm_var.robust"])) &
                                           tau <= it.ob["dm_tau"]+(tdm*sqrt(it.ob["dm_var.robust"])))

    it.ob["dm_r2"] = summary(dmlm)$adj.r.squared

    #===================================     WLS with Covariate    =================================#

    reglm = lm(Y ~ Tr + Z1, data = dat.long, weights = n)
    it.ob["reg_tau"] = coefficients(reglm)["Tr"]
    it.ob["reg_var.robust"] = vcovHC(reglm, type = "HC1")["Tr","Tr"]
    it.ob["reg_covg.robust"] = as.numeric(tau >= it.ob["reg_tau"]-(treg*sqrt(it.ob["reg_var.robust"])) &
                                              tau <= it.ob["reg_tau"]+(treg*sqrt(it.ob["reg_var.robust"])))

    it.ob["reg_r2"] = summary(reglm)$adj.r.squared


    #======================================= Horvitz-Thompson =======================================#

    #Horvitz Thompson w Middleton & Aronow variance est
    htout <- des_raj_dif(Y = dat.long$Y, Tr = dat.long$Tr, P = dat.long$P, n = dat.long$n, type = "ht")

    it.ob["ht_tau"] = htout$tauhat

    it.ob["ht_var.ma"] = htout$varhat
    it.ob["ht_covg.ma"] = as.numeric(tau >= it.ob["ht_tau"]-(z*sqrt(it.ob["ht_var.ma"])) &
                                       tau <= it.ob["ht_tau"]+(z*sqrt(it.ob["ht_var.ma"])))

    #==================================== Des Raj Difference w Cov ====================================#

    #Middleton & Aronow - des raj difference with only n
    drdout <- des_raj_dif(Y = dat.long$Y, Tr = dat.long$Tr, P = dat.long$P, n = dat.long$n, Z = dat.long$Z1, type = "adj")
    it.ob["drd_tau"] = drdout$tauhat
    it.ob["drd_var.ma"] = drdout$varhat
    it.ob["drd_covg.ma"] = as.numeric(tau >= drdout$tauhat-(z*sqrt(drdout$varhat)) &
                                        tau <= drdout$tauhat+(z*sqrt(drdout$varhat)))

    #==================================== Su & Ding w Cov ====================================#

    dat.long <- dat.long %>%
      ungroup() %>%
      mutate(Ytilde = n*Y*(2*M/N),
             ntilde = n - mean(n),
             xtilde = Z1 - mean(Z1))

    treg2 <- qt(p = 1-(.05/2), df = 2*M-4)

    # recommended Su & Ding estimator for cluster randomized trials
    htnlm <- lm(Ytilde ~ Tr + ntilde + xtilde + Tr:ntilde + Tr:xtilde, data = dat.long)

    it.ob["htadj_tau"] = coefficients(htnlm)["Tr"]

    #robust OLS
    it.ob["htadj_var.robust"] = vcovHC(htnlm, type = "HC1")["Tr","Tr"]
    it.ob["htadj_covg.robust"] = as.numeric(tau >= it.ob["htadj_tau"]-(treg2*sqrt(it.ob["htadj_var.robust"])) &
                                              tau <= it.ob["htadj_tau"]+(treg2*sqrt(it.ob["htadj_var.robust"])))


    #========================  Fixed Effect Estimator (harmonic mean weights) =======================#

    felm = lm(Y ~ Tr + P, data = dat.long, weights = n)
    it.ob["fe_tau"] = coefficients(felm)["Tr"]
    it.ob["fe_var.robust"] = vcovHC(felm, type = "HC1")["Tr","Tr"]
    it.ob["fe_covg.robust"] = as.numeric(tau >= it.ob["fe_tau"]-(tfe*sqrt(it.ob["fe_var.robust"])) &
                                           tau <= it.ob["fe_tau"]+(tfe*sqrt(it.ob["fe_var.robust"])))

    it.ob["fe_r2"] = summary(felm)$adj.r.squared

    #======================================  Fixed Effect + COV =====================================#

    felmcov = lm(Y ~ Tr + Z1 + P, data = dat.long, weights = n)
    it.ob["feadj_tau"] =coefficients(felmcov)["Tr"]
    it.ob["feadj_var.robust"] = vcovHC(felmcov, type = "HC1")["Tr","Tr"]
    it.ob["feadj_covg.robust"] = as.numeric(tau >= it.ob["feadj_tau"]-(tfecov*sqrt(it.ob["feadj_var.robust"])) &
                                           tau <= it.ob["feadj_tau"]+(tfecov*sqrt(it.ob["feadj_var.robust"])))

    it.ob["feadj_r2"] = summary(felmcov)$adj.r.squared

    #======================================  Random Effect + COV =====================================#
    relm <- lmer(Y ~ Tr + Z1 + (1|P), data = dat.long, weights = n)

    it.ob["readj_tau"] = fixef(relm)["Tr"]
    it.ob["readj_var.reg"] = vcov(relm)["Tr","Tr"]

    #use the confidence interval calculation built into LMER (does something complicated)
    #confinttemp = suppressMessages(confint(relm))
    # that gave WAY to narrow of intervals... do it by hand instead. Not clear what the degrees of freedom should be, so
    # use the same as ols
    it.ob["readj_covg.reg"] = as.numeric(tau >= it.ob["readj_tau"]-(treg*sqrt(it.ob["readj_var.reg"])) &
                                              tau <= it.ob["readj_tau"]+(treg*sqrt(it.ob["readj_var.reg"])))

    #=====================  Imai, King, Nall estimator (arithmetic mean weights) ====================#

    #Imai weighted estimate - this is also the same as WLS with cluster weight as the total inds in the pair
    paired <- pair(Y=dat.long$Y, Tr=dat.long$Tr, Z=dat.long$n, P = dat.long$P, n =dat.long$n)

    #imai weighted estimate with arithmethic mean weights
    waest = weighted_effect_est(ordered = paired$ordered, n_ordered = paired$n_ordered,
                                weight = "arithmetic")

    it.ob["ws_tau"] = waest$tau

    #ikn variance
    it.ob["ws_var.ikn"] = waest$vhat_ikn
    it.ob["ws_covg.ikn"] = as.numeric(tau >= it.ob["ws_tau"]-(tfe*sqrt(it.ob["ws_var.ikn"])) &
                                        tau <= it.ob["ws_tau"]+(tfe*sqrt(it.ob["ws_var.ikn"])))

   #======================================  PLOOP - no adj =====================================#

    lpout = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, n =dat.long$n, pred = p_simple)

    it.ob["p_unadj_tau"] = lpout["tauhat"]
    it.ob["p_unadj_var.loop"] = lpout["varhat"]
    it.ob["p_unadj_covg.loop"] = as.numeric(tau >= lpout["tauhat"]-(treg*sqrt(lpout["varhat"])) &
                                            tau <= lpout["tauhat"]+(treg*sqrt(lpout["varhat"])))

    #======================================  PLOOP - cov - OLS interp =====================================#

    lpout = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, Z = dat.long$Z1, n =dat.long$n, pred = p_ols_interp)

    it.ob["p_ols_interp_tau"] = lpout["tauhat"]
    it.ob["p_ols_interp_var.loop"] = lpout["varhat"]
    it.ob["p_ols_interp_covg.loop"] = as.numeric(tau >= lpout["tauhat"]-(treg*sqrt(lpout["varhat"])) &
                                                   tau <= lpout["tauhat"]+(treg*sqrt(lpout["varhat"])))

  #======================================  PLOOP - cov - OLS po =====================================#

    lpout = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, Z = dat.long$Z1, n =dat.long$n, pred = p_ols_po)

    it.ob["p_ols_po_tau"] = lpout["tauhat"]
    it.ob["p_ols_po_var.loop"] = lpout["varhat"]
    it.ob["p_ols_po_covg.loop"] = as.numeric(tau >= lpout["tauhat"]-(treg*sqrt(lpout["varhat"])) &
                                                   tau <= lpout["tauhat"]+(treg*sqrt(lpout["varhat"])))


    #======================================  PLOOP - cov - OLS v12 =====================================#

    lpout = p_loop(Y=dat.long$Y, Tr=dat.long$Tr, P = dat.long$P, Z = dat.long$Z1, n =dat.long$n, pred = p_ols_v12)

    it.ob["p_ols_v12_tau"] = lpout["tauhat"]
    it.ob["p_ols_v12_var.loop"] = lpout["varhat"]
    it.ob["p_ols_v12_covg.loop"] = as.numeric(tau >= lpout["tauhat"]-(treg*sqrt(lpout["varhat"])) &
                                               tau <= lpout["tauhat"]+(treg*sqrt(lpout["varhat"])))


    it.ob
  }

  return(sim.out)


}


###################################### FULL SIM WRAPPER #########################################


run_prct_sim <- function(n.dat=100, n.treat= 1000, seed = 276,                          # simulation parameters
                         M = 20,                                                        # number of pairs
                         nimin = 15, nimax = 35, mkimin = 0, mkimax = 0,                # parameters for generating n
                         cov_clustv = 0,                                                # cluster-level variance explained by covariate
                         alpha = 5, pairv = 1, clusterv = 0,                            # parameters for generating intercept of control potential outcomes
                         indv = 1,                                                      # variance of individual error term
                         nbeta = 0,                                                     # parameter for effect of cluster size on control outcome
                         base_tau = 1, tau_var = 0, tau_fslope = 0, tau_fintercept = 0, # parameters for generating tau
                         tau_gslope = 0,
                         loc_shift = 0,                                                 # location shift
                         simfun = all_est_sims                                          # which simulation function to use (estimators)
                         ){

  sim.list <- foreach(i = 1:n.dat) %dopar% {

    # generate cluster data
    clust.dat <- gen_clust_dat(M, nimin, nimax, mkimin, mkimax, cov_clustv,
                               alpha, pairv, clusterv, indv, nbeta, base_tau,
                              tau_var, tau_fslope, tau_fintercept, tau_gslope,
                              loc_shift,
                              seed = seed + i)

    # calculate true ate
    n <- clust.dat$dat %>% select(n1, n2) %>% as.matrix()
    true.ate = sum(clust.dat$tau*n)/sum(n)

    sim.out <- simfun(n.sim = n.treat,
                              seed = seed + i,
                              dat.orig = clust.dat$dat,
                              po_mt = clust.dat$po_t,
                              po_mc = clust.dat$po_c,
                              tau = true.ate)

   list(sim.result = sim.out, tau = true.ate, dat = clust.dat)

  }

  return(sim.list)

}


############################ CLEAN SIMULATION DATA FOR RESULTS ##################################

clean_sim_list <- function(sim.list,
                           methods = c("dm", "ht","ws", "fe",
                                       "drd", "htadj",
                                       "p_unadj","p_unadjw", "p_unadjmi", "p_unadjtruth","p_nadj"),
                           labs = c("hajek", "horvitz-thompson",
                                    "imai", "fixed-effects",
                                    "des raj","adj ht",
                                    "ploop MI","ploop wMI", "ploop no loo MI", "ploop MI truth","ploop n")){


  # create data set with all of the results combined
  sim.dat.full <- foreach(i = 1:length(sim.list), .combine=rbind) %do% {

    sim.list[[i]]$sim.result %>%
      data.frame() %>%
      mutate(data.set = paste0("dat",i))

  }

  # create data set with one row per data generation
  sim.dat <- foreach(i = 1:length(sim.list), .combine=rbind) %do% {

    sim.list[[i]]$sim.result %>%
      data.frame() %>%
      reframe(across(everything(), .fns = list("e" = mean, "v"= var), .names = "{.col}_{.fn}")) %>%
      pivot_longer(everything(),
                   names_to=c("method", "est", "property"),
                   names_pattern = "(\\w*)_(.*)_(e|v)") %>%
      separate(est, into = c("est","est_type"), fill = "right") %>%
      pivot_wider(names_from = property, values_from = value) %>%
      mutate(true.ate = case_when(est == "tau" ~ sim.list[[i]]$tau,
                                  TRUE ~ NA_real_),
             bias = e - true.ate,
             mse = bias^2 + v)

  }

  # re-label estimation methods
  sim.dat$method = factor(sim.dat$method, levels = methods, labels = labs)

  # create more cleaned data for results
  clean.dat <- sim.dat %>%
    group_by(method, est, est_type) %>%
    summarize(across(everything(), ~mean(.x, na.rm = T)))


  return(list(dat.gen.results = sim.dat, results.sum = clean.dat,
              full.results = sim.dat.full))

}




#######################################################################
#solve quadratic formula
# https://dk81.github.io/dkmathstats_site/rmath-quad-formula-r.html

quad <- function(a, b, c) {

  discriminant <- (b^2) - (4*a*c)

  if(discriminant < 0) {
    return(paste0("This quadratic equation has no real numbered roots."))
  }
  else if(discriminant > 0) {
    x_int_plus <- (-b + sqrt(discriminant)) / (2*a)
    x_int_neg <- (-b - sqrt(discriminant)) / (2*a)

    return(list("root1" = x_int_plus,
                "root2" = x_int_neg))
  }
  else #discriminant = 0  case
    x_int <- (-b) / (2*a)
    return(list("root" = x_int))
}
