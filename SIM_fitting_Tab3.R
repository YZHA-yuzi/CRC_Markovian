# -------------------  R codes for the simulation study ------------------ #
# NOTES:
# Fitting models to simulated data across various scenarios
# please make sure the working directory is the folder where 
# files "FUNs.R" is saved
# !!!!! Fitting the M*_{4,6} model is computationally expensive,
# the run-time of this R script could be long !!!!!!
############################################################################

## load packages
library(openCR)
library(statmod)
library(bbmle)
library(ars)
library(reshape2)
library(lme4)
library(MCMCpack)
library(bayestestR)
library(dplyr)

## source self-defined R functions
source("FUNs.R")

## create a folder to store simulated datasets
if(!base::dir.exists(paths = "./inter_res")){
  base::dir.create("inter_res")
}


############# SCENARIO III ##############
# simulation results are summarized in Table 3
# data were generated under the M_{4,6}
# 6-stream case: N = 100, 200, pc = 0.85
# Fitted Model:
# (a) 1st-order Markov model
# (b) M_{T-2,T} true model
# (c) AIC selected model, M_{T-2, T} and M*_{T-2,T}
#########################################

args <-  commandArgs(trailingOnly = TRUE)
s_index = eval( parse(text=args[1]) )
B_index = eval( parse(text=args[2]) )
## **s_index** = 1, 2
## (index of simulation scenarios presented in Table 1)
## **B_index** = 1, ..., 1000/A, 
## (where A is the number of simulations ran in one R job)
## e.g., When set A = 100, the total of 1,000 simulation can be
## ran by implementing parallel computation. 
## As a result, 10 jobs are submitted to run simultaneously.  

set.seed(1234 + B_index)

num.s = 6 
# the number of posterior samples used for calculating credible intervals
n.post = 500 
print(paste0("T = ", num.s))
print(paste0("Data mimicing snowshoe hare data were generated based on the model M_{", 
             num.s - 2, ",", num.s, "}"))
print(paste0("# of posterior samples for credible intervals = ", n.post))

## read in simulated datasets
load(paste0("./data/dat_sparse_M", num.s - 2, num.s, ".rda")) # dat.sim.list

nstreams = num.s; prof.mat <- get_allhist(nstreams)
colnames(prof.mat) <- paste0("Y", 1:nstreams)
prof.char <- apply(prof.mat, 1, paste0, collapse = "")
s_index.use <- s_index; dat.sim.all <- dat.sim.list$dat[[s_index.use]]
N.true <- dat.sim.list$parm$N[s_index.use]

print(paste0("Sce = ", s_index.use, 
             "; N = ", N.true, "; pc = ", dat.sim.list$parm$pc[s_index.use],
             "; nc = ", mean(rowSums(dat.sim.all))))
nsims.curr = 1  # "nsims.curr" is the A mentioned before, 
# !!!!!!!! For this specific simulation implemented using this R script, 
# we set A = 1, since fitting the model M*{4,6} is 
# computationally expensive. !!!!!!!! 
B.vec = (1+nsims.curr*(B_index-1)):(nsims.curr*B_index)


# ---------- I. Initialize Vectors/Matrices for Storing Results --------- #
res.all <- list("Yang" = list(), "M_psi" = list(), 
                "All_F" = list(), "AIC_sel" = list())
# Yang = 1st-order Markov model
# M_psi = M_{T-2, T}
# All_F = all candidate models which include M_{T-2,T} and M*_{T-2,T}
# AIC_sel = AIC-selected mdoel

hist.sel.num = paste0("1", paste0(rep(0, nstreams - 2), collapse = ""), "1")
hist.sel.den = paste0("1", paste0(rep(0, nstreams - 1), collapse = ""))
hist.key = paste0(paste0(rep(0, nstreams - 1), collapse = ""), "1")
# ------ I. END Initialize Vectors/Matrices for Storing Results --------- #


# -------------------- II. BEGIN ESTIMATION ---------------- #
count = 1
start = proc.time()[3]
for(isim in B.vec){
  ## 0. prepare data 
  dat.full <- data.frame(prof = prof.char, 
                         count = as.numeric(dat.sim.all[isim, ]))
  if(dat.full$count[(2^nstreams - 1)] == 0){
    dat.full$count[(2^nstreams - 1)] <- 1
  }
  dat.aggre.i <- as.data.frame(matrix(dat.full$count, nrow = 1))
  colnames(dat.aggre.i) <- paste0("n", dat.full$prof)
  
  nc = sum(dat.full$count)
  
  ## 1. fit 1-order Markov model = Yang's model
  assump.vec = c("MM1b")
  col.yang <- c("Model", "Nhat", "SE", "lci", "uci", "Nhat_alter", "AIC_alter")
  res.yang <- as.data.frame(matrix(NA, ncol = length(col.yang), 
                                   nrow = 1:length(assump.vec)))
  colnames(res.yang) <- col.yang
  res.yang[["Model"]] <- assump.vec
  for(k in 1:1:length(assump.vec)){
    assump.use = assump.vec[k] 
    dat.yang.i <- get.dat.yang(dat = dat.aggre.i, 
                               assump = assump.use, 
                               nstreams = nstreams)
    Nhat.yang <- get.Nhat.yang(dat.yang = dat.yang.i, nstreams = nstreams, 
                               nc = nc, assump = assump.use)
    parm.hat.yang <- get.parm.yang(Nhat = Nhat.yang, dat.yang = dat.yang.i,
                                   nstreams = nstreams, nc = nc,
                                   assump = assump.use)
    var.yang <- get.var.yang(Nhat = Nhat.yang, parmhat = parm.hat.yang,
                             nstreams = nstreams, assump = assump.use)
    CI.yang <- get.CI.yang(Nhat = Nhat.yang, varhat = var.yang)
    CI.yang[1] <- max(CI.yang[1], nc)
    res.yang[k, c("Nhat", "SE", "lci", "uci")] <- 
      c(Nhat.yang, sqrt(var.yang), CI.yang)
    
    # fit Yang's model using alternative framework
    res.k <- fit.alter.yang(nstreams = nstreams,
                            dat_aggre = dat.aggre.i, 
                            assump = assump.use)
    res.yang[k, c("Nhat_alter", "AIC_alter")] <- c(res.k$Nhat, res.k$AIC)
  }
  res.yang[["Ntrue"]] <- N.true
  res.yang[["sim"]] <- isim
  res.all$Yang[[count]] <- res.yang
  
  ### 2. Fit F4_cons and F4 model ###
  max.j = nstreams - 2
  max.j.1 <- "4days"
  model.can <- expand.grid(c(T, F), max.j.1)
  colnames(model.can) <- c("cons", "days")
  
  m.name.mat <- data.frame("days" = paste0("F",substr(model.can$days, 1, 1)),
                           "cons" = ifelse(model.can$cons, "_cons", ""))
  m.name.vec <- apply(m.name.mat, 1, paste0, collapse = "")
  model.can[["model_name"]] <- m.name.vec
  
  model.can[["numparm"]] <- rep(NA, nrow(model.can))
  re.mat <- as.data.frame(matrix(NA, ncol = 4, nrow = nrow(model.can)))
  colnames(re.mat) <- c("model", "Nhat", "pc", "AIC")
  re.mat$model <- m.name.vec
  re.can.i <- list()
  for(i in 1:nrow(model.can)){
    re.alter.i <- fit.markovian.assump(nstreams = nstreams, 
                                 dat_aggre = dat.aggre.i,
                                 assump = as.character(model.can$days[i]),
                                 obscons = model.can$cons[i], 
                                 hessian = T)
    model.can$numparm[i] <- re.alter.i$num_parm
    re.mat$Nhat[i] <- re.alter.i$Nhat
    re.mat$pc[i] <- re.alter.i$pc
    re.mat$AIC[i] <- re.alter.i$AIC
    re.can.i[[i]] <- re.alter.i
  } # end loop over all possible candidate models 
  re.mat$numparm <- model.can$numparm
  re.mat[["Ntrue"]] <- N.true
  re.mat[["sim"]] <- isim
  res.all[["All_F"]][[count]] <- re.mat
  
  ### Fit the selected model to get 95% CI and BC estimators
  sel.index <- which.min(re.mat$AIC)
  m.sel <- model.can$model_name[sel.index]
  
  ### Step 3: Extract the selected model and get fitted sel counts
  re.alter.sel <- re.can.i[[sel.index]]
  ## GET fitted cell counts ###
  est.mat.i <- data.frame("parm_name_curr" = names(re.alter.sel$parm),
                          "value" = as.numeric(re.alter.sel$parm))
  parm.mat.sim <- left_join(re.alter.sel$parm.name, est.mat.i,
                            by = c("parm_name_curr" = "parm_name_curr"))
  ## 3. convert parameters to capture probabilities
  parm.mat.sim.use <- parm.mat.sim[,c("ID", "S_curr", "pre_hist",
                                      "parm_name", "value")]
  prob.true.sim <- f_parm_prob(nstreams = nstreams,
                               parm.mat = parm.mat.sim.use)
  nfitted.i <- re.alter.sel$Nhat*prob.true.sim$prob
  dat.full[["count_fitted"]] <- nfitted.i[-c(2^nstreams)]
  
  if(m.sel == paste0("F", nstreams-2)){
    ### computed Nhat based on observed cell counts ###
    n.num.i <- dat.full$count[dat.full$prof == hist.sel.num]
    n.den.i <- dat.full$count[dat.full$prof == hist.sel.den]
    
    psi.i <- n.num.i/(n.num.i + n.den.i)
    psi.bc.i <- (n.num.i + 1)/(1 + n.num.i + n.den.i)
    psi.bcJ.i <- (n.num.i + 0.5)/(1 + n.num.i + n.den.i)
    psi.bc.comb.i <- (n.num.i + 0.75)/(1 + n.num.i + n.den.i)
    
    psi.vec.all <- c(psi.bcJ.i, psi.bc.i, psi.bc.comb.i, psi.i)
    nc.fitted <- sum(dat.full$count)
    nlast.fitted <- dat.full$count[dat.full$prof == hist.key]
  }else{
    ### computed Nhat based on fitted cell counts ###
    n.num.i <- dat.full$count_fitted[dat.full$prof == hist.sel.num]
    n.den.i <- dat.full$count_fitted[dat.full$prof == hist.sel.den]
    
    psi.i <- n.num.i/(n.num.i + n.den.i)
    psi.bc.i <- (n.num.i + 1)/(1 + n.num.i + n.den.i)
    psi.bcJ.i <- (n.num.i + 0.5)/(1 + n.num.i + n.den.i)
    psi.bc.comb.i <- (n.num.i + 0.75)/(1 + n.num.i + n.den.i)
    
    psi.vec.all <- c(psi.bcJ.i, psi.bc.i, psi.bc.comb.i, psi.i)
    nc.fitted <- sum(dat.full$count_fitted)
    nlast.fitted <- dat.full$count_fitted[dat.full$prof == hist.key]
  }
  
  prior.name <- c("equprob")
  est.name.vec <- c("fitted_bcJ", "fitted_bcB", "fitted_bcC", "fitted_psi")
  est.name.vec.forci <- c("fitted_bcJ", "fitted_bcB")
  re.sel.mat <- as.data.frame(matrix(NA, ncol = 5, 
                                     nrow = length(est.name.vec)))
  colnames(re.sel.mat) <- c("Estimator", "Nhat", "SE", "lci", "uci")
  re.sel.mat[["psi"]] <- psi.vec.all[1:length(est.name.vec)]
  comp.ci = T
  for(kk in 1:length(est.name.vec)){
    est.name.curr <- est.name.vec[kk]
    if(est.name.curr == "fitted_bcB"){
      re.sel.mat$Estimator[kk] <- "BC_B"
      re.sel.mat$Nhat[kk] <- 
        get_Nhat_psi(nc = nc.fitted, nlast = nlast.fitted, psi = psi.bc.i)
      a.i = 1; b.i = 0
    }else if(est.name.curr == "fitted_bcJ"){
      re.sel.mat$Estimator[kk] <- "BC_J"
      re.sel.mat$Nhat[kk] <- 
        get_Nhat_psi(nc = nc.fitted, nlast = nlast.fitted, psi = psi.bcJ.i)
      a.i = 0.5; b.i = 0.5
    }else if(est.name.curr == "fitted_psi"){
      re.sel.mat$Estimator[kk] <- "NO_BC"
      re.sel.mat$Nhat[kk] <- 
        get_Nhat_psi(nc = nc.fitted, nlast = nlast.fitted, psi = psi.i)
      a.i = 0; b.i = 0
    }else if(est.name.curr == "fitted_bcC"){
      re.sel.mat$Estimator[kk] <- "BC_C"
      re.sel.mat$Nhat[kk] <- 
        get_Nhat_psi(nc = nc.fitted, nlast = nlast.fitted, psi = psi.bc.comb.i)
      a.i = 0.75; b.i = 0.25
    }
    if(comp.ci){
      if(est.name.curr %in% est.name.vec.forci){
        if(!(m.sel == paste0("F", max.j))){
          ## get 95% CI ##
          # start = proc.time()[3]
          re.alter.f1.CI.i <- 
            Dir_CI_alter_fitted(df_aggre = dat.aggre.i, n.post = n.post,
                                assump = as.character(model.can$days[sel.index]), 
                                nstreams = nstreams,
                                obscons = model.can$cons[sel.index], 
                                hessian = TRUE,
                                Dir_prior = prior.name,
                                a = a.i, b = b.i,
                                hist.sel.num = hist.sel.num,
                                hist.sel.den = hist.sel.den)
          # proc.time()[3] - start
          re.sel.mat$SE[kk] <- sd(re.alter.f1.CI.i$post)
          inter.i <- re.alter.f1.CI.i$CI
          re.sel.mat$lci[kk] <- inter.i[1]
          re.sel.mat$uci[kk] <- inter.i[2]
        }else{
          ## get 95% CI for the M_{T-2,T} model ##
          ## No need to compute fitted cell counts
          re.alter.f1.CI.i <- 
            Dir_CI_alter_Fmax(df_aggre = dat.aggre.i, n.post = n.post,
                              nstreams = nstreams,
                              Dir_prior = prior.name,
                              a = a.i, b = b.i,
                              hist.sel.num = hist.sel.num,
                              hist.sel.den = hist.sel.den)
          re.sel.mat$SE[kk] <- sd(re.alter.f1.CI.i$post)
          inter.i <- re.alter.f1.CI.i$CI
          re.sel.mat$lci[kk] <- inter.i[1]
          re.sel.mat$uci[kk] <- inter.i[2]
        }
      } # END compute 95% CI
    } # END compute 95% CI
  } # end loop over different BC estimators
  re.sel.mat[["Ntrue"]] <- N.true
  re.sel.mat[["sim"]] <- isim
  re.sel.mat[["Model"]] <- m.sel
  res.all[["AIC_sel"]][[count]] <- re.sel.mat
  
  ### Fit true model ####
  if(m.sel == "F4"){
    n.num.i <- dat.full$count[dat.full$prof == hist.sel.num]
    n.den.i <- dat.full$count[dat.full$prof == hist.sel.den]
    
    psi.i <- n.num.i/(n.num.i + n.den.i)
    psi.bc.i <- (n.num.i + 1)/(1 + n.num.i + n.den.i)
    psi.bcJ.i <- (n.num.i + 0.5)/(1 + n.num.i + n.den.i)
    psi.bc.comb.i <- (n.num.i + 0.75)/(1 + n.num.i + n.den.i)
    psi.vec.all <- c(psi.bc.comb.i, psi.i)
    
    nc.fitted <- sum(dat.full$count)
    nlast.fitted <- dat.full$count[dat.full$prof == hist.key]
    prior.name <- c("equprob")
    est.name.vec <- c("fitted_bcC", "fitted_psi")
    re.true.mat <- as.data.frame(matrix(NA, ncol = 5, 
                                        nrow = length(est.name.vec)))
    colnames(re.true.mat) <- c("Estimator", "Nhat", "SE", "lci", "uci")
    re.true.mat[["psi"]] <- psi.vec.all[1:length(est.name.vec)]
    comp.ci = T
    for(kk in 1:length(est.name.vec)){
      est.name.curr <- est.name.vec[kk]
      if(est.name.curr == "fitted_bcB"){
        re.true.mat$Estimator[kk] <- "BC_B"
        re.true.mat$Nhat[kk] <- 
          get_Nhat_psi(nc = nc.fitted, nlast = nlast.fitted, psi = psi.bc.i)
        a.i = 1; b.i = 0
      }else if(est.name.curr == "fitted_bcJ"){
        re.true.mat$Estimator[kk] <- "BC_J"
        re.true.mat$Nhat[kk] <- 
          get_Nhat_psi(nc = nc.fitted, nlast = nlast.fitted, psi = psi.bcJ.i)
        a.i = 0.5; b.i = 0.5
      }else if(est.name.curr == "fitted_psi"){
        re.true.mat$Estimator[kk] <- "NO_BC"
        re.true.mat$Nhat[kk] <- 
          get_Nhat_psi(nc = nc.fitted, nlast = nlast.fitted, psi = psi.i)
        a.i = 0; b.i = 0
      }else if(est.name.curr == "fitted_bcC"){
        re.true.mat$Estimator[kk] <- "BC_C"
        re.true.mat$Nhat[kk] <- 
          get_Nhat_psi(nc = nc.fitted, nlast = nlast.fitted, psi = psi.bc.comb.i)
        a.i = 0.75; b.i = 0.25
      }
      if(comp.ci){
        #### Fit the proposed model using the closed-form MLE ####
        re.alter.f1.CI.i <- 
          Dir_CI_alter_Fmax(df_aggre = dat.aggre.i, n.post = n.post,
                            nstreams = nstreams,
                            Dir_prior = prior.name,
                            a = a.i, b = b.i,
                            hist.sel.num = hist.sel.num,
                            hist.sel.den = hist.sel.den)
        re.true.mat$SE[kk] <- sd(re.alter.f1.CI.i$post)
        inter.i <- re.alter.f1.CI.i$CI
        re.true.mat$lci[kk] <- inter.i[1]
        re.true.mat$uci[kk] <- inter.i[2]
      } # END compute 95% CI
    } # end loop over different BC estimators
    re.true.mat <- rbind.data.frame(re.sel.mat[1:2,colnames(re.true.mat)], 
                                    re.true.mat)
    re.true.mat[["Ntrue"]] <- N.true
    re.true.mat[["sim"]] <- isim
    re.true.mat[["Model"]] <- "psi"
  }else if(m.sel == "F4_cons"){
    
    n.num.i <- dat.full$count[dat.full$prof == hist.sel.num]
    n.den.i <- dat.full$count[dat.full$prof == hist.sel.den]
    
    psi.i <- n.num.i/(n.num.i + n.den.i)
    psi.bc.i <- (n.num.i + 1)/(1 + n.num.i + n.den.i)
    psi.bcJ.i <- (n.num.i + 0.5)/(1 + n.num.i + n.den.i)
    psi.bc.comb.i <- (n.num.i + 0.75)/(1 + n.num.i + n.den.i)
    psi.vec.all <- c(psi.bcJ.i, psi.bc.i, psi.bc.comb.i, psi.i)
    
    nc.fitted <- sum(dat.full$count)
    nlast.fitted <- dat.full$count[dat.full$prof == hist.key]
    prior.name <- c("equprob")
    est.name.vec <- c("fitted_bcJ", "fitted_bcB", "fitted_bcC", "fitted_psi")
    re.true.mat <- as.data.frame(matrix(NA, ncol = 5, 
                                        nrow = length(est.name.vec)))
    colnames(re.true.mat) <- c("Estimator", "Nhat", "SE", "lci", "uci")
    re.true.mat[["psi"]] <- psi.vec.all[1:length(est.name.vec)]
    comp.ci = T
    for(kk in 1:length(est.name.vec)){
      est.name.curr <- est.name.vec[kk]
      if(est.name.curr == "fitted_bcB"){
        re.true.mat$Estimator[kk] <- "BC_B"
        re.true.mat$Nhat[kk] <- 
          get_Nhat_psi(nc = nc.fitted, nlast = nlast.fitted, psi = psi.bc.i)
        a.i = 1; b.i = 0
      }else if(est.name.curr == "fitted_bcJ"){
        re.true.mat$Estimator[kk] <- "BC_J"
        re.true.mat$Nhat[kk] <- 
          get_Nhat_psi(nc = nc.fitted, nlast = nlast.fitted, psi = psi.bcJ.i)
        a.i = 0.5; b.i = 0.5
      }else if(est.name.curr == "fitted_psi"){
        re.true.mat$Estimator[kk] <- "NO_BC"
        re.true.mat$Nhat[kk] <- 
          get_Nhat_psi(nc = nc.fitted, nlast = nlast.fitted, psi = psi.i)
        a.i = 0; b.i = 0
      }else if(est.name.curr == "fitted_bcC"){
        re.true.mat$Estimator[kk] <- "BC_C"
        re.true.mat$Nhat[kk] <- 
          get_Nhat_psi(nc = nc.fitted, nlast = nlast.fitted, psi = psi.bc.comb.i)
        a.i = 0.75; b.i = 0.25
      }
      if(comp.ci){
        #### Fit the proposed model using the closed-form MLE ####
        # start = proc.time()[3]
        re.alter.f1.CI.i <- 
          Dir_CI_alter_Fmax(df_aggre = dat.aggre.i, n.post = n.post,
                            nstreams = nstreams,
                            Dir_prior = prior.name,
                            a = a.i, b = b.i,
                            hist.sel.num = hist.sel.num,
                            hist.sel.den = hist.sel.den)
        # proc.time()[3] - start
        re.true.mat$SE[kk] <- sd(re.alter.f1.CI.i$post)
        inter.i <- re.alter.f1.CI.i$CI
        re.true.mat$lci[kk] <- inter.i[1]
        re.true.mat$uci[kk] <- inter.i[2]
      } # END compute 95% CI
    } # end loop over different BC estimators
    
    re.true.mat[["Ntrue"]] <- N.true
    re.true.mat[["sim"]] <- isim
    re.true.mat[["Model"]] <- "psi"
  }
  
  re.true.mat.1 <- re.true.mat
  re.true.mat.1 <- re.true.mat.1[,!colnames(re.true.mat.1) %in% c("SE", "psi")]
  re.psi.j <- subset(re.true.mat.1, Estimator == "BC_J")
  re.psi.b <- subset(re.true.mat.1, Estimator == "BC_B")
  re.psi.bj <- subset(re.true.mat.1, Estimator == "BC_C")
  re.psi.bj$uci <- re.psi.j$uci
  re.psi.bj$lci <- re.psi.b$lci
  re.true.mat.1[re.true.mat.1$Estimator == "BC_C", ] <- re.psi.bj
  re.true.mat.1 <- subset(re.true.mat.1, Estimator %in% c("BC_C", "NO_BC"))
  re.true.mat.1$Estimator <- c("BC", "unadjusted")
  
  res.all[["M_psi"]][[count]] <- re.true.mat.1
  count = count + 1
  # cat(count, ",")
}
proc.time()[3] - start
save(res.all, file = paste0("./inter_res/Res_sparse_M", num.s-2, num.s, 
                            "_N", N.true, "_B", B_index, ".rda"))


