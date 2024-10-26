############################################################################
# -----------------  R codes for the real data applications -------------- #
# NOTES:
# summarize simulation results 
# please make sure the working directory is the folder where 
# files "FUNs.R" is saved.
# Three real CRC datasets are analyzed
# 1. 4-stream Possums data collected in 1983
# 2. 5-stream Microtus data collected in 1981
# 3. 6-stream snowshare hare data 
# !!!!! Fitting the M*_{4,6} model is computationally expensive,
# the run-time of this R script could be long !!!!!!
############################################################################
############################################################################

# load packages 
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

args <-  commandArgs(trailingOnly = TRUE)
index_s = eval( parse(text=args[1]) )

## **index_s** = 1, 2, 3
## (index of the number of trapping occasions T)

# -------------------- BEGIN REAL DATA ANALYSIS ---------------- #
name.mat <- data.frame("type" = c("Snowshoe_hare", "Mouse", "Possum"),
                       "nstreams" = 6:4)
dat.all.list <- replicate(3, list(), simplify = F)
names(dat.all.list) <- name.mat$type
for(i in 1:3){
  dat.all.list[[i]] <- 
    readxl::read_xlsx(path = "./data/ThreeRealDatasets.xlsx",
                      sheet = names(dat.all.list)[[i]])
}
# the number of posterior samples used for calculating credible intervals
n.post = 500

start = proc.time()[3]
dat.his <- dat.all.list[[name.mat$type[index_s]]]
nstreams = nchar(as.character(dat.his$hist[1]))

hist.sel.num = paste0("1", paste0(rep(0, nstreams - 2), collapse = ""), "1")
hist.sel.den = paste0("1", paste0(rep(0, nstreams - 1), collapse = ""))
hist.key = paste0(paste0(rep(0, nstreams - 1), collapse = ""), "1")

hist.all <- get_allhist(nstreams)
dat.fake <- data.frame("prof" = apply(hist.all, 1, paste, collapse = ""))
dat.full <- left_join(dat.fake, dat.his, by = c("prof" = "hist"))
dat.full$count[is.na(dat.full$count)] <- 0
dat_aggre <- as.data.frame(matrix(dat.full$count, nrow = 1))
colnames(dat_aggre) <- paste0("n", dat.full$prof)

# Fit 1st-order Markov model
nc = sum(dat_aggre)
assump.vec = c("MM1b")
res.yang <- as.data.frame(matrix(NA, ncol = 4, nrow = length(assump.vec)))
colnames(res.yang) <- c("Nhat", "SE", "lci", "uci")
rownames(res.yang) <- assump.vec
for(k in 1:length(assump.vec)){
  assump.use = assump.vec[k] 
  dat.yang.i <- get.dat.yang(dat = dat_aggre, 
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
  res.yang[k, ] <- c(Nhat.yang, sqrt(var.yang), CI.yang)
}

# Fit 1st-order Markov model using the proposed modeling framework
res.yang.alter <- as.data.frame(matrix(NA, ncol = 2, nrow = length(assump.vec)))
colnames(res.yang.alter) <- c("Nhat", "AIC")
rownames(res.yang.alter) <- assump.vec
for(k in 1:length(assump.vec)){
  assump.use = assump.vec[k] 
  res.k <- fit.alter.yang(nstreams = nstreams,
                          dat_aggre = dat_aggre, 
                          assump = assump.use)
  res.yang.alter[k, ] <- c(res.k$Nhat, res.k$AIC)
}
res.yang.all <- list("Yang" = data.frame("nc" = nc, res.yang),
                     "Yang_fitbyalter" = data.frame("nc" = nc, res.yang.alter))


# Fit candidate models M_{T-2,T} and M*_{T-2,T}
max.j = nstreams - 2
if(max.j > 1){
  max.j.1 <- c(paste0(max.j, "days"))
}else{
  max.j.1 <- "1day"
}
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
                                     dat_aggre = dat_aggre,
                                     assump = as.character(model.can$days[i]),
                                     obscons = model.can$cons[i], 
                                     hessian = T)
  model.can$numparm[i] <- re.alter.i$num_parm
  re.mat$Nhat[i] <- re.alter.i$Nhat
  re.mat$pc[i] <- re.alter.i$pc
  re.mat$AIC[i] <- re.alter.i$AIC
  re.can.i[[i]] <- re.alter.i
  cat(i, ",")
} # end loop over all possible candidate models 
re.mat$numparm <- model.can$numparm


# Fit the selected model to get 95% CI and BC estimators
sel.index = which.min(re.mat$AIC)
m.sel <- model.can$model_name[sel.index]

### Extract the selected model and get fitted sel counts
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

prior.name <- c("equprob")
est.name.vec <- c("fitted_bcJ", "fitted_bcB", "fitted_bcC", "fitted_psi")
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
    if(!(m.sel == paste0("F", max.j))){
      ## get 95% CI ##
      re.alter.f1.CI.i <- 
        Dir_CI_alter_fitted(df_aggre = dat_aggre, n.post = n.post,
                            assump = model.can$days[sel.index], 
                            nstreams = nstreams,
                            obscons = model.can$cons[sel.index], 
                            hessian = TRUE,
                            Dir_prior = prior.name,
                            a = a.i, b = b.i,
                            hist.sel.num = hist.sel.num,
                            hist.sel.den = hist.sel.den)
      re.sel.mat$SE[kk] <- sd(re.alter.f1.CI.i$post)
      inter.i <- re.alter.f1.CI.i$CI
      re.sel.mat$lci[kk] <- inter.i[1]
      re.sel.mat$uci[kk] <- inter.i[2]
    }else{
      #### Fmax model using the closed-form MLE ####
      ### computed Nhat based on observed cell counts ###
      re.alter.f1.CI.i <- 
        Dir_CI_alter_Fmax(df_aggre = dat_aggre, n.post = n.post,
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
} # end loop over different BC estimators
re.sel.bcnew <- 
  data.frame("Estimator" = "BC_new",
             "Nhat" = re.sel.mat$Nhat[re.sel.mat$Estimator == "BC_C"],
             "SE" = re.sel.mat$SE[re.sel.mat$Estimator == "BC_C"],
             "lci" = re.sel.mat$lci[re.sel.mat$Estimator == "BC_B"],
             "uci" = re.sel.mat$uci[re.sel.mat$Estimator == "BC_J"],
             "psi" = re.sel.mat$psi[re.sel.mat$Estimator == "BC_C"])
re.sel.mat <- rbind.data.frame(re.sel.mat, re.sel.bcnew)
re.sel.mat <- subset(re.sel.mat, Estimator %in% c("BC_new", "NO_BC"))
re.sel.mat$Estimator <- c("unadjusted", "BC")
re.app.all <- list("Yang" = res.yang.all,
                   "Alter" = re.sel.mat,
                   "Alter_can" = re.mat)
save(re.app.all,
     file = paste0("./inter_res/Res_App_T", nstreams, ".rda"))
proc.time()[3] - start 



