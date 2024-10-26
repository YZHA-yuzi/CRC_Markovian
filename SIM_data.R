############################################################################
# -------------------  R codes for the simulation study ------------------ #
# NOTES:
# generating simulated datasets
# please make sure the working directory is the folder where 
# files "FUNs.R" 
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
if(!base::dir.exists(paths = "./data")){
  base::dir.create("data")
}

############# SCENARIO I ##############
# simulation results are summarized in Table 1 
# data were generated under the M_{T-2,T}, for T = 3, 4
# 3-stream case: N = 200, 500, pc = 0.541
# 4-stream case: N = 200, 500, pc = 0.571
# Fitted Model:
# (a) 1st-order Markov model
# (b) M_{T-2,T} true model
# (c) AIC selected model, M_{T-2, T} and M*_{T-2,T}
#########################################
set.seed(1234)
# the number of trapping occasions, T = 3, 4 
nstreams.vec <- c(3, 4)
# the true population size N = 200, 500
N.true.vec <- c(200, 500)
# the total  number of simulations
nsims = 1000
for(nstreams in nstreams.vec){
  if(nstreams == 3){
    count = 1; dat.list <- list()
    
    psi.vec <- 0.1; psi = psi.vec[1]
    parm.names <- c("p1", "p21", "p21bar",
                    "p312", "p312bar", "p31bar2", "p31bar2bar")
    parm.setup <- as.data.frame(matrix(NA, nrow = length(psi.vec)*
                                         length(N.true.vec),
                           ncol = 2^nstreams - 1 + 3))
    colnames(parm.setup) <- c("Sce", "N", parm.names, "pc")
    
    p1 = 0.4; p21 = 0.4; p21bar = 0.15; 
    p312 = 0.35; p31bar2 = 0.25
    # the M_{1, 3} model imposes only one untestable constraint: 
    # p31bar2bar = p312bar 
    p312bar = psi
    parm.vec <- c(p1, p21, p21bar, p312, p312bar, p31bar2, psi)
    # compute p_{hi} based on specified parameters
    prob.vec.true <- comp.prob(parm.vec)
    
    for(N in N.true.vec){
      
      parm.setup$Sce[count] <- paste0("S", count)
      parm.setup$N[count] <- N
      parm.setup[count, parm.names] <- parm.vec
      parm.setup$pc[count] <- (1 - prob.vec.true[8])
      
      # generate data based on the multinominal dist
      dat.full.sim <- t(rmultinom(nsims, size = N,
                                  prob = as.numeric(prob.vec.true)))
      dat.full.sim <- as.data.frame(dat.full.sim)
      
      prof.full.mat <- expand.grid(rep(list(c(1,0)), nstreams))
      for(i in nstreams:1){
        prof.full.mat = prof.full.mat[order(prof.full.mat[,i],
                                            decreasing = T), ]
      }
      colnames(prof.full.mat) <- paste0("X", 1:nstreams)
      prof.char.vec <- paste0("n",
                              apply(prof.full.mat[-(2^nstreams),],
                                    1, paste, collapse = ""))
      dat.sim <- dat.full.sim[,-(2^nstreams)]
      colnames(dat.sim) <- prof.char.vec
      dat.list[[count]] <- dat.sim
      count = count + 1
    } # end loop over N
    dat.sim.list <- list("dat" = dat.list,
                         "parm" = parm.setup)
    save(dat.sim.list, file = paste0("./data/dat_M13.rda"))
  }else if(nstreams == 4){
    
    ## get some estimated parameters from real data 
    y.sel = "1983" # 1983, 4 trapping occasions
    dat <- openCR::FebpossumCH[[y.sel]][,,1]
    dat.his <- apply(dat, 1, function(x) base::paste(x, collapse = ""))
    dat.his <- as.data.frame(table(dat.his))
    colnames(dat.his) <- c("prof", "count")
    hist.all <- get_allhist(nstreams)
    dat.fake <- data.frame("prof" = apply(hist.all, 1, paste, collapse = ""))
    dat.full <- left_join(dat.fake, dat.his, by = c("prof" = "prof"))
    dat.full$count[is.na(dat.full$count)] <- 0
    dat_aggre <- as.data.frame(matrix(dat.full$count, nrow = 1))
    colnames(dat_aggre) <- paste0("n", dat.full$prof)
    re.alter.i <- fit.markovian.assump(nstreams = nstreams, 
                                 dat_aggre = dat_aggre,
                                 assump = "2days", obscons = F, hessian = T)
    est.mat.i <- data.frame("parm_name_curr" = names(re.alter.i$parm),
                            "value" = as.numeric(re.alter.i$parm))
    parm.mat.sim <- left_join(re.alter.i$parm.name, est.mat.i,
                              by = c("parm_name_curr" = "parm_name_curr"))
    parm.mat.sim.use <- parm.mat.sim[,c("ID", "S_curr", "pre_hist",
                                        "parm_name", "value")]
    
    count = 1; dat.list <- list()
    psi.vec = c(0.3)
    psi = psi.vec[1]
    parm.setup <- 
      as.data.frame(matrix(NA, nrow = length(psi.vec)*length(N.true.vec),
                           ncol = 2^nstreams - 1 + 3))
    parm.names <- parm.mat.sim.use$parm_name
    colnames(parm.setup) <- c("Sce", "N", parm.names, "pc")
    probhist.setup <- 
      as.data.frame(matrix(NA, nrow = length(psi.vec)*length(N.true.vec),
                           ncol = 2^nstreams))
    colnames(probhist.setup) <- 
      paste0("p", apply(get_allhist_all(4), 1, paste, collapse = ""))
    
    for(N in N.true.vec){
      
      ## modify estimated parameters 
      parm.mat.sim.use.1 <- parm.mat.sim.use
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$value > 0.6] <-
        parm.mat.sim.use.1$value[parm.mat.sim.use.1$value > 0.6]*0.5
      parm.mat.sim.use.1$value <- round(parm.mat.sim.use.1$value, 2)
      
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p1"] <- 0.15
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p2_0"] <- 0.18
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p3_00"] <- 0.12
      
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p4_010"] <- 0.45
      
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p2_1"] <- 0.55
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p3_11"] <- 0.35
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p3_01"] <- 0.15
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p4_111"] <- 0.4
      
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p4_000"] <- psi
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p4_100"] <- psi
      
      ## convert parameters to capture probabilities
      prob.true.sim <- f_parm_prob(nstreams = nstreams,
                                   parm.mat = parm.mat.sim.use.1)
      parm.setup$Sce[count] <- paste0("S", count)
      parm.setup$N[count] <- N
      parm.setup[count, parm.names] <- parm.mat.sim.use.1$value
      parm.setup$pc[count] <- (1 - prob.true.sim$prob[2^nstreams])
      probhist.setup[count, ] <- as.numeric(prob.true.sim$prob)
      dat.full.sim <- t(rmultinom(nsims, size = N,
                                  prob = as.numeric(prob.true.sim$prob)))
      dat.full.sim <- as.data.frame(dat.full.sim)
      
      prof.full.mat <- expand.grid(rep(list(c(1,0)), nstreams))
      for(i in nstreams:1){
        prof.full.mat = prof.full.mat[order(prof.full.mat[,i],
                                            decreasing = T), ]
      }
      colnames(prof.full.mat) <- paste0("X", 1:nstreams)
      prof.char.vec <- paste0("n",
                              apply(prof.full.mat[-(2^nstreams),],
                                    1, paste, collapse = ""))
      dat.sim <- dat.full.sim[,-(2^nstreams)]
      colnames(dat.sim) <- prof.char.vec
      dat.list[[count]] <- dat.sim
      count = count + 1
    } # loop over N
    probhist.setup <- data.frame("Sce" = parm.setup$Sce,
                                 probhist.setup)
    dat.sim.list <- list("dat" = dat.list,
                         "parm" = parm.setup,
                         "probhist" = probhist.setup)
    save(dat.sim.list, file = paste0("./data/dat_M24.rda"))
  }
} # end loop over T 




############# SCENARIO II ##############
# simulation results are summarized in Table 2
# data were generated under the 1st-order Markov model, for T = 3, 4
# 3-stream case: N = 200, 500, pc = 0.488
# 4-stream case: N = 200, 500, pc = 0.590
# Fitted Model:
# (a) 1st-order Markov model
# (b) M_{T-2,T} true model
# (c) AIC selected model, M_{T-2, T} and M*_{T-2,T}
#########################################
set.seed(1234)
# the number of trapping occasions, T = 3, 4 
nstreams.vec <- c(3, 4)
# the true population size N = 200, 500
N.true.vec <- c(200, 500)
# the total  number of simulations
nsims = 1000
for(nstreams in nstreams.vec){
  if(nstreams == 3){
    count = 1; dat.list <- list()
    
    psi.vec <- 0.2; psi = psi.vec[1]
    parm.names <- c("p1", "p21", "p21bar",
                    "p312", "p312bar", "p31bar2", "p31bar2bar")
    parm.setup <- as.data.frame(matrix(NA, nrow = length(psi.vec)*
                                         length(N.true.vec),
                                       ncol = 2^nstreams - 1 + 3))
    colnames(parm.setup) <- c("Sce", "N", parm.names, "pc")
    
    p1 = psi; p21bar = psi; p312bar = psi
    p21 = 0.4; 
    p312 = p21; p31bar2 = p21
    
    parm.vec <- c(p1, p21, p21bar, p312, p312bar, p31bar2, psi)
    # compute p_{hi} based on specified parameters
    prob.vec.true <- comp.prob(parm.vec)
    
    for(N in N.true.vec){
      
      parm.setup$Sce[count] <- paste0("S", count)
      parm.setup$N[count] <- N
      parm.setup[count, parm.names] <- parm.vec
      parm.setup$pc[count] <- (1 - prob.vec.true[8])
      
      # generate data based on the multinominal dist
      dat.full.sim <- t(rmultinom(nsims, size = N,
                                  prob = as.numeric(prob.vec.true)))
      dat.full.sim <- as.data.frame(dat.full.sim)
      
      prof.full.mat <- expand.grid(rep(list(c(1,0)), nstreams))
      for(i in nstreams:1){
        prof.full.mat = prof.full.mat[order(prof.full.mat[,i],
                                            decreasing = T), ]
      }
      colnames(prof.full.mat) <- paste0("X", 1:nstreams)
      prof.char.vec <- paste0("n",
                              apply(prof.full.mat[-(2^nstreams),],
                                    1, paste, collapse = ""))
      dat.sim <- dat.full.sim[,-(2^nstreams)]
      colnames(dat.sim) <- prof.char.vec
      dat.list[[count]] <- dat.sim
      count = count + 1
    } # end loop over N
    dat.sim.list <- list("dat" = dat.list,
                         "parm" = parm.setup)
    save(dat.sim.list, file = paste0("./data/dat_T3_MarkovModel.rda"))
  }else if(nstreams == 4){
    
    ## get some estimated parameters from real data 
    y.sel = "1983" # 1983, 4 trapping occasions
    dat <- openCR::FebpossumCH[[y.sel]][,,1]
    dat.his <- apply(dat, 1, function(x) base::paste(x, collapse = ""))
    dat.his <- as.data.frame(table(dat.his))
    colnames(dat.his) <- c("prof", "count")
    hist.all <- get_allhist(nstreams)
    dat.fake <- data.frame("prof" = apply(hist.all, 1, paste, collapse = ""))
    dat.full <- left_join(dat.fake, dat.his, by = c("prof" = "prof"))
    dat.full$count[is.na(dat.full$count)] <- 0
    dat_aggre <- as.data.frame(matrix(dat.full$count, nrow = 1))
    colnames(dat_aggre) <- paste0("n", dat.full$prof)
    re.alter.i <- fit.markovian.assump(nstreams = nstreams, 
                                 dat_aggre = dat_aggre,
                                 assump = "2days", obscons = F, hessian = T)
    est.mat.i <- data.frame("parm_name_curr" = names(re.alter.i$parm),
                            "value" = as.numeric(re.alter.i$parm))
    parm.mat.sim <- left_join(re.alter.i$parm.name, est.mat.i,
                              by = c("parm_name_curr" = "parm_name_curr"))
    parm.mat.sim.use <- parm.mat.sim[,c("ID", "S_curr", "pre_hist",
                                        "parm_name", "value")]
    
    count = 1; dat.list <- list()
    psi.vec = c(0.2)
    psi = psi.vec[1]
    parm.setup <- 
      as.data.frame(matrix(NA, nrow = length(psi.vec)*length(N.true.vec),
                           ncol = 2^nstreams - 1 + 3))
    parm.names <- parm.mat.sim.use$parm_name
    colnames(parm.setup) <- c("Sce", "N", parm.names, "pc")
    probhist.setup <- 
      as.data.frame(matrix(NA, nrow = length(psi.vec)*length(N.true.vec),
                           ncol = 2^nstreams))
    colnames(probhist.setup) <- 
      paste0("p", apply(get_allhist_all(4), 1, paste, collapse = ""))
    
    for(N in N.true.vec){
      ## modify estimated parameters 
      parm.mat.sim.use.1 <- parm.mat.sim.use
      
      grp.1 <- c("p1", "p2_0",
                 paste0("p3_", sapply(get_allhist_all(1), 
                                      paste0, collapse = ""), "0"),
                 paste0("p4_", apply(get_allhist_all(2), 
                                     1, paste0, collapse = ""), "0"))
      grp.2 <- c("p2_1",
                 paste0("p3_", sapply(get_allhist_all(1), 
                                      paste0, collapse = ""), "1"),
                 paste0("p4_", apply(get_allhist_all(2), 
                                     1, paste0, collapse = ""), "1"))
      
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name %in% grp.1] <- 
        psi
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name %in% grp.2] <- 
        0.35
      
      ## convert parameters to capture probabilities
      prob.true.sim <- f_parm_prob(nstreams = nstreams,
                                   parm.mat = parm.mat.sim.use.1)
      parm.setup$Sce[count] <- paste0("S", count)
      parm.setup$N[count] <- N
      parm.setup[count, parm.names] <- parm.mat.sim.use.1$value
      parm.setup$pc[count] <- (1 - prob.true.sim$prob[2^nstreams])
      probhist.setup[count, ] <- as.numeric(prob.true.sim$prob)
      dat.full.sim <- t(rmultinom(nsims, size = N,
                                  prob = as.numeric(prob.true.sim$prob)))
      dat.full.sim <- as.data.frame(dat.full.sim)
      
      prof.full.mat <- expand.grid(rep(list(c(1,0)), nstreams))
      for(i in nstreams:1){
        prof.full.mat = prof.full.mat[order(prof.full.mat[,i],
                                            decreasing = T), ]
      }
      colnames(prof.full.mat) <- paste0("X", 1:nstreams)
      prof.char.vec <- paste0("n",
                              apply(prof.full.mat[-(2^nstreams),],
                                    1, paste, collapse = ""))
      dat.sim <- dat.full.sim[,-(2^nstreams)]
      colnames(dat.sim) <- prof.char.vec
      dat.list[[count]] <- dat.sim
      count = count + 1
      
    } # loop over N
    probhist.setup <- data.frame("Sce" = parm.setup$Sce,
                                 probhist.setup)
    dat.sim.list <- list("dat" = dat.list,
                         "parm" = parm.setup,
                         "probhist" = probhist.setup)
    save(dat.sim.list, file = paste0("./data/dat_T4_MarkovModel.rda"))
  }
} # end loop over T 




############# SCENARIO III ##############
# simulation results are summarized in Table 3
# data were generated under the M_{4,6}
# 6-stream case: N = 100, 200, pc = 0.85
# Fitted Model:
# (a) 1st-order Markov model
# (b) M_{T-2,T} true model
# (c) AIC selected model, M_{T-2, T} and M*_{T-2,T}
#########################################
# NOTES:
# To mimic the snowshoe hare data analyzed in the manuscript, 
# we first obtain parameters by fitting the proposed M_{4,6} model
# to the snowshoe hare data

### read in snowshoe hare data 
dat.full <- readxl::read_xlsx(path = "./data/ThreeRealDatasets.xlsx", 
                              sheet = "Snowshoe_hare")
dat_aggre <- as.data.frame(matrix(dat.full$count, nrow = 1))
colnames(dat_aggre) <- paste0("n", dat.full$hist)

nstreams = 6
re.alter.i <- fit.markovian.assump(nstreams = nstreams, 
                             dat_aggre = dat_aggre,
                             assump = "4days", obscons = F, hessian = T)
est.mat.i <- data.frame("parm_name_curr" = names(re.alter.i$parm),
                        "value" = as.numeric(re.alter.i$parm))
parm.mat.sim <- left_join(re.alter.i$parm.name, est.mat.i,
                          by = c("parm_name_curr" = "parm_name_curr"))
parm.mat.sim.use <- parm.mat.sim[,c("ID", "S_curr", "pre_hist",
                                    "parm_name", "value")]
## II. (1) generate data under the assumption psi = p612bar3bar4bar5bar 
set.seed(1234)
nstreams = 6
nsims = 1000
N.true.vec = c(100, 200)
parm.setup <- 
  as.data.frame(matrix(NA, nrow = 1*length(N.true.vec),
                       ncol = 2^nstreams - 1 + 3))
parm.names <- parm.mat.sim.use$parm_name
colnames(parm.setup) <- c("Sce", "N", parm.names, "pc")
dat.list <- list()
count = 1
for(N in N.true.vec){
  ## modify estimated parameters 
  parm.mat.sim.use.1 <- parm.mat.sim.use
  
  ## convert parameters to capture probabilities
  prob.true.sim <- f_parm_prob(nstreams = nstreams,
                               parm.mat = parm.mat.sim.use.1)
  parm.setup$Sce[count] <- paste0("S", count)
  parm.setup$N[count] <- N
  parm.setup[count, parm.names] <- parm.mat.sim.use.1$value
  parm.setup$pc[count] <- (1 - prob.true.sim$prob[2^nstreams])
  
  dat.full.sim <- t(rmultinom(nsims, size = N,
                              prob = as.numeric(prob.true.sim$prob)))
  dat.full.sim <- as.data.frame(dat.full.sim)
  
  prof.full.mat <- expand.grid(rep(list(c(1,0)), nstreams))
  for(i in nstreams:1){
    prof.full.mat = prof.full.mat[order(prof.full.mat[,i],
                                        decreasing = T), ]
  }
  colnames(prof.full.mat) <- paste0("X", 1:nstreams)
  prof.char.vec <- paste0("n",
                          apply(prof.full.mat[-(2^nstreams),],
                                1, paste, collapse = ""))
  dat.sim <- dat.full.sim[,-(2^nstreams)]
  colnames(dat.sim) <- prof.char.vec
  dat.list[[count]] <- dat.sim
  count = count + 1
} # end loop over N
dat.sim.list <- list("dat" = dat.list,
                     "parm" = parm.setup)
save(dat.sim.list, file = paste0("./data/dat_sparse_M46.rda"))




############# SCENARIO IV ##############
# simulation results are summarized in Table S1
# SPARSE data were generated under the M_{T-2,T} for T = 3, 4
# 3-stream case: N = 100, 200, pc = 0.464
# 4-stream case: N = 100, 200, pc = 0.509
# Fitted Model:
# (a) 1st-order Markov model
# (b) M_{T-2,T} true model
# (c) AIC selected model, M_{T-2, T} and M*_{T-2,T}
#########################################
set.seed(1234)
# the number of trapping occasions, T = 3, 4 
nstreams.vec <- c(3, 4)
# the true population size N = 100, 200
N.true.vec <- c(100, 200)
# the total  number of simulations
nsims = 1000
for(nstreams in nstreams.vec){
  if(nstreams == 3){
    count = 1; dat.list <- list()
    
    psi.vec <- 0.1; psi = psi.vec[1]
    parm.names <- c("p1", "p21", "p21bar",
                    "p312", "p312bar", "p31bar2", "p31bar2bar")
    parm.setup <- as.data.frame(matrix(NA, nrow = length(psi.vec)*
                                         length(N.true.vec),
                                       ncol = 2^nstreams - 1 + 3))
    colnames(parm.setup) <- c("Sce", "N", parm.names, "pc")
    
    p1 = 0.3; p21 = 0.7; p21bar = 0.15; 
    p312 = 0.35; p31bar2 = 0.25
    # the M_{1, 3} model imposes only one untestable constraint: 
    # p31bar2bar = p312bar 
    p312bar = psi
    parm.vec <- c(p1, p21, p21bar, p312, p312bar, p31bar2, psi)
    # compute p_{hi} based on specified parameters
    prob.vec.true <- comp.prob(parm.vec)
    
    for(N in N.true.vec){
      
      parm.setup$Sce[count] <- paste0("S", count)
      parm.setup$N[count] <- N
      parm.setup[count, parm.names] <- parm.vec
      parm.setup$pc[count] <- (1 - prob.vec.true[8])
      
      # generate data based on the multinominal dist
      dat.full.sim <- t(rmultinom(nsims, size = N,
                                  prob = as.numeric(prob.vec.true)))
      dat.full.sim <- as.data.frame(dat.full.sim)
      
      prof.full.mat <- expand.grid(rep(list(c(1,0)), nstreams))
      for(i in nstreams:1){
        prof.full.mat = prof.full.mat[order(prof.full.mat[,i],
                                            decreasing = T), ]
      }
      colnames(prof.full.mat) <- paste0("X", 1:nstreams)
      prof.char.vec <- paste0("n",
                              apply(prof.full.mat[-(2^nstreams),],
                                    1, paste, collapse = ""))
      dat.sim <- dat.full.sim[,-(2^nstreams)]
      colnames(dat.sim) <- prof.char.vec
      dat.list[[count]] <- dat.sim
      count = count + 1
    } # end loop over N
    dat.sim.list <- list("dat" = dat.list,
                         "parm" = parm.setup)
    save(dat.sim.list, file = paste0("./data/dat_sparse_M13.rda"))
  }else if(nstreams == 4){
    
    ## get some estimated parameters from real data 
    y.sel = "1983" # 1983, 4 trapping occasions
    dat <- openCR::FebpossumCH[[y.sel]][,,1]
    dat.his <- apply(dat, 1, function(x) base::paste(x, collapse = ""))
    dat.his <- as.data.frame(table(dat.his))
    colnames(dat.his) <- c("prof", "count")
    hist.all <- get_allhist(nstreams)
    dat.fake <- data.frame("prof" = apply(hist.all, 1, paste, collapse = ""))
    dat.full <- left_join(dat.fake, dat.his, by = c("prof" = "prof"))
    dat.full$count[is.na(dat.full$count)] <- 0
    dat_aggre <- as.data.frame(matrix(dat.full$count, nrow = 1))
    colnames(dat_aggre) <- paste0("n", dat.full$prof)
    re.alter.i <- fit.markovian.assump(nstreams = nstreams, 
                                 dat_aggre = dat_aggre,
                                 assump = "2days", obscons = F, hessian = T)
    est.mat.i <- data.frame("parm_name_curr" = names(re.alter.i$parm),
                            "value" = as.numeric(re.alter.i$parm))
    parm.mat.sim <- left_join(re.alter.i$parm.name, est.mat.i,
                              by = c("parm_name_curr" = "parm_name_curr"))
    parm.mat.sim.use <- parm.mat.sim[,c("ID", "S_curr", "pre_hist",
                                        "parm_name", "value")]
    
    count = 1; dat.list <- list()
    psi.vec = c(0.2)
    psi = psi.vec[1]
    parm.setup <- 
      as.data.frame(matrix(NA, nrow = length(psi.vec)*length(N.true.vec),
                           ncol = 2^nstreams - 1 + 3))
    parm.names <- parm.mat.sim.use$parm_name
    colnames(parm.setup) <- c("Sce", "N", parm.names, "pc")
    probhist.setup <- 
      as.data.frame(matrix(NA, nrow = length(psi.vec)*length(N.true.vec),
                           ncol = 2^nstreams))
    colnames(probhist.setup) <- 
      paste0("p", apply(get_allhist_all(4), 1, paste, collapse = ""))
    
    for(N in N.true.vec){
      
      ## modify estimated parameters 
      parm.mat.sim.use.1 <- parm.mat.sim.use
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$value > 0.6] <-
        parm.mat.sim.use.1$value[parm.mat.sim.use.1$value > 0.6]*0.5
      parm.mat.sim.use.1$value <- round(parm.mat.sim.use.1$value, 2)
      
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p1"] <- 0.15
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p2_0"] <- 0.18
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p3_00"] <- 0.12
      
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p4_010"] <- 0.45
      
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p2_1"] <- 0.55
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p3_11"] <- 0.35
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p3_01"] <- 0.15
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p4_111"] <- 0.4
      
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p4_000"] <- psi
      parm.mat.sim.use.1$value[parm.mat.sim.use.1$parm_name == "p4_100"] <- psi
      
      ## convert parameters to capture probabilities
      prob.true.sim <- f_parm_prob(nstreams = nstreams,
                                   parm.mat = parm.mat.sim.use.1)
      parm.setup$Sce[count] <- paste0("S", count)
      parm.setup$N[count] <- N
      parm.setup[count, parm.names] <- parm.mat.sim.use.1$value
      parm.setup$pc[count] <- (1 - prob.true.sim$prob[2^nstreams])
      probhist.setup[count, ] <- as.numeric(prob.true.sim$prob)
      dat.full.sim <- t(rmultinom(nsims, size = N,
                                  prob = as.numeric(prob.true.sim$prob)))
      dat.full.sim <- as.data.frame(dat.full.sim)
      
      prof.full.mat <- expand.grid(rep(list(c(1,0)), nstreams))
      for(i in nstreams:1){
        prof.full.mat = prof.full.mat[order(prof.full.mat[,i],
                                            decreasing = T), ]
      }
      colnames(prof.full.mat) <- paste0("X", 1:nstreams)
      prof.char.vec <- paste0("n",
                              apply(prof.full.mat[-(2^nstreams),],
                                    1, paste, collapse = ""))
      dat.sim <- dat.full.sim[,-(2^nstreams)]
      colnames(dat.sim) <- prof.char.vec
      dat.list[[count]] <- dat.sim
      count = count + 1
    } # loop over N
    probhist.setup <- data.frame("Sce" = parm.setup$Sce,
                                 probhist.setup)
    dat.sim.list <- list("dat" = dat.list,
                         "parm" = parm.setup,
                         "probhist" = probhist.setup)
    save(dat.sim.list, file = paste0("./data/dat_sparse_M24.rda"))
  }
} # end loop over T 



