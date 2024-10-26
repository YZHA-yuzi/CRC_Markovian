############################################################################
# --------- Self-defined R functions used for the simualtion study ------- #
############################################################################

# A function to compute the plug-in estimator in Equation (11)
get_Nhat_psi <- function(nc, nlast, psi){
  return((nc - nlast) + nlast/psi)
}


# A function to compute prob of each capture profile based on parameters vector
# p1, p21, p21bar, p312, p312bar, p31bar2, psi T = 3
comp.prob <- function(x){
  p1 = x[1]; p21 = x[2]; p21bar = x[3];
  p312 = x[4]; p312bar = x[5]; p31bar2 = x[6]; psi = x[7]
  
  p111 = p1*p21*p312; p110 = p1*p21*(1-p312)
  p101 = p1*(1-p21)*p312bar; p100 = p1*(1-p21)*(1-p312bar)
  p011 = (1-p1)*p21bar*p31bar2; p010 = (1-p1)*p21bar*(1-p31bar2)
  p001 = (1-p1)*(1-p21bar)*psi; p000 = (1-p1)*(1-p21bar)*(1-psi)
  prob.vec = c(p111, p110, p101, p100, p011, p010, p001, p000)
  return(prob.vec)
}


# A function to generate all possible OBSERVED capture history for a given T
get_allhist <- function(nstreams){
  df <- expand.grid(rep(list(c(1,0)), nstreams))[-(2^nstreams), ]
  for(i in nstreams:1){
    df = df[order(df[,i], decreasing = T), ]
  }
  return(df)
}


# A function to generate all possible capture history for a given T
get_allhist_all <- function(nstreams){
  df <- expand.grid(rep(list(c(1,0)), nstreams))
  for(i in nstreams:1){
    df = df[order(df[,i], decreasing = T), ]
  }
  return(df)
}


# A function from parameters to capture probabilities
f_parm_prob <- function(nstreams, parm.mat){
  # parm.mat:
  # a data framework obtained from f_prob_parm 
  parm.dic.use <- parm.mat
  hist.use <- get_allhist_all(nstreams = nstreams)
  parmid.forhist <- matrix(NA, ncol = nstreams, nrow = 2^nstreams)
  parmname.forhist <- matrix(NA, ncol = nstreams, nrow = 2^nstreams)
  sign.forhist <- matrix(NA, ncol = nstreams, nrow = 2^nstreams)
  for(j in 1:nstreams){
    if(j == 1){
      index <- hist.use[,j] == 1
      sign.forhist[index, j] <- "+"
      sign.forhist[!index, j] <- "1-"
      parmid.forhist[,j] <- 1 
      parmname.forhist[,j] <- "p1"
    }else if(j == 2){
      hist.pre.j <- as.character(hist.use[,j-1])
      parm.j <- paste0(paste0("p", j, "_"), hist.pre.j)
      index <- hist.use[,j] == 1
      sign.forhist[index, j] <- "+"
      sign.forhist[!index, j] <- "1-"
      parmid.forhist[,j] <- parm.dic.use$ID[
        sapply(parm.j, function(x) which(parm.dic.use$parm_name == x))]
      parmname.forhist[,j] <- parm.j
    }else{
      hist.pre.j <- apply(hist.use[,1:(j-1)], 1, paste, collapse = "")
      parm.j <- paste0(paste0("p", j, "_"), hist.pre.j)
      index <- hist.use[,j] == 1
      sign.forhist[index, j] <- "+"
      sign.forhist[!index, j] <- "1-"
      parmid.forhist[,j] <- parm.dic.use$ID[
        sapply(parm.j, function(x) which(parm.dic.use$parm_name == x))]
      parmname.forhist[,j] <- parm.j
    }
  }
  
  ### get expression of capture probabilities of different capture histories
  expr.vec <- c()
  for(i in 1:(2^nstreams)){
    expr0 <- paste0(sign.forhist[i, ], "parms[", parmid.forhist[i, ], "]")
    expr0 <- sapply(expr0, function(x) paste0("(", x, ")"))
    expr0 <- paste(expr0, collapse = "*")
    expr0 <- paste0("max(", expr0, ", 1e-10)")
    expr.vec <- c(expr.vec, expr0)
  }
  expr.vec.use <- expr.vec
  parms <- parm.dic.use$value
  
  prob.re <- sapply(1:length(parms), function(l) 
    eval(parse(text = expr.vec[l])))
  
  prob.mat <- data.frame(hist.use, "prob" = c(prob.re, 1-sum(prob.re)))
  colnames(prob.mat) <- c(paste0("S", 1:nstreams), "prob")
  return(prob.mat)
}

#### A function to fit proposed model while 
#### imposing Markovian-type constraints
## e.g., psi = p612bar3bar4bar5bar = p61bar23bar4bar5bar = p6123bar4bar5bar 
fit.markovian.assump <- function(nstreams, dat_aggre, assump, 
                                 obscons = FALSE, hessian = FALSE){
  ## INPUT: 
  ## nstreams: the number of data streams
  ## dat_aggre: aggregate data sorted by 111, 110, ..., 001
  ## assump: "4days", "3days", "2days", "1day"
  ## hessian = T, don't compute hessian matrix = no SE
  ## obscons = T, also impose other testable constraints 
  parm.dic <- NULL
  parm.dic <- rbind(parm.dic, c("p1", NA))
  for(i in 1:(nstreams-1)){
    name.i = paste0("p", i+1)
    if(i == 1){
      hist.i <- as.character(get_allhist_all(i))
      parm.dic <- rbind(parm.dic, cbind(rep(name.i, 2^i), hist.i))
    }else{
      hist.i <- apply(get_allhist_all(i), 1, paste, collapse = "")
      parm.dic <- rbind(parm.dic, cbind(rep(name.i, 2^i), hist.i))
    }
  }
  parm.dic <- as.data.frame(parm.dic)
  colnames(parm.dic) <- c("S_curr", "pre_hist")
  rownames(parm.dic) <- NULL
  parm.dic <- data.frame("ID" = 1:nrow(parm.dic), parm.dic)
  parm.dic[["parm_name"]] <- 
    c("p1", apply(parm.dic[-1, ], 1, function(x) paste0(x[2], "_", x[3])))
  
  parm.dic$parm_name_curr <- parm.dic$parm_name
  naprt = as.numeric(substr(assump, 1, 1))
  if(obscons){
    con.s <- (nstreams:1)[(nstreams:1) - naprt > 1]
    for(ss in 1:length(con.s)){
      con.s.curr <- con.s[ss]
      if((con.s.curr - (naprt + 1)) == 1){
        pre.hist.curr <- get_allhist_all(con.s.curr - (naprt + 1))
        pre.hist.curr.1 <- as.character(pre.hist.curr)
      }else{
        pre.hist.curr <- get_allhist_all(con.s.curr - (naprt + 1))
        pre.hist.curr.1 <- apply(pre.hist.curr, 1, paste, collapse = "")
      }
      ## other testable constraints ##
      if(naprt == 1){ mid.hist.curr <- get_allhist_all(naprt) }else{
        mid.hist.curr <- apply(get_allhist_all(naprt), 1, paste, collapse = "")
      }
      for(smid in 1:(length(mid.hist.curr)-1)){
        parm.dic$parm_name_curr[parm.dic$S_curr == paste0("p", con.s.curr) & 
                                  parm.dic$pre_hist %in% 
                                  paste0(pre.hist.curr.1,
                                         paste0(mid.hist.curr[smid], collapse = ""))] <- 
          paste0("theta_", ss, "_", smid)
      }
      parm.dic$parm_name_curr[parm.dic$S_curr == paste0("p", con.s.curr) & 
                                parm.dic$pre_hist %in% 
                                paste0(pre.hist.curr.1,
                                       paste0(rep("0", naprt), collapse = ""))] <- 
        paste0("psi_", ss)
    }
  }else{
    con.s <- (nstreams:1)[(nstreams:1) - naprt > 1]
    for(ss in 1:length(con.s)){
      con.s.curr <- con.s[ss]
      if((con.s.curr - (naprt + 1)) == 1){
        pre.hist.curr <- get_allhist_all(con.s.curr - (naprt + 1))
        pre.hist.curr.1 <- as.character(pre.hist.curr)
      }else{
        pre.hist.curr <- get_allhist_all(con.s.curr - (naprt + 1))
        pre.hist.curr.1 <- apply(pre.hist.curr, 1, paste, collapse = "")
      }
      parm.dic$parm_name_curr[parm.dic$S_curr == paste0("p", con.s.curr) & 
                                parm.dic$pre_hist %in% 
                                paste0(pre.hist.curr.1,
                                       paste0(rep("0", naprt), collapse = ""))] <- 
        paste0("psi_", ss)
    }
  }
  
  parm.vec.curr <- unique(parm.dic$parm_name_curr)
  id.vec <- lapply(parm.dic$parm_name_curr, 
                   function(x) which(parm.vec.curr == x))
  parm.dic$ID.curr <- do.call(c, id.vec)
  
  parm.dic.use <- parm.dic
  hist.use <- get_allhist_all(nstreams = nstreams)
  parmid.forhist <- matrix(NA, ncol = nstreams, nrow = 2^nstreams)
  parmname.forhist <- matrix(NA, ncol = nstreams, nrow = 2^nstreams)
  sign.forhist <- matrix(NA, ncol = nstreams, nrow = 2^nstreams)
  for(j in 1:nstreams){
    if(j == 1){
      index <- hist.use[,j] == 1
      sign.forhist[index, j] <- "+"
      sign.forhist[!index, j] <- "1-"
      parmid.forhist[,j] <- 1 
      parmname.forhist[,j] <- "p1"
    }else if(j == 2){
      hist.pre.j <- as.character(hist.use[,j-1])
      parm.j <- paste0(paste0("p", j, "_"), hist.pre.j)
      index <- hist.use[,j] == 1
      sign.forhist[index, j] <- "+"
      sign.forhist[!index, j] <- "1-"
      parmid.forhist[,j] <- parm.dic.use$ID.curr[
        sapply(parm.j, function(x) which(parm.dic.use$parm_name == x))]
      parmname.forhist[,j] <- parm.j
    }else{
      hist.pre.j <- apply(hist.use[,1:(j-1)], 1, paste, collapse = "")
      parm.j <- paste0(paste0("p", j, "_"), hist.pre.j)
      index <- hist.use[,j] == 1
      sign.forhist[index, j] <- "+"
      sign.forhist[!index, j] <- "1-"
      parmid.forhist[,j] <- parm.dic.use$ID.curr[
        sapply(parm.j, function(x) which(parm.dic.use$parm_name == x))]
      parmname.forhist[,j] <- parm.j
    }
  }
  
  ### get expression of capture probabilities of different capture histories
  expr.vec <- c()
  for(i in 1:(2^nstreams)){
    expr0 <- paste0(sign.forhist[i, ], "parms[", parmid.forhist[i, ], "]")
    expr0 <- sapply(expr0, function(x) paste0("(", x, ")"))
    expr0 <- paste(expr0, collapse = "*")
    expr0 <- paste0("max(", expr0, ", 1e-10)")
    expr.vec <- c(expr.vec, expr0)
  }
  expr.vec.use <- expr.vec
  parm.names <- parm.vec.curr
  
  ll_nstreams <- function(parms){
    ll_vec <- c()
    for(l in 1:(length(expr.vec)-1)){
      lamb.l <- eval(parse(text = expr.vec[l]))
      ll_vec <- c(ll_vec, df_aggre[1, l]*log(lamb.l))
    }
    nc.forll <- sum(df_aggre)
    ll = sum(ll_vec) - 
      nc.forll*log(1 - eval(parse(text = expr.vec[length(expr.vec)])))
    return(-1*ll)
  }
  parnames(ll_nstreams) <- parm.names
  
  nparms = length(parm.names)
  initial_parm = rep(0.1, nparms)
  lower_parm = rep(1e-6, nparms)
  upper_parm = rep(1 - 1e-6, nparms)
  names(initial_parm) <- names(lower_parm) <- names(upper_parm) <- 
    parm.names
  
  fit_alter <- mle2(ll_nstreams, start = initial_parm, 
                    data = list(df_aggre = dat_aggre,
                                expr.vec = expr.vec.use),
                    method="L-BFGS-B", lower= lower_parm, upper = upper_parm,
                    skip.hessian = hessian)
  
  parms <- coef(fit_alter)
  pc.alter <- 1 - eval(parse(text = expr.vec[2^nstreams])) 
  nc.dat <- sum(dat_aggre)
  Nhat <- nc.dat/pc.alter
  parm.hat <- coef(fit_alter)
  return(list("Nhat" = Nhat, "pc" = pc.alter,
              "parm" = parm.hat, "parm.name" = parm.dic.use,
              "AIC" = AIC(fit_alter), "logl" = logLik(fit_alter),
              "num_parm" = length(parm.hat)))
}



get.status <- function(x, assump = "MM2b"){
  # x is a vector representing a capture history
  if(assump == "MM2b"){
    re <- ifelse(x[1] == 1, "c", "a")
    for(i in 2:length(x)){
      if(x[i] == 1){
        re.i <- "c"
      }else if(x[i] == 0 & sum(x[1:(i-1)]) != 0){
        re.i <- "b"
      }else if(x[i] == 0 & sum(x[1:(i-1)]) == 0){
        re.i <- "a"
      }
      re <- c(re, re.i)
    }
    return(c("a", re))
  }else if(assump == "MM1b"){
    return(c("0", as.character(as.numeric(x))))
  }
}


# A function to construct data used for Yang's method
get.dat.yang <- function(dat, assump = "MM2b", nstreams){
  # 1. get all possible capture histories stored in a matrix
  prof.full.mat <- expand.grid(rep(list(c(1,0)), nstreams))
  for(i in nstreams:1){
    prof.full.mat = prof.full.mat[order(prof.full.mat[,i], 
                                        decreasing = T), ]
  }
  colnames(prof.full.mat) <- paste0("X", 1:nstreams)
  prof.obs <- prof.full.mat[-(2^nstreams),]
  # 2. get status at different time points using notations in Yang 
  status.yang <- as.data.frame(matrix(NA, ncol = nstreams + 1, 
                                      nrow = 2^nstreams - 1))
  colnames(status.yang) <- paste0("t", c(0:nstreams))
  status.yang.1 <- status.yang[,-1]
  for(i in 1:nrow(prof.obs)){
    s.i <- get.status(x = prof.obs[i,], assump = assump)
    status.yang[i, ] <- s.i
    for(l in 2:(nstreams+1)){
      status.yang.1[i, l-1] <- paste(s.i[(l-1):l], collapse = "")
    }
  }
  status.yang$count <- as.numeric(dat)
  status.yang.1$count <- as.numeric(dat)
  if(assump == "MM2b"){
    status.name <- c("aa", "ac", "bb", "bc", "cb", "cc")
  }else if(assump == "MM1b"){
    status.name <- c("00", "01", "10", "11")
  }
  
  # 3. get counts for different status
  dat.yang <- as.data.frame(matrix(NA, ncol = length(status.name), 
                                   nrow = nstreams))
  colnames(dat.yang) <- status.name
  for(s in status.name){
    for(i in 1:nstreams){
      sel.index = which(status.yang.1[[paste0("t", i)]] == s)
      dat.yang[[s]][i] <- sum(status.yang.1$count[sel.index])
    }
  }
  dat.yang.final <- as.data.frame(matrix(colSums(dat.yang), nrow = 1))
  colnames(dat.yang.final) <- status.name
  return(dat.yang.final)
}

# A function to get point estimator proposed in Yang 
get.Nhat.yang <- function(dat.yang, nstreams, nc, assump = "MM2b"){
  if(assump == "MM2b"){
    y.dat <- c(nc, nstreams, sum(dat.yang[["ac"]]), 
               sum(dat.yang[,c("aa", "ac")]))
    names(y.dat) <- c("M", "T", "n_ac", "n_adot")
    fN.yang <- function(x, y){
      M = y[["M"]]; num.t = y[["T"]]; 
      n_ac = y[["n_ac"]]; n_adot = y[["n_adot"]]
      re = (1 - M/x) - (1 - n_ac/(num.t*(x-M) + n_adot))^num.t
      return(re)
    }
    pos.range <- c(y.dat[["M"]], 10*y.dat[["M"]])
    check.use <- sapply(pos.range, fN.yang, y = y.dat)
    if(all(check.use < 0)){
      Nhat.yang <- NA
    }else{
      Nhat.yang <- uniroot(fN.yang, y = y.dat, 
                           interval = c(y.dat[["M"]], 10*y.dat[["M"]]))$root
    }
    return(Nhat.yang)
  }else if(assump == "MM1b"){
    y.dat <- c(nc, nstreams, sum(dat.yang[["01"]]), 
               sum(dat.yang[,c("01", "00")]))
    names(y.dat) <- c("M", "T", "n_01", "n_0dot")
    fN.yang <- function(x, y){
      M = y[["M"]]; num.t = y[["T"]]; 
      n_01 = y[["n_01"]]; n_0dot = y[["n_0dot"]]
      re = (1 - M/x) - (1 - n_01/(num.t*(x-M) + n_0dot))^num.t
      return(re)
    }
    pos.range <- c(y.dat[["M"]], 10*y.dat[["M"]])
    check.use <- sapply(pos.range, fN.yang, y = y.dat)
    if(all(check.use < 0)){
      Nhat.yang <- NA
    }else{
      Nhat.yang <- uniroot(fN.yang, y = y.dat, 
                           interval = c(y.dat[["M"]], 10*y.dat[["M"]]))$root
    }
    return(Nhat.yang)
  }
}

### A function to get estimated parameters in Yang's model 
get.parm.yang <- function(Nhat, dat.yang, 
                          nstreams, nc, assump = "MM2b"){
  if(assump == "MM2b"){
    n_ac = sum(dat.yang[["ac"]])
    n_adot = sum(dat.yang[,c("aa", "ac")])
    n_bc = sum(dat.yang[["bc"]])
    n_bdot = sum(dat.yang[,c("bb", "bc")])
    n_cc = sum(dat.yang[["cc"]])
    n_cdot = sum(dat.yang[,c("cb", "cc")])
    p_cc = n_cc/n_cdot
    p_bc = n_bc/n_bdot
    p_ac = n_ac/(nstreams*(Nhat - nc) + n_adot)
    re = data.frame("p_ac" = p_ac, "p_bc" = p_bc, "p_cc" = p_cc)
    return(re)
  }else if(assump == "MM1b"){
    n_01 = dat.yang[["01"]]
    n_11 = dat.yang[["11"]]
    n_0dot = sum(dat.yang[,c("00","01")])
    n_1dot = sum(dat.yang[,c("10","11")])
    p_01 = n_01/(nstreams*(Nhat - nc) + n_0dot)
    p_11 = n_11/n_1dot
    # p_ac = p_bc = p_01, and p_cc = p_11
    re = data.frame("p_01" = p_01, "p_11" = p_11)
    return(re)
  }
}

### A function to get variance of Nhat obtained using Yang's model
get.var.yang <- function(Nhat, parmhat, nstreams, 
                         assump = "MM2b"){
  # assump: MM2b (p_ac, p_bc, p_cc)
  # MM1b (p_ac = p_bc, p_cc)
  if(assump == "MM2b"){
    Q <- (1 - parmhat[["p_ac"]])^nstreams
    p_aa <- 1 - parmhat[["p_ac"]]
    re.var <- Nhat/(-1 + 1/Q + 
                      ((nstreams^2)*(1 - p_aa)/(1 - Q))*(1 - 1/p_aa))
    return(re.var)
  }else if(assump == "MM1b"){
    ## eqn. (10) in Yang's paper
    p_01 <- parmhat[["p_01"]]
    p_00 <- 1 - p_01
    p_11 <- parmhat[["p_11"]]
    p_10 <- 1 - p_11
    comp.mid <- rep(NA, nstreams)
    for(l in 1:nstreams){
      comp.mid[l] = (p_10 + (p_11 - p_01)^(l-1)*p_01)/(p_01 + p_10)
    }
    Q <- (1 - p_01)^nstreams
    re.var <- Nhat/(-1 + 1/Q + 
                      ((nstreams^2)/sum(comp.mid))*(1 - 1/p_00))
    return(re.var)
  }
}

### A function to CI of Nhat
### normal approximation, log-transformation
get.CI.yang <- function(Nhat, varhat){
  C = exp(1.96*sqrt(log(1 + varhat/Nhat^2)))
  return(c(Nhat/C, Nhat*C))
}

#### A function to fit Yang's model using the alternative framework
fit.alter.yang <- function(nstreams, dat_aggre, assump){
  ## INPUT: 
  ## nstreams: the number of data streams
  ## dat_aggre: aggregate data sorted by 111, 110, ..., 001
  ## assump: "MM1b" or "MM2b"
  parm.dic <- NULL
  parm.dic <- rbind(parm.dic, c("p1", NA))
  for(i in 1:(nstreams-1)){
    name.i = paste0("p", i+1)
    if(i == 1){
      hist.i <- as.character(get_allhist_all(i))
      parm.dic <- rbind(parm.dic, cbind(rep(name.i, 2^i), hist.i))
    }else{
      hist.i <- apply(get_allhist_all(i), 1, paste, collapse = "")
      parm.dic <- rbind(parm.dic, cbind(rep(name.i, 2^i), hist.i))
    }
  }
  parm.dic <- as.data.frame(parm.dic)
  colnames(parm.dic) <- c("S_curr", "pre_hist")
  rownames(parm.dic) <- NULL
  parm.dic <- data.frame("ID" = 1:nrow(parm.dic), parm.dic)
  parm.dic[["parm_name"]] <- 
    c("p1", apply(parm.dic[-1, ], 1, function(x) paste0(x[2], "_", x[3])))
  
  parm.type.1 <- c("p1", paste0("p", 2:nstreams ,"_",
                                sapply(1:(nstreams-1), function(x) 
                                  paste(rep(0, x), collapse = ""))))
  parm.type.2 <- c("p2_1")
  parm.type.3 <- c()
  for(j in 3:nstreams){
    k = (j-2)
    if(k == 1){
      hist.k <- as.character(get_allhist_all(k))
      index.k <- ! sapply(strsplit(hist.k, ""), function(x) all(x == "0"))
      parm.type.2 <- c(parm.type.2, paste0(paste0("p", j, "_"), hist.k, "1"))
      parm.type.3 <- c(parm.type.3, 
                       paste0(paste0("p", j, "_"), hist.k[index.k], "0"))
    }else{
      hist.k <- apply(get_allhist_all(k), 1, paste, collapse = "")
      index.k <- ! sapply(strsplit(hist.k, ""), function(x) all(x == "0"))
      parm.type.2 <- c(parm.type.2, paste0(paste0("p", j, "_"), hist.k, "1"))
      parm.type.3 <- c(parm.type.3, 
                       paste0(paste0("p", j, "_"), hist.k[index.k], "0"))
    }
  }
  
  if(assump == "MM1b"){
    parm.dic$ID[parm.dic$parm_name %in% parm.type.1] <- 1
    parm.dic$ID[parm.dic$parm_name %in% parm.type.2] <- 2
    parm.dic$ID[parm.dic$parm_name %in% parm.type.3] <- 1
  }else if(assump == "MM2b"){
    parm.dic$ID[parm.dic$parm_name %in% parm.type.1] <- 1
    parm.dic$ID[parm.dic$parm_name %in% parm.type.2] <- 2
    parm.dic$ID[parm.dic$parm_name %in% parm.type.3] <- 3
  }
  
  parm.dic.use <- parm.dic
  hist.use <- get_allhist_all(nstreams = nstreams)
  parmid.forhist <- matrix(NA, ncol = nstreams, nrow = 2^nstreams)
  parmname.forhist <- matrix(NA, ncol = nstreams, nrow = 2^nstreams)
  sign.forhist <- matrix(NA, ncol = nstreams, nrow = 2^nstreams)
  for(j in 1:nstreams){
    if(j == 1){
      index <- hist.use[,j] == 1
      sign.forhist[index, j] <- "+"
      sign.forhist[!index, j] <- "1-"
      parmid.forhist[,j] <- 1 
      parmname.forhist[,j] <- "p1"
    }else if(j == 2){
      hist.pre.j <- as.character(hist.use[,j-1])
      parm.j <- paste0(paste0("p", j, "_"), hist.pre.j)
      index <- hist.use[,j] == 1
      sign.forhist[index, j] <- "+"
      sign.forhist[!index, j] <- "1-"
      parmid.forhist[,j] <- parm.dic.use$ID[
        sapply(parm.j, function(x) which(parm.dic.use$parm_name == x))]
      parmname.forhist[,j] <- parm.j
    }else{
      hist.pre.j <- apply(hist.use[,1:(j-1)], 1, paste, collapse = "")
      parm.j <- paste0(paste0("p", j, "_"), hist.pre.j)
      index <- hist.use[,j] == 1
      sign.forhist[index, j] <- "+"
      sign.forhist[!index, j] <- "1-"
      parmid.forhist[,j] <- parm.dic.use$ID[
        sapply(parm.j, function(x) which(parm.dic.use$parm_name == x))]
      parmname.forhist[,j] <- parm.j
    }
  }
  
  ### get expression of capture probabilities of different capture histories
  expr.vec <- c()
  for(i in 1:(2^nstreams)){
    expr0 <- paste0(sign.forhist[i, ], "parms[", parmid.forhist[i, ], "]")
    expr0 <- sapply(expr0, function(x) paste0("(", x, ")"))
    expr0 <- paste(expr0, collapse = "*")
    expr0 <- paste0("max(", expr0, ", 1e-10)")
    expr.vec <- c(expr.vec, expr0)
  }
  expr.vec.use <- expr.vec
  
  if(assump == "MM1b"){
    parm.names <- c("p1", "p21")
  }else if(assump == "MM2b"){
    parm.names <- c("p1", "p21", "p312bar")
  }
  ll_nstreams <- function(parms){
    ll_vec <- c()
    for(l in 1:(length(expr.vec)-1)){
      lamb.l <- eval(parse(text = expr.vec[l]))
      ll_vec <- c(ll_vec, df_aggre[1, l]*log(lamb.l))
    }
    nc.forll <- sum(df_aggre)
    ll = sum(ll_vec) - 
      nc.forll*log(1 - eval(parse(text = expr.vec[length(expr.vec)])))
    return(-1*ll)
  }
  parnames(ll_nstreams) <- parm.names
  
  nparms = length(parm.names)
  initial_parm = rep(0.1, nparms)
  lower_parm = rep(1e-6, nparms)
  upper_parm = rep(1 - 1e-6, nparms)
  names(initial_parm) <- names(lower_parm) <- names(upper_parm) <- 
    parm.names
  
  fit_alter <- mle2(ll_nstreams, start = initial_parm, 
                    data = list(df_aggre = dat_aggre,
                                expr.vec = expr.vec.use),
                    method="L-BFGS-B", lower= lower_parm, upper = upper_parm,
                    skip.hessian = T)
  pc.alter <- (1 - (1 - coef(fit_alter)[["p1"]])^nstreams)
  nc.dat <- sum(dat_aggre)
  Nhat <- nc.dat/pc.alter
  parm.hat <- coef(fit_alter)
  return(list("Nhat" = Nhat, "parm" = parm.hat,
              "AIC" = AIC(fit_alter)))
}




### get dirichlet-based CI based on fitted cell counts
Dir_CI_alter_fitted <- function(df_aggre, n.post, 
                                assump = "1day", nstreams = 6,
                                obscons = FALSE, hessian = FALSE,
                                Dir_prior = "equprob", 
                                prob.cond.input = NULL,
                                a = 0, b = 0,
                                hist.sel.num, hist.sel.den){
  ## INPUT: 
  ## Dir_prior: type of Dirichlet prior
  ## "optimal": weigthed average Fienberg (1973)
  ## "improper": all 0s
  ## "Jeffreys": 0.5s
  ## "equprob": all 1/(2^K-1)s
  ## a = 1, b = 0
  ## a = 0.5, b = 0.5, Jeffreys'
  ## a = 0, b = 0, no bias-correction
  
  loc1 = which(colnames(df_aggre) == paste0("n", hist.sel.num))
  loc2 = which(colnames(df_aggre) == paste0("n", hist.sel.den))
  
  if(Dir_prior == "improper"){
    d_p = 0
    d_post = as.numeric(df_aggre) + d_p
  }else if(Dir_prior == "Jeffreys"){
    d_p = 0.5
    d_post = as.numeric(df_aggre) + d_p
  }else if(Dir_prior == "equprob"){
    d_p = 1/(2^nstreams-1)
    d_post = as.numeric(df_aggre) + d_p
  }
  
  # generate simulated data by assigning Dirichlet priors
  pstarcond <- rdirichlet(n.post, alpha = d_post) # p111, ..., p001
  ## define a function to find pc
  f <- function(nc, parm, assump, obscons, hessian, nstreams, a, b){
    dat_simed <- as.data.frame(t(nc*parm))
    re <- fit.markovian.assump(nstreams = nstreams, dat_aggre = dat_simed,
                         assump = assump, obscons = obscons,
                         hessian = hessian)
    
    ## GET fitted cell counts ###
    est.mat.i <- data.frame("parm_name_curr" = names(re$parm),
                            "value" = as.numeric(re$parm))
    parm.mat.sim <- left_join(re$parm.name, est.mat.i,
                              by = c("parm_name_curr" = "parm_name_curr"))
    ## 3. convert parameters to capture probabilities
    parm.mat.sim.use <- parm.mat.sim[,c("ID", "S_curr", "pre_hist",
                                        "parm_name", "value")]
    prob.true.sim <- f_parm_prob(nstreams = nstreams,
                                 parm.mat = parm.mat.sim.use)
    nfitted.i <- re$Nhat*prob.true.sim$prob
    
    ### computed Nhat based on fitted cell counts ###
    n.num.i <- nfitted.i[loc1] # "n100001"
    n.den.i <- nfitted.i[loc2] # "n100000"
    
    psi.i <- (n.num.i + a)/(a + b + n.num.i + n.den.i)
    nc.fitted <- sum(nfitted.i[-c(2^nstreams)])
    nlast.fitted <- nfitted.i[(2^nstreams-1)]
    Nhat.sim <- get_Nhat_psi(nc = nc.fitted, 
                             nlast = nlast.fitted, psi = psi.i)
    pcj <- nc.fitted/Nhat.sim
    return(pcj)
  }
  # start = proc.time()[3]
  pcj <- apply(pstarcond, 1, f,
               nc = as.numeric(sum(df_aggre)),
               assump = assump, obscons = obscons, hessian = hessian,
               nstreams = nstreams,
               a = a, b = b)
  Nnew <- round(rowSums(df_aggre)/pcj)
  # not condition on nc
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) rbinom(1, x[1], x[2]))
  Npost <- ncnew/pcj
  nc <- as.numeric(sum(df_aggre))
  Npost[Npost < nc] <- nc
  re <- list(post = Npost, median = median(Npost, na.rm = T),
             mean = mean(Npost, na.rm = T),
             CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T)),
             "pc" = pcj)
  return(re)
}



### get dirichlet-based CI based on observed cell counts
### under the F4 model
Dir_CI_alter_Fmax <- function(df_aggre, n.post, 
                              nstreams = 6, Dir_prior = "equprob", 
                              a = 0, b = 0,
                              hist.sel.num, hist.sel.den){
  ## INPUT: 
  ## Dir_prior: type of Dirichlet prior
  ## "optimal": weigthed average Fienberg (1973)
  ## "improper": all 0s
  ## "Jeffreys": 0.5s
  ## "equprob": all 1/(2^K-1)s
  ## a = 1, b = 0
  ## a = 0.5, b = 0.5, Jeffreys'
  ## a = 0, b = 0, no bias-correction
  
  loc1 = which(colnames(df_aggre) == paste0("n", hist.sel.num))
  loc2 = which(colnames(df_aggre) == paste0("n", hist.sel.den))
  
  if(Dir_prior == "improper"){
    d_p = 0
    d_post = as.numeric(df_aggre) + d_p
    # cat(Dir_prior)
  }else if(Dir_prior == "Jeffreys"){
    d_p = 0.5
    d_post = as.numeric(df_aggre) + d_p
  }else if(Dir_prior == "equprob"){
    d_p = 1/(2^nstreams-1)
    d_post = as.numeric(df_aggre) + d_p
  }
  
  # generate simulated data by assigning Dirichlet priors
  pstarcond <- rdirichlet(n.post, alpha = d_post) # p111, ..., p001
  ## define a function to find pc
  f <- function(nc, parm, nstreams, a, b){
    dat_simed <- as.data.frame(t(nc*parm))
    ### computed Nhat based on fitted cell counts ###
    n.num.i <- dat_simed[loc1] # "n100001"
    n.den.i <- dat_simed[loc2] # "n100000"
    
    psi.i <- (n.num.i + a)/(a + b + n.num.i + n.den.i)
    nc.fitted <- sum(dat_simed[-c(2^nstreams)])
    nlast.fitted <- dat_simed[(2^nstreams-1)]
    Nhat.sim <- get_Nhat_psi(nc = nc.fitted, 
                             nlast = nlast.fitted, psi = psi.i)
    pcj <- as.numeric(nc.fitted/Nhat.sim)
    return(pcj)
  }
  # start = proc.time()[3]
  pcj <- apply(pstarcond, 1, f, nc = as.numeric(sum(df_aggre)),
               nstreams = nstreams, a = a, b = b)
  Nnew <- round(rowSums(df_aggre)/pcj)
  # not condition on nc
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) rbinom(1, x[1], x[2]))
  Npost <- ncnew/pcj
  nc <- as.numeric(sum(df_aggre))
  Npost[Npost < nc] <- nc
  re <- list(post = Npost, median = median(Npost, na.rm = T),
             mean = mean(Npost, na.rm = T),
             CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T)),
             "pc" = pcj)
  return(re)
}


# A function to summarize results ######
get.sum <- function(re.i, model.name, est.name){
  keep.index = is.finite(re.i$Nhat)
  Ntrue = re.i$Ntrue[1]
  sumre.l <-
    data.frame("Model" = model.name,
               "Estimator" = est.name,
               "true" = Ntrue,
               "est" = mean(re.i$Nhat[keep.index]),
               "bias" = mean(re.i$Nhat[keep.index]) - Ntrue,
               "rel.bias" = (mean(re.i$Nhat[keep.index]) - Ntrue)/Ntrue,
               "emp.sd" = sd(re.i$Nhat[keep.index]),
               "cover" = mean(re.i$lci[keep.index] <= Ntrue &
                                re.i$uci[keep.index] >= Ntrue),
               "width" = mean((re.i$uci - re.i$lci)[keep.index]),
               "med.width" = median((re.i$uci - re.i$lci)[keep.index]),
               "mis.low" = mean(re.i$uci[keep.index] <= Ntrue),
               "mis.high" = mean(re.i$lci[keep.index] >= Ntrue),
               "RMSE" = sqrt(mean((re.i$Nhat[keep.index] - Ntrue)^2)),
               "num.keep" = sum(keep.index))
  return(sumre.l)
}

# A function to compute the number of constraints defined by a specified model
num.constraints <- function(nstreams, j.assumed, testable = F){
  re.num.1 = sum(sapply(1:(nstreams - j.assumed - 1), function(x) 2^x-1))
  if(!testable){
    # compute the number of constraints imposed in the M_{j,T} model
    return(re.num.1)
  }else{
    # compute the number of constraints imposed in the M*_{j,T} model
    re.num.2 = (2^j.assumed - 1)*sum(sapply(1:(nstreams - j.assumed - 1), 
                                            function(x) 2^x-1))
    return(re.num.2 + re.num.1)
  }
}

