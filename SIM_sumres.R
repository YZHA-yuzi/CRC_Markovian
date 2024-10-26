############################################################################
# -------------------  R codes for the simulation study ------------------ #
# NOTES:
# summarize simulation results 
# please make sure the working directory is the folder where 
# files "FUNs.R" is saved.
############################################################################

library(writexl)

## source self-defined R functions
source("FUNs.R")

## create a folder to store tables and figures to present
## simulation and real data analysis results
if(!base::dir.exists(paths = "./TabFigs")){
  base::dir.create("TabFigs")
}

tab.sum.excel <- replicate(4, list(), simplify = F)
names(tab.sum.excel) <- c(paste0("Tab", 1:3), "TabS1")

#########################################
# generate Tables 1, 2, and S1 
#########################################
tab.name.vec <- c(paste0("Tab", 1:2), "TabS1")
for(i.tab in 1:length(tab.name.vec)){
  tab.name <- tab.name.vec[i.tab]
  tab.sum.all <- list()
  count1 = 1
  for(num.s in c(3, 4)){
    ## read in simulated datasets
    if(tab.name == "Tab1"){
      load(paste0("./data/dat_M", num.s - 2, num.s, ".rda")) # dat.sim.list
      gen.m.s = paste0("M", num.s - 2, num.s)
    }else if(tab.name == "Tab2"){
      load(paste0("./data/dat_T", num.s, "_MarkovModel.rda")) # dat.sim.list
      gen.m.s = paste0("1st-order Markov model")
    }else if(tab.name == "TabS1"){
      load(paste0("./data/dat_sparse_M", num.s - 2, num.s, ".rda")) # dat.sim.list
      gen.m.s = paste0("M", num.s - 2, num.s)
    }
    N.true.vec <- dat.sim.list$parm$N
    count = 1
    for(N.true in N.true.vec){
      re.yang.all <- re.psi.all <- re.aic.all <- list()
      for(B_index in 1:10){
        if(tab.name == "Tab1"){
          load(paste0("./inter_res/Res_M", num.s-2, num.s, "_N", N.true,
                      "_B", B_index, ".rda"))
        }else if(tab.name == "Tab2"){
          load(paste0("./inter_res/Res_T", num.s, "_Markov_N", N.true,
                      "_B", B_index, ".rda"))
        }else if(tab.name == "TabS1"){
          load(paste0("./inter_res/Res_sparse_M", num.s-2, num.s, 
                      "_N", N.true, "_B", B_index, ".rda"))
        }
        re.yang.all[[B_index]] <-
          do.call(rbind.data.frame,
                  lapply(res.all[["Yang"]],
                         function(x) x[x$Model == "MM1b", ]))
        re.psi.all[[B_index]] <-
          do.call(rbind.data.frame, res.all[["M_psi"]])
        re.aic.all[[B_index]] <-
          do.call(rbind.data.frame, res.all[["AIC_sel"]])
      } # end loop over batches
      
      re.yang.all <- do.call(rbind.data.frame, re.yang.all)
      re.psi.all <- do.call(rbind.data.frame, re.psi.all)
      re.aic.all <- do.call(rbind.data.frame, re.aic.all)
      
      tab.sum.i <-
        rbind.data.frame(
          get.sum(re.i = re.yang.all, model.name = "1st-order Markov model", 
                  est.name = "MLE"),
          do.call(rbind.data.frame,
                  lapply(c("unadjusted", "BC"),
                         function(x)
                           get.sum(re.i = subset(re.psi.all, Estimator == x),
                                   model.name = paste0("M", num.s - 2, num.s), 
                                   est.name = x))),
          do.call(rbind.data.frame,
                  lapply(c("BC"),
                         function(x)
                           get.sum(re.i = subset(re.aic.all, Estimator == x),
                                   model.name = "AIC selection", est.name = x))))
      tab.sum.i[["true_pc"]] <- dat.sim.list$parm[count, "pc"]
      tab.sum.i[["true_nc"]] <- dat.sim.list$parm[count, "pc"]*
        dat.sim.list$parm[count, "N"]
      tab.sum.i <- data.frame("gen_model" = gen.m.s,
                              "Scenario" = paste0("T", num.s, "_S", count),
                              tab.sum.i)
      tab.sum.all[[count1]] <- tab.sum.i
      count = count + 1
      count1 = count1 + 1
    }
  }
  tab.sum.mat <- do.call(rbind.data.frame, tab.sum.all)
  tab.sum.final <- tab.sum.mat[,c("gen_model", "Scenario", "true", "true_pc",
                                  "Model", "Estimator", "rel.bias", 
                                  "cover", "width", "num.keep")]
  tab.sum.final[,c("cover")] <- 
    format(round(tab.sum.final[,c("cover")]*100, 1),
           nsmalls = 1)
  tab.sum.final[,c("true_pc", "rel.bias", "width")] <-
    format(round(tab.sum.final[,c("true_pc", "rel.bias", "width")], 3),
           nsmalls = 3)
  tab.sum.excel[[tab.name]] <- tab.sum.final
} # end loop over tables




#########################################
# generate Table 3 
#########################################
tab.name.vec <- c("Tab3")
for(i.tab in 1:length(tab.name.vec)){
  tab.name <- tab.name.vec[i.tab]
  tab.sum.all <- list()
  count1 = 1
  for(num.s in 6){
    ## read in simulated datasets
    load(paste0("./data/dat_sparse_M", num.s - 2, num.s, ".rda")) # dat.sim.list
    gen.m.s = paste0("M", num.s - 2, num.s)
    N.true.vec <- dat.sim.list$parm$N
    count = 1
    for(N.true in N.true.vec){
      re.yang.all <- re.psi.all <- re.aic.all <- list()
      for(B_index in 1:1000){
        load(paste0("./inter_res/Res_sparse_M", num.s-2, num.s, 
                    "_N", N.true, "_B", B_index, ".rda"))
        re.yang.all[[B_index]] <-
          do.call(rbind.data.frame,
                  lapply(res.all[["Yang"]],
                         function(x) x[x$Model == "MM1b", ]))
        re.psi.all[[B_index]] <-
          do.call(rbind.data.frame, res.all[["M_psi"]])
        re.aic.all[[B_index]] <-
          do.call(rbind.data.frame, res.all[["AIC_sel"]])
      } # end loop over batches
      
      re.yang.all <- do.call(rbind.data.frame, re.yang.all)
      re.psi.all <- do.call(rbind.data.frame, re.psi.all)
      re.aic.all <- do.call(rbind.data.frame, re.aic.all)
      
      tab.sum.i <-
        rbind.data.frame(
          get.sum(re.i = re.yang.all, model.name = "1st-order Markov model", 
                  est.name = "MLE"),
          do.call(rbind.data.frame,
                  lapply(c("unadjusted", "BC"),
                         function(x)
                           get.sum(re.i = subset(re.psi.all, Estimator == x),
                                   model.name = paste0("M", num.s - 2, num.s), 
                                   est.name = x))),
          do.call(rbind.data.frame,
                  lapply(c("BC"),
                         function(x)
                           get.sum(re.i = subset(re.aic.all, Estimator == x),
                                   model.name = "AIC selection", est.name = x))))
      tab.sum.i[["true_pc"]] <- dat.sim.list$parm[count, "pc"]
      tab.sum.i[["true_nc"]] <- dat.sim.list$parm[count, "pc"]*
        dat.sim.list$parm[count, "N"]
      tab.sum.i <- data.frame("gen_model" = gen.m.s,
                              "Scenario" = paste0("T", num.s, "_S", count),
                              tab.sum.i)
      tab.sum.all[[count1]] <- tab.sum.i
      count = count + 1
      count1 = count1 + 1
    }
  }
  tab.sum.mat <- do.call(rbind.data.frame, tab.sum.all)
  tab.sum.final <- tab.sum.mat[,c("gen_model", "Scenario", "true", "true_pc",
                                  "Model", "Estimator", "rel.bias", 
                                  "cover", "width", "num.keep")]
  tab.sum.final[,c("cover")] <- 
    format(round(tab.sum.final[,c("cover")]*100, 1),
           nsmalls = 1)
  tab.sum.final[,c("true_pc", "rel.bias", "width")] <-
    format(round(tab.sum.final[,c("true_pc", "rel.bias", "width")], 3),
           nsmalls = 3)
  tab.sum.excel[[tab.name]] <- tab.sum.final
} # end loop over tables
writexl::write_xlsx(tab.sum.excel, path = "./TabFigs/Tabs123S1.xlsx")
