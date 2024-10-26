############################################################################
# -----------------  R codes for the real data applications -------------- #
# NOTES:
# summarize real data analysis results to generate Figure 1  
# please make sure the working directory is the folder where 
# files "FUNs.R" is saved.
############################################################################

# load packages 
library(ggplot2)
library(readxl)
library(ggpubr)

#### Summarize real data analysis results ####
name.mat <- data.frame("type" = c("Snowshoe_hare", "Mouse", "Possum"),
                       "nstreams" = 6:4)
tab.appres <- list()
for(num_s in 1:3){
  nstreams = name.mat$nstreams[num_s]
  load(paste0("./inter_res/Res_App_T", nstreams, ".rda"))
  dat <- readxl::read_xlsx(path = "./data/ThreeRealDatasets.xlsx",
                           sheet = name.mat$type[num_s])
  tab.1 <- re.app.all$Alter[,c("Estimator","Nhat", "lci", "uci")]
  tab.1 <- 
    data.frame(tab.1, 
               "Model" = re.app.all$Alter_can$model[which.min(re.app.all$Alter_can$AIC)])
  tab.2 <- data.frame("Estimator" = "Yang_MM1b", 
                      re.app.all$Yang$Yang[1,c("Nhat", "lci", "uci")],
                      "Model" = "Yang_MLE")
  tab.comb <- data.frame("Scenario" = paste0("T=", nstreams), 
                         "nc" = sum(dat$count),
                         rbind.data.frame(tab.1, tab.2))
  tab.comb$Estimator <- factor(tab.comb$Estimator, 
                               levels = c("BC", "unadjusted", "Yang_MM1b"))
  tab.comb <- tab.comb[order(tab.comb$Estimator), ]
  tab.comb <- data.frame(tab.comb, 
                         "num_zerocounts" = sum(dat$count == 0),
                         "prop_zerocounts" = mean(dat$count == 0))
  tab.appres[[num_s]] <- tab.comb
}
tab.appres <- do.call(rbind.data.frame, tab.appres)


# -------------------- BEGIN GENERATING FIGURE 1 ---------------- #
#### Plot for T = 6 
tab.forpl1 <- subset(tab.appres, Scenario == "T=6")
xlab.vec.1 <- c("BC" = expression(paste("M*"["T-2,T"], 
                                        "\n (bias-corrected)", sep = "")),
                "unadjusted" = expression(paste("M*"["T-2,T"], 
                                                "\n (unadjusted)", sep = "")),
                "Yang_MM1b" = "1st-order Markov")
til.1 <- expression(paste("Snowshoe hare data (T=6), ", 
                          "n"["c"], "=68, ", 
                          "%zero counts=48", sep = ""))
pl.1 <- ggplot(tab.forpl1, 
               aes(x = Estimator, y = Nhat, ymin = lci, ymax = uci)) + 
  geom_pointrange(position = position_dodge(width=0.40)) + 
  geom_errorbar(position = position_dodge(width=0.40), 
                width = 0.1, linewidth = 0.8) + 
  scale_x_discrete(labels = xlab.vec.1) + 
  theme_bw() + 
  labs(title = til.1,
       x = "Estimator", 
       y = "Estimated N and 95% CI") +
  theme(axis.text.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(size = rel(1.3)),
        axis.title.x = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)))

#### Plot for T = 5
tab.forpl2 <- subset(tab.appres, Scenario == "T=5")
xlab.vec.2 <- c("BC" = expression(paste("M*"["T-2,T"], 
                                            "\n (bias-corrected)", sep = "")),
                "unadjusted" = expression(paste("M*"["T-2,T"], 
                                           "\n (unadjusted)", sep = "")),
                "Yang_MM1b" = "1st-order Markov")
til.2 <- expression(paste("Mouse data (T=5), ", 
                          "n"["c"], "=55, ", 
                          "%zero counts=32", sep = ""))
pl.2 <- ggplot(tab.forpl2, 
               aes(x = Estimator, y = Nhat, ymin = lci, ymax = uci)) + 
  geom_pointrange(position = position_dodge(width=0.40)) + 
  geom_errorbar(position = position_dodge(width=0.40), 
                width = 0.1, size = 0.8) + 
  scale_x_discrete(labels = xlab.vec.2) + 
  theme_bw() + 
  labs(title = til.2,
       x = "Estimator", 
       y = "Estimated N and 95% CI") +
  theme(axis.text.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(size = rel(1.3)),
        axis.title.x = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)))

#### Plot for T = 4
tab.forpl3 <- subset(tab.appres, Scenario == "T=4")
xlab.vec.3 <- c("BC" = expression(paste("M"["T-2,T"], 
                                            "\n (bias-corrected)", sep = "")),
                "unadjusted" = expression(paste("M"["T-2,T"], 
                                           "\n (unadjusted)", sep = "")),
                "Yang_MM1b" = "1st-order Markov")
til.3 <- expression(paste("Possum data (T=4), ", 
                          "n"["c"], "=138, ", 
                          "%zero counts=0", sep = ""))
pl.3 <- ggplot(tab.forpl3, 
               aes(x = Estimator, y = Nhat, ymin = lci, ymax = uci)) + 
  geom_pointrange(position = position_dodge(width=0.40)) + 
  geom_errorbar(position = position_dodge(width=0.40), 
                width = 0.1, size = 0.8) + 
  scale_x_discrete(labels = xlab.vec.3) + 
  theme_bw() + 
  labs(title = til.3,
       x = "Estimator", 
       y = "Estimated N and 95% CI") +
  theme(axis.text.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(size = rel(1.3)),
        axis.title.x = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)))

pl.comb <- ggpubr::ggarrange(pl.1, pl.2, pl.3, ncol = 2,
                             nrow = 2, labels = c("(A)", "(B)", "(C)"))
ggsave(filename = paste0("./TabFigs/Fig1.png"),
       height = 9, width = 12, units = "in", dpi = 600,
       plot = pl.comb, bg = "white")
