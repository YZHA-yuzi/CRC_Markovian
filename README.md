README
================

In this document, we illustrate the modeling framework proposed in the
paper titled “New capture-recapture models of behavioral response for
estimating the size of a closed animal population” with a real dataset
which was collected over four consecutive nights in 1983 in New Zealand
for capturing brushtail possums. We also describe how to produce
simulation results included in the paper in this document.

# Analyze the four-sample possum data using the proposed modeling framework

Read in self-defined functions

``` r
source("FUNs.R")
```

Read in possum data

``` r
dat.his <- readxl::read_xlsx(path = "./data/ThreeRealDatasets.xlsx", 
                             sheet = "Possum")
```

Prepare data used for the model fitting

The possum data is summarized in a frequency table with
![2^4-1=15](https://latex.codecogs.com/png.latex?2%5E4-1%3D15 "2^4-1=15")
categories.

``` r
# the number of trapping occasions T = 4
nstreams = 4

hist.sel.num = paste0("1", paste0(rep(0, nstreams - 2), collapse = ""), "1")
hist.sel.den = paste0("1", paste0(rep(0, nstreams - 1), collapse = ""))
hist.key = paste0(paste0(rep(0, nstreams - 1), collapse = ""), "1")

hist.all <- get_allhist(nstreams)
dat.fake <- data.frame("prof" = apply(hist.all, 1, paste, collapse = ""))
dat.full <- left_join(dat.fake, dat.his, by = c("prof" = "hist"))
dat.full$count[is.na(dat.full$count)] <- 0
dat_aggre <- as.data.frame(matrix(dat.full$count, nrow = 1))
colnames(dat_aggre) <- paste0("n", dat.full$prof)
dat_aggre
```

    ##   n1111 n1110 n1101 n1100 n1011 n1010 n1001 n1000 n0111 n0110 n0101 n0100 n0011
    ## 1    11     5     5     4     6    19     6    29     3     5     4    19     1
    ##   n0010 n0001
    ## 1    20     1

To implement the proposed model selection precedure, we fit obtain AIC
of the two candidate models:
![M\_{2,4}](https://latex.codecogs.com/png.latex?M_%7B2%2C4%7D "M_{2,4}")
and
![M^\*\_{2,4}](https://latex.codecogs.com/png.latex?M%5E%2A_%7B2%2C4%7D "M^*_{2,4}")
using the R function `fit.markovian.assump`.

This function contains 4 arguments:

- `nstreams` specifies the number of time points, i.e.,
  ![T](https://latex.codecogs.com/png.latex?T "T") in the manuscript.
- `dat_aggre` is the frequency table summarizing the data, the number of
  cells in this table is
  ![2^T-1](https://latex.codecogs.com/png.latex?2%5ET-1 "2^T-1").
- `assump` is a character string which specifies the desired
  Markovian-type assumption. For example, `assump` = “2days” indicates
  ![j = 2](https://latex.codecogs.com/png.latex?j%20%3D%202 "j = 2").
  This assumption implies that after 2 time points, animals forget their
  responses to the prior trapping occasions.
- `obscons` takes a logical value, `obscons = FALSE` indicates only
  untestable constraints defined in Equation (3) are imposed, while
  `obscons = TRUE` indicates both untestable constraints defined in
  Equation (3) and testable constraints defined in Equation (7) are
  imposed.  
- `hessian` takes a logical value, `hessian = TRUE` indicates NO hessian
  matrix is computed when maximizing the conditional likelihood derived
  based on model (1) in the manuscript, while `hessian = FALSE`
  indicates hessian matrix is computed. The default value is `FALSE`.

The output of this function contains:

- `Nhat` = the maximum likelihood estimator (MLE) of
  ![N](https://latex.codecogs.com/png.latex?N "N") derived by
  numerically maximizing the conditional likelihood derived based on
  model (1) under the desired constraints.
- `pc` = MLE of
  ![p_c(\boldsymbol{\theta})](https://latex.codecogs.com/png.latex?p_c%28%5Cboldsymbol%7B%5Ctheta%7D%29 "p_c(\boldsymbol{\theta})")
  (i.e., the probability of being captured at least once) derived by
  numerically maximizing the conditional likelihood derived based on
  model (1) under the desired constraints.
- `AIC` = the AIC value of the fitted model.
- `num_parm` = the number of parameters of the fitted model.

Compute the AIC of the
![M\_{2, 4}](https://latex.codecogs.com/png.latex?M_%7B2%2C%204%7D "M_{2, 4}")
model for the possum data

``` r
# compute the number of constraints imposed by the M_{2,4} model
num.constraints(nstreams = 4, j.assumed = 2, testable = F)
re.m24 <- fit.markovian.assump(nstreams = nstreams, 
                               dat_aggre = dat_aggre,
                               assump = "2days",
                               obscons = F, hessian = T)
aic.m24 <- re.m24$AIC
aic.m24
```

    ## [1] 676.2008

Compute the AIC of the
![M^\*\_{2, 4}](https://latex.codecogs.com/png.latex?M%5E%2A_%7B2%2C%204%7D "M^*_{2, 4}")
model for the possum data

``` r
# compute the number of constraints imposed by the M*_{2,4} model
num.constraints(nstreams = 4, j.assumed = 2, testable = T)
re.m24star <- fit.markovian.assump(nstreams = nstreams, 
                               dat_aggre = dat_aggre,
                               assump = "2days",
                               obscons = T, hessian = T)
aic.m24star <- re.m24star$AIC
aic.m24star
```

    ## [1] 680.3874

The
![M\_{2, 4}](https://latex.codecogs.com/png.latex?M_%7B2%2C%204%7D "M_{2, 4}")
has the lower AIC compared to the
![M^\*\_{2,4}](https://latex.codecogs.com/png.latex?M%5E%2A_%7B2%2C4%7D "M^*_{2,4}"),
thus the
![M\_{2, 4}](https://latex.codecogs.com/png.latex?M_%7B2%2C%204%7D "M_{2, 4}")
is selected.

Obtain bias-corrected estimator in Equation (13) and 95% credible
interval based on the selected
![M\_{2,4}](https://latex.codecogs.com/png.latex?M_%7B2%2C4%7D "M_{2,4}")
model for the possum data.

The function `get_Nhat_psi` is used to obtain bias-corrected estimator.
Note that, since the
![M\_{2,4}](https://latex.codecogs.com/png.latex?M_%7B2%2C4%7D "M_{2,4}")
model only impose one untestable constraint, observed cell counts are
used.

``` r
# nc = the number of uniquely identified animals
# nlast = n_{00...01} 
# psi = (m1+0.75)/(m0+m1+1)
m0 = dat_aggre[["n1000"]]
m1 = dat_aggre[["n1001"]]
Nhat.bc <- ceiling(get_Nhat_psi(nc = sum(dat_aggre), 
                        nlast = dat_aggre[["n0001"]], 
                        psi = (m1+0.75)/(m0+m1+1)))
Nhat.bc
```

    ## [1] 143

As discussed in the manuscript, the lower bound of the credible interval
is computed using the Beta(1, 0) prior, while the uppoer bound is
computed using the Beta(0.5, 0.5) prior.

The function `Dir_CI_alter_Fmax` is used to obtain 95% credible interval
associated the bias-corrected estimator.

- The argument `dat_aggre` is the frequency table summarizing the data,
  the number of cells in this table is
  ![2^T-1](https://latex.codecogs.com/png.latex?2%5ET-1 "2^T-1").
- `nstreams` specifies the number of time points, i.e.,
  ![T](https://latex.codecogs.com/png.latex?T "T") in the manuscript.
- The argument `n.post` specifies the number of posterior samples
  generated to compute credible intervals.
- The arguments `a` and `b` are parameters in Beta priors. When
  computing the credible interval, we set `a=1` and `b=0` to obtain
  lower bound and `a=0.5` and `b=0.5` to obtain upper bound.
- `hist.sel.num` specifies the capture profile of
  ![m_1](https://latex.codecogs.com/png.latex?m_1 "m_1") in Equation
  (12). For example `hist.sel.num="1001"` for the possum data.
- `hist.sel.num` specifies the capture profile of
  ![m_0](https://latex.codecogs.com/png.latex?m_0 "m_0") in Equation
  (12). For example `hist.sel.den="0001"` for the possum data.

``` r
set.seed(2345)
# n.post = the number of posterior samples 
# hist.sel.num specifies the capture profile of m1
# hist.sel.den specifies the capture profile of m0
CI.lower <- ceiling(Dir_CI_alter_Fmax(df_aggre = dat_aggre, n.post = 500,
                        nstreams = nstreams,
                        a = 1, b = 0,
                        hist.sel.num = "1001",
                        hist.sel.den = "1000")[["CI"]][1])

CI.upper <- ceiling(Dir_CI_alter_Fmax(df_aggre = dat_aggre, n.post = 500,
                        nstreams = nstreams,
                        Dir_prior = "equprob",
                        a = 0.5, b = 0.5,
                        hist.sel.num = "1001",
                        hist.sel.den = "1000")[["CI"]][2])
CI <- c(CI.lower, CI.upper)
CI
res.mat <- data.frame("Nhat" = Nhat.bc, "lci" = CI[1], "uci" = CI[2])
pander(res.mat, style = "rmarkdown", split.table = Inf)
```

    ## [1] 138 164

| Nhat | lci | uci |
|:----:|:---:|:---:|
| 143  | 138 | 164 |

The
![M\_{2, 4}](https://latex.codecogs.com/png.latex?M_%7B2%2C%204%7D "M_{2, 4}")
yieled an estimate of 143 (95% CI: 138, 164).

# Produce simulation results and real data analysis results included in the manuscript

## Data generation

Simulated data used in simulation studies can be generated by running
the R script `SIM_data.R`.

## R scripts

- `FUNs.R` (self-defined functions)
- `SIM_data.R` (generate simulated data)
- `SIM_fitting_Tab1.R` (conduct simulations presented in Table 1)
- `SIM_fitting_Tab2.R` (conduct simulations presented in Table 2)
- `SIM_fitting_Tab3.R` (conduct simulations presented in Table 3)
- `SIM_fitting_TabS1.R` (conduct simulations presented in Table S1)
- `SIM_sumres.R` (summarize simulation results)
- `RealApp.R` (conduct real data analysis)
- `RealApp_Fig1.R` (generate Figure 1)

## Bash shell scripts

- `fit_Tab1.sh`
- `fit_Tab2.sh`
- `fit_Tab3.sh`
- `fit_TabS1.sh`
- `fit_RealApp.sh`

## Instructions for conducting the simulation study

There are two options to run those provided R scripts: (1) RStudio or
(2) Bash Shell scripts. When running R codes in RStudio, please create a
new R project and make sure R scripts are stored at the same folder
where the created R project `*.Rproj` is; when running R codes using
provided bash shell scripts, please set the working directory to the
folder where R scripts and bash shell scripts are located at.

Run-time was approximated using a computer with Apple M2 Pro and 32 Gb
memory.

Steps to produce simulation results presented in Section 4 are given
below.

**Step 1**. Generate data

Run `SIM_data.R` in RStudio; or execute the bash script `sim_data.sh`
using the command `bash sim_data.sh` in terminal.

All simulated datasets will be saved to the folder named `data` under
the working directory.

**Step 2**. Conduct simulation studies under various simulation
scenarios

Run `SIM_fitting_Tab1.R`, `SIM_fitting_Tab2.R`, `SIM_fitting_Tab3.R`,
and `SIM_fitting_TabS1.R` in RStudio (arguments `num_index`, `s_index`,
and `B_index` must be specified); or execute bash scripts `fit_Tab1.sh`,
`fit_Tab2.sh`, `fit_Tab3.sh`, and `fit_TabS1.sh`. The computation time
varies by simulation scenarios and the number of simulations ran in the
for-loop when the parallel computation is implemented.

Simulation results will be saved to a sub-directory named `inter_res`.
Tables 1 to 3, and S1 will be generated and saved to the folder
`TabFigs` by running `SIM_sumres.R` in RStudio.

## Instructions for conducting real data applications

Run `RealApp.R` in RStudio or execute bash scripts `fit_RealApp.sh`;
results will be saved to the sub-directory `inter_res`. To generate
Figure 1 which will be saved to the folder `TabFigs`, run
`RealApp_Fig1.R` in RStudio.
