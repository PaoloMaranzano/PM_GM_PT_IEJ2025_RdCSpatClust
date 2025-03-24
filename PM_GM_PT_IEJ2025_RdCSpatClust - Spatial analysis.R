########################################################################################################
##########                            Replication code for                                    ##########
##########    "The Geography of the Italian Citizenship Income:                               ##########
##########    the Role of Poverty and Inequality in Determining Spatial Heterogeneity across  ##########
##########    the Italian Municipalities" (Italian Economic Journal, 2025)                    ##########                                                    ##########
##########            Paolo Maranzano (UniMiB)                                                ##########
##########            Gianluca Monturano (UniBa)                                              ##########
##########            Pasquale Tridico (UniRoma3)                                             ##########
##########    Part I: estimating the SCR-Poisson models                                       ##########
########################################################################################################

##### Setup
setwd("~/GitHub/PM_GM_PT_IEJ2025_RdCSpatClust")

##### Libraries
library(tidyverse)
library(sf)
library(spatialreg)
library(readxl)
library(spdep)
library(ggpubr)

##### Source SCR codes from Sugasawa 2021
source("SCR-main/SCR-function.R", encoding = 'UTF-8')
source("SCR-main/AuxFn - SCR_OptimEstim.R", encoding = 'UTF-8')

##### Dataset
load("PM_GM_PT_IEJ2025_RdCSpatClust - Data.Rdata")

##### Shapefile
load("Shape_comuni.RData")

##### Specifications
Specs <- read_excel("Meta_DB_RDC.xlsx", sheet = "Spec Paolo")


#########################################################################
########## Clustered spatial regression: automatic computation ##########
#########################################################################
# Years of interest
Years_RDC <- 2019:2022

# Models specifications
Specs_idx <- 1:5

# Clustering parameters
N_grps <- 1:50
N_ngbh <- c(25,50,100)
Models <- c("poisson")

# Paramater grid
Grid_models <- expand.grid(Years_RDC,Specs_idx,N_grps,N_ngbh,Models)
fit <- vector(mode = "list", length = dim(Grid_models)[1])
timing <- matrix(data = NA ,nrow = dim(Grid_models)[1], ncol = 2)


#############################################
############### Main for loop ###############
#############################################

for (m in 1:dim(Grid_models)[1]) {
  # m <- 1

  t1 <- Sys.time()
  cat(paste0("Iteration: ",m," of ",dim(Grid_models)[1], ": ",
             "Year = ",Grid_models[m,1]," -- Spec = ",Grid_models[m,2]," -- ",
             "Nghb = ",Grid_models[m,4]," -- NumGrps = ",Grid_models[m,3],
             " started at in ",t1,"\n"))

  if (Grid_models[m,1] != 2022) {
    shape <- data_sp %>%
      filter(Year == Grid_models[m,1]) %>%
      dplyr::select(ISTAT_Code_LAU)
  } else {
    shape <- data_sp %>%
      filter(Year == 2021) %>%
      dplyr::select(ISTAT_Code_LAU)
  }

  fit[[m]] <- SCR_OptimEstim(Dataset = RDC_interp_masked,
                             Grid_models = Grid_models,
                             optim_idx = m,
                             Specs = Specs)$fit_optim

  t2 <- Sys.time()
  timing[m,] <- c(as.character(t1),as.character(t2))
  cat(paste0("Iteration: ",m," of ",dim(Grid_models)[1], ": ",
             "Year = ",Grid_models[m,1]," -- Spec = ",Grid_models[m,2]," -- ",
             "Nghb = ",Grid_models[m,4]," -- NumGrps = ",Grid_models[m,3],
             " completed in ",round(t2-t1,2)," secs.","\n"))

  ##### Save partial results
  n_blocks <- 20
  if (m %in% round(seq(from = 0, to = dim(Grid_models)[1], length.out = n_blocks),0)) {
    save(Grid_models,timing,fit, file = paste0("PartOut_",Grid_models[m,1],"_",m,".RData"))
  }

  gc()

}

##### Bind information criteria for each model
ICs <- matrix(data = NA ,nrow = dim(Grid_models)[1], ncol = 3)
for (m in 1:dim(Grid_models)[1]) {
  # fit_red[[m]] <- fit[[m]]$fit_optim
  ICs[m,] <- c(fit[[m]]$fit_optim$BIC,fit[[m]]$fit_optim$AIC,fit[[m]]$fit_optim$HQ)
}
ICs <- cbind(Grid_models,ICs)
colnames(ICs) <- c("Year","Spec","N_grps","N_ngbh","Response","BIC","AIC","HQ")

##### Store complete results
save(Grid_models,ICs,timing,fit,file = paste0("OutputComplete_Specs_",
                                              paste(Specs_idx,collapse = ""),"_",
                                              paste(Years_RDC,collapse = ""),".RData"))
