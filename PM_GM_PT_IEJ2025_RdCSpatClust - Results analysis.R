########################################################################################################
##########                            Replication code for                                    ##########
##########    "The Geography of the Italian Citizenship Income:                               ##########
##########    the Role of Poverty and Inequality in Determining Spatial Heterogeneity across  ##########
##########    the Italian Municipalities" (Italian Economic Journal, 2025)                    ##########                                                    ##########
##########            Paolo Maranzano (UniMiB)                                                ##########
##########            Gianluca Monturano (UniBa)                                              ##########
##########            Pasquale Tridico (UniRoma3)                                             ##########
##########      Part II: analysis of the empirical results                                    ##########
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





#############################################
############### Data analysis ###############
#############################################

##### Binding
ICs <- ICs %>%
  mutate(Spec_lab = paste0("Spec.",Spec),
         N_ngbh_lab = paste0("Neighbors = ",N_ngbh),
         N_ngbh_lab = factor(N_ngbh_lab,
                             levels = c("Neighbors = 25","Neighbors = 50","Neighbors = 100"),
                             labels = c("Neighbors = 25","Neighbors = 50","Neighbors = 100"),ordered = T))

##### Select the year of interest
YR <- 2020

cols <- c("Spec.1" = "black", "Spec.2" = "#3399FF", "Spec.3" = "#FF3333",
          "Spec.4" = "#006633","Spec.5" = "#0000FF", "Spec.6" = "#FF9933")


##### Idnetificazione dei modelli con minori ICs
# Indicazione 5: BIC suggerisce modelli più parsimoniosi (minor numero di gruppi) rispetto HQ e AIC
# Indicazione 6: Secondo BIC, i modelli migliori sono tra 20, 25 e 26 gruppi
#       --> Parsimonia perchè le performance sono molto simili
# Indicazione 7: con il passare del tempo, il numero ottimo di gruppi aumenta (quasi raddoppia)
#       --> Aumento dell'eterogeneità degli effetti
ICs_min <- ICs %>%
  group_by(Year,Spec,Response,N_ngbh) %>%
  mutate(minBIC = min(BIC),minAIC = min(AIC),minHQ = min(HQ)) %>%
  ungroup() %>%
  add_rownames()

ICs_min[which(ICs_min$BIC == ICs_min$minBIC),] %>% arrange(BIC) %>% View()
ICs_min[which(ICs_min$AIC == ICs_min$minAIC),] %>% arrange(AIC)
ICs_min[which(ICs_min$HQ == ICs_min$minHQ),] %>% arrange(HQ)

##### Performance plots
### 1. Information criteria
# Indicazione 1: Ngbh = 25 perchè abbiamo i minori ICs, per tutti i ICs
p1 <- ICs %>%
  filter(N_grps >= 5) %>%
  pivot_longer(cols = c("BIC"), names_to = "IC", values_to = "Value") %>%
  ggplot() +
  geom_line(mapping = aes(x = N_grps, y = Value, col = Spec_lab)) +
  geom_point(mapping = aes(x = N_grps, y = Value, col = Spec_lab)) +
  ggh4x::facet_grid2(rows = vars(Year), cols = vars(as.factor(N_ngbh_lab)), scales = "free", independent = "y") +
  geom_vline(data = ICs_min %>% filter(BIC == minBIC),
             mapping = aes(xintercept = N_grps, col = Spec_lab), linewidth = 1.2) +
  # ggh4x::facet_grid2(rows = vars(Year), cols = vars(as.factor(N_ngbh))) +
  scale_color_manual("Specification", values = cols) +
  labs(title = "Bayesian Information Criterion (BIC) by number of clusters, neighbors, and year",
       x = "Number of clusters", y = "BIC") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 20),
        axis.title = element_text(size = 14))
ggpubr::ggexport(plotlist = list(p1), filename = "BIC_values.png", width = 2000, height = 1400, res = 150)

### 2. One-group-more gain
# Indicazione 2: dopo 15 gruppi il guadagno diventa marginale intorno allo 0
ICs %>%
  pivot_longer(cols = c("BIC","AIC","HQ"), names_to = "IC", values_to = "Value") %>%
  group_by(Spec_lab,N_ngbh_lab,Response,IC) %>%
  mutate(gain = c(NA,diff(Value))) %>%
  filter(N_grps >= 5) %>%
  ggplot() +
  geom_line(mapping = aes(x = N_grps, y = gain, col = as.factor(Spec_lab)), linewidth = 1.1) +
  geom_hline(mapping = aes(yintercept = 0), col = "black", linewidth = 1) +
  ggh4x::facet_grid2(rows = vars(Year), cols = vars(N_ngbh_lab), scales = "free", independent = "y") +
  scale_color_manual("Specification", values = cols) +
  labs(title = "One-group-more gain for Bayesian Information Criterion (BIC) by number of clusters, neighbors, and year",
       x = "Number of clusters", y = "BIC") +
  theme_bw()

### 3. Gain w.r.t. the pooled (K=1) regression
# Indicazione 3: Rispetto al pooled, tutte le clustered regression fanno meglio
# Indicazione 4: il guadagno diventa nullo dopo 15 gruppi
ICs %>%
  filter(Year == YR) %>%
  pivot_longer(cols = c("BIC","AIC","HQ"), names_to = "IC", values_to = "Value") %>%
  group_by(Year,Spec,N_ngbh,Response,IC) %>%
  mutate(Index = Value/Value[N_grps==1]*100) %>%
  ggplot() +
  geom_line(mapping = aes(x = N_grps, y = Index, col = as.factor(Spec)), linewidth = 1.1) +
  ggh4x::facet_grid2(rows = vars(IC), cols = vars(N_ngbh), scales = "free", independent = "y") +
  theme_bw()






##### Ristima del modello ottimo
# Grid model e specifications
Years_RDC <- 2019:2022
Specs_idx <- c(1:6)
N_grps <- 1:50
N_ngbh <- c(25,50,100)
Models <- c("poisson")
Grid_models <- expand.grid(Years_RDC,Specs_idx,N_grps,N_ngbh,Models)
shape <- data_sp %>%
  filter(Year == YR) %>%
  dplyr::select(ISTAT_Code_LAU)

ICs_min %>%
  filter(BIC == minBIC,Spec==3, N_ngbh == 25) %>%
  View()


# Select the optimal model parameters
optim_idx <- which(Grid_models$Var3 == 25 & Grid_models$Var4 == 25 & Grid_models$Var2 == 3)


optim_idx <- ICs_min %>%
  filter(BIC == minBIC,Spec==3, N_ngbh == 25) %>%
  dplyr::select(rowname) %>%
  pull()

##### Define colors scale by variable
Cols_beta <- vector(mode = "list", length = 5)
Cols_beta_m <- vector(mode = "list", length = 4)
Cols_beta_m_j <- vector(mode = "list", length = 3)

##### Graphical analysis of clusters
for (spec_num in 1:5) {

  cat(paste0("Specification ",spec_num," of ",5,": estimating the optimal model for each year","\n"))
  # spec_num <- 1

  estimate <- TRUE
  optim_idx <- ICs_min %>%
    filter(BIC == minBIC, N_ngbh == 25, Spec == spec_num)
  n_mods <- dim(optim_idx)[1]
  fit_optim <- vector(mode = "list", length = n_mods)
  Years <- c(2019,2020,2021,2022)
  variables <- c("Gini_con_0_lag1","Redditi_Procapite_lag1","share_of_poverty_lag1")
  variables_label <- c("Gini index (lagged)","Per capita income (lagged)","Share of poverty (lagged)")
  plot_res <- vector(mode = "list", length = n_mods)


  for (m in 1:n_mods) {

    # m <- 1

    ##### Estimation
    cat(paste0("Model ",m," of ",n_mods,": estimating","\n"))
    fit_optim[[m]] <- SCR_OptimEstim(Dataset = RDC_interp_masked,
                                     Grid_models = ICs,
                                     optim_idx = as.numeric(optim_idx[m,"rowname"]),
                                     Specs = Specs)
    mod <- fit_optim[[m]]

    ##### Estimating optimal model for each specification and year + mapping
    Coefs <- vector(mode = "list", length = length(variables))
    for (j in 1:length(variables)) {

      # j <- 1

      # Variable position
      var_idx <- which(colnames(fit_optim[[m]]$X_sp_opt)==variables[j])-1

      if (!is_empty(var_idx)) {
        # Attach clustering results
        if (variables[j] == "Redditi_Procapite_lag1") {
          # Redditi pro capite are measured in €, we report estimates as € each 1000 inhabitants
          nnn <- mod$Y_sp_opt %>%
            mutate(Clust_opt = mod$fit_optim$group,
                   NumFam = exp(mod$Offset_sp_opt$log_NucleiFamTot),
                   Betas_opt = mod$fit_optim$sBeta[,var_idx]*100*1000,
                   SEBetas_opt = sqrt(mod$fit_optim$VCov[,,mod$fit_optim$group][var_idx,var_idx,])*100*1000,
                   tBetas_opt = mod$fit_optim$sBeta[,var_idx]/sqrt(mod$fit_optim$VCov[,,mod$fit_optim$group][var_idx,var_idx,]))
          nnn2 <- mod$X_sp_opt %>%
            dplyr::select(ISTAT_Code_LAU,any_of(variables)) %>% st_drop_geometry()
        } else {
          # Gini and share of poverty are 0-100 percentages
          nnn <- mod$Y_sp_opt %>%
            mutate(Clust_opt = mod$fit_optim$group,
                   NumFam = exp(mod$Offset_sp_opt$log_NucleiFamTot),
                   Betas_opt = mod$fit_optim$sBeta[,var_idx]*100,
                   SEBetas_opt = sqrt(mod$fit_optim$VCov[,,mod$fit_optim$group][var_idx,var_idx,])*100,
                   tBetas_opt = mod$fit_optim$sBeta[,var_idx]/sqrt(mod$fit_optim$VCov[,,mod$fit_optim$group][var_idx,var_idx,]))
          nnn2 <- mod$X_sp_opt %>%
            dplyr::select(ISTAT_Code_LAU,any_of(variables)) %>% st_drop_geometry()
        }
        new_class <- data.frame(Clust_opt = names(mod$fit_optim$Beta[var_idx,]),
                                rank = rank(mod$fit_optim$Beta[var_idx,])) %>%
          mutate(Clust_opt = as.numeric(gsub("\\D", "", Clust_opt))) %>%
          arrange(rank)
        nnn <- left_join(x = nnn, y = new_class, by = "Clust_opt")
        nnn <- left_join(x = nnn, y = nnn2, by = "ISTAT_Code_LAU")
      }

      nnn <- nnn %>%
        mutate(Year = Years[m]) %>%
        mutate(Var = variables_label[j], .after = "Year")

      Cols_beta_m_j[[j]] <- nnn

    }

    Cols_beta_m[[m]] <- Cols_beta_m_j

  }

  Cols_beta[[spec_num]] <- Cols_beta_m

}



LineWidth <- c("Sign. Positive 0.1%" = 1.2,
               "Sign. Positive 1%" = 1.1,
               "Sign. Positive 5%" = 1.1,
               "Sign. Positive 10%" = 1,
               "Sign. Negative 0.1%" = 1.2,
               "Sign. Negative 1%" = 1.1,
               "Sign. Negative 5%" = 1.1,
               "Sign. Negative 10%" = 1,
               "Non signif." = 0.5)
LineCol <- c("Sign. Positive 0.1%" = "red",
             "Sign. Positive 1%" = "red",
             "Sign. Positive 5%" = "red",
             "Sign. Positive 10%" = "red",
             "Sign. Negative 0.1%" = "darkgreen",
             "Sign. Negative 1%" = "darkgreen",
             "Sign. Negative 5%" = "darkgreen",
             "Sign. Negative 10%" = "darkgreen",
             "Non signif." = "black")

CoefAggr_List <- vector(mode = "list", length = 5)
for (spec_num in 1:5) {
  for (j in 1:length(variables)) {

    # spec_num <- 1
    # j <- 1

    cat(paste0("Specification ",spec_num," of 5: mapping ",variables_label[j],": ",j," of ",length(variables_label),"\n"))
    Beta_sf <- bind_rows(Cols_beta[[spec_num]][[1]][[j]],
                         Cols_beta[[spec_num]][[2]][[j]],
                         Cols_beta[[spec_num]][[3]][[j]],
                         Cols_beta[[spec_num]][[4]][[j]])

    Beta_sf_grp <- Beta_sf %>%
      group_by(Year,tBetas_opt) %>%
      summarise(geometry = sf::st_union(geometry)) %>%
      ungroup() %>%
      mutate(Signif = case_when(
        # tBetas_opt >= qnorm(0.999) ~ "Sign. Positive 0.1%",
        # tBetas_opt >= qnorm(0.99) & tBetas_opt < qnorm(0.999) ~ "Sign. Positive 1%",
        # tBetas_opt >= qnorm(0.95) & tBetas_opt < qnorm(0.995) ~ "Sign. Positive 5%",
        tBetas_opt >= qnorm(0.95) ~ "Sign. Positive 5%",
        # tBetas_opt >= qnorm(0.90) & tBetas_opt < qnorm(0.975) ~ "Sign. Positive 10%",
        # tBetas_opt <= qnorm(0.0005) ~ "Sign. Negative 0.1%",
        # tBetas_opt <= qnorm(0.005) & tBetas_opt > qnorm(0.0005) ~ "Sign. Negative 1%",
        # tBetas_opt <= qnorm(0.025) & tBetas_opt > qnorm(0.005) ~ "Sign. Negative 5%",
        tBetas_opt <= qnorm(0.05) ~ "Sign. Negative 5%",
        # tBetas_opt <= qnorm(0.10) & tBetas_opt > qnorm(0.025) ~ "Sign. Negative 10%",
        TRUE ~ "Non signif."
      )
      )

    p <- Beta_sf %>%
      ggplot() +
      geom_sf(mapping = aes(fill = round(Betas_opt,1)), color = NA) +
      geom_sf(data = Beta_sf_grp, aes(linewidth = Signif, col = Signif), alpha = 0) +
      facet_wrap(~ Year) +
      coord_sf(ylim = c(36.5, 47.1), expand = FALSE) +
      scale_fill_gradient2(latex2exp::TeX(paste0("\\beta"," coefficient (% \\Delta)")),
                           low = "darkgreen",mid = "white", high = "red", midpoint = 0) +
      scale_linewidth_manual("Significance", values = LineWidth) +
      scale_color_manual("Significance", values = LineCol) +
      labs(title = paste0("Specification n.",spec_num,": ",variables_label[j])) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold",size = 20)) +
      guides(color = guide_legend(nrow = 3))
    ggpubr::ggexport(plotlist = list(p), filename = paste0("Plots_R1/","Spec",spec_num,"_",variables[j],".png"), width = 2000, height = 1400, res = 150)

    if (spec_num == 1) {
      Beta_sf_Aggr <- Beta_sf %>%
        st_drop_geometry() %>%
        group_by(Var,rank,Year) %>%
        rename(Predictor = any_of(variables[j])) %>%
        summarise(
          Nuclei_pesati_RdC = mean(Nuclei_pesati_RdC,na.rm=T),
          Betas_opt = mean(Betas_opt,na.rm=T),
          SEBetas_opt = mean(SEBetas_opt,na.rm=T),
          tBetas_opt = mean(tBetas_opt,na.rm=T),
          Predictor = mean(Predictor,na.rm=T)
        ) %>%
        mutate(Signif = case_when(
          tBetas_opt >= qnorm(0.95) ~ "Sign. Positive 5%",
          tBetas_opt <= qnorm(0.05) ~ "Sign. Negative 5%",
          TRUE ~ "Non signif.")
        )
      CoefAggr_List[[spec_num]][[j]] <- Beta_sf_Aggr
      p2 <- Beta_sf_Aggr %>%
        ggplot(mapping = aes(x = Predictor, y = Betas_opt, col = Signif)) +
        geom_errorbar(mapping = aes(x = Predictor, y = Betas_opt, ymin = Betas_opt-qnorm(0.975)*SEBetas_opt, ymax = Betas_opt+qnorm(0.975)*SEBetas_opt)) +
        geom_point(size = 2) +
        geom_hline(yintercept = 0) +
        facet_wrap(~ Year) +
        theme_bw() +
        labs(x = paste0("Average ",variables_label[j]),
             y = latex2exp::TeX(paste0("\\beta"," coefficient (% \\Delta) for ",variables_label[j])),
             title = paste0("Cluster-wise estimated coefficient by ",variables_label[j]))
      ggpubr::ggexport(plotlist = list(p2), filename = paste0("Plots_R1/AvgByGrp_","Spec",spec_num,"_",variables[j],".png"), width = 1200, height = 1000, res = 150)
    }

  }

}

spec_num <- 1
p3 <- bind_rows(CoefAggr_List[[spec_num]][[1]],
                CoefAggr_List[[spec_num]][[2]],
                CoefAggr_List[[spec_num]][[3]]) %>%
  ggplot(mapping = aes(x = Predictor, y = Betas_opt, col = Signif)) +
  geom_errorbar(mapping = aes(x = Predictor, y = Betas_opt, ymin = Betas_opt-qnorm(0.975)*SEBetas_opt, ymax = Betas_opt+qnorm(0.975)*SEBetas_opt)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0) +
  scale_color_manual("Significance", values = LineCol) +
  ggh4x::facet_grid2(rows = vars(Year), cols = vars(Var), scales = "free") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold",size = 18),
        axis.title = element_text(size = 16), strip.text = element_text(size = 16)) +
  labs(x = paste0("Average predictor"),
       y = latex2exp::TeX(paste0("\\beta"," coefficient (% \\Delta)")),
       title = paste0("Spec ",spec_num,": ","Confidence interval at 5% for the cluster-wise estimated coefficient by socio-economic drivers"))
ggpubr::ggexport(plotlist = list(p3), filename = paste0("Plots_R1/AvgByGrp_","Spec",spec_num,"_",variables[j],"_Pooled.png"), width = 2000, height = 1400, res = 150)



########################################################################################
########## Tabella per esercizio empirico sulla distribuzione dei 23 miliardi ##########
########################################################################################
spec_num <- 1
# Gini
j <- 1
Tab_sf <- bind_rows(Cols_beta[[spec_num]][[1]][[j]],
                    Cols_beta[[spec_num]][[2]][[j]],
                    Cols_beta[[spec_num]][[3]][[j]],
                    Cols_beta[[spec_num]][[4]][[j]])
save(Tab_sf,file = "Dati_Es23mld.RData")



###############################################################
########## Pooled Poisson regression (no clustering) ##########
###############################################################
library(xtable)
library(texreg)
source("SCR-main/AuxFn - PoolingPoisson_Estim.R", encoding = 'UTF-8')
fit_pooled <- vector(mode = "list", length = 5)
fit_pooled[[1]] <- PoolingPoisson_Estim(Dataset = RDC_interp, Year = 2019, Specs = Specs)$fit_optim
fit_pooled[[2]] <- PoolingPoisson_Estim(Dataset = RDC_interp, Year = 2020, Specs = Specs)$fit_optim
fit_pooled[[3]] <- PoolingPoisson_Estim(Dataset = RDC_interp, Year = 2021, Specs = Specs)$fit_optim
fit_pooled[[4]] <- PoolingPoisson_Estim(Dataset = RDC_interp, Year = 2022, Specs = Specs)$fit_optim
fit_pooled[[5]] <- PoolingPoisson_Estim(Dataset = RDC_interp, Year = 2019:2022, Specs = Specs)$fit_optim

### Store in latex
Mods_gof <- GLM_GoF(GLModel_list = fit_pooled)$GLModel_GoF
Mods_gof <- round(Mods_gof,3)
texreg_table <- texreg(l = fit_pooled,
                       file = "Output.tex",
                       caption = "Specification 1 of Poisson regression: pooled estimates (no spatial clustering) of regression coefficients for year-specific regression (columns 2 to 5) and time-pooled (all data from 2019 to 2022) regression (column 6)",
                       label = "Tab:Pooled_Poisson",
                       caption.above = TRUE,
                       custom.gof.rows = list(
                         "Resid.SD" = Mods_gof[1,],
                         "Num. \\ obs" = Mods_gof[2,],
                         "Num. \\ pars." = Mods_gof[3,],
                         "Log Likelihood" = Mods_gof[4,],
                         "AIC" = Mods_gof[5,],
                         "BIC" = Mods_gof[6,],
                         "Nagelkerke's Pseudo R$^2$" = Mods_gof[7,]
                       ),
                       custom.model.names = c("2019","2020","2021","2022","2019:2022")
)
