PoolingPoisson_Estim <- function(Dataset,Year,Specs) {
  ##### Define the dataset
  # Year <- 2019
  # Dataset <- RDC_interp
  # Specs <- Specs
  s <- 1
  
  # Response variable
  Y_name <- pull(Specs[which(Specs[,s+1] == "Response"),1])
  Y_df <- Dataset %>%
    filter(Anno %in% Year) %>%
    dplyr::select(ISTAT_Code_LAU = istat.x, Anno, all_of(Y_name))
  Y_df <- Y_df[which(complete.cases(Y_df)),]
  
  # Offset
  Offset_name <- "log_NucleiFamTot"
  Offset_df <- Dataset %>%
    rename(ISTAT_Code_LAU = istat.x) %>%
    filter(Anno %in% Year,
           ISTAT_Code_LAU %in% Y_df$ISTAT_Code_LAU) %>%
    dplyr::select(ISTAT_Code_LAU, Anno, all_of(Offset_name))
  Offset_df <- Offset_df[which(complete.cases(Offset_df)),]
  
  # Covariates
  Xlag_name <- pull(Specs[which(Specs[,s+1] == "Covariate Lagged"),1])
  X_lagged <- Dataset %>%
    rename(ISTAT_Code_LAU = istat.x) %>%
    filter(Anno %in% (Year - 1),
           ISTAT_Code_LAU %in% Y_df$ISTAT_Code_LAU) %>%
    mutate(Anno = Anno + 1) %>%
    dplyr::select(ISTAT_Code_LAU, Anno, all_of(Xlag_name))
  colnames(X_lagged)[-c(1,2)] <- paste0(colnames(X_lagged)[-c(1,2)],"_lag1")
  X_name <- pull(Specs[which(Specs[,s+1] == "Covariate"),1])
  X_cont <- Dataset %>%
    rename(ISTAT_Code_LAU = istat.x) %>%
    filter(Anno %in% Year,
           ISTAT_Code_LAU %in% Y_df$ISTAT_Code_LAU) %>%
    dplyr::select(ISTAT_Code_LAU, Anno, all_of(X_name))
  X_df <- full_join(x = X_cont, y = X_lagged, by = c("ISTAT_Code_LAU","Anno"))
  X_df <- X_df[which(complete.cases(X_df)),]
  
  YX_df <- inner_join(x = X_df, y = Y_df, by = c("ISTAT_Code_LAU","Anno"))
  YX_df <- inner_join(x = YX_df, y = Offset_df, by = c("ISTAT_Code_LAU","Anno"))
  summary(YX_df)
  
  X_df_complete <- YX_df[which(complete.cases(YX_df)),] %>%
    dplyr::select(colnames(X_df))
  Y_df_complete <- YX_df[which(complete.cases(YX_df)),] %>%
    dplyr::select(colnames(Y_df))
  Offset_df_complete <- YX_df[which(complete.cases(YX_df)),] %>%
    dplyr::select(colnames(Offset_df))
  dim(X_df_complete)
  dim(Y_df_complete)
  dim(Offset_df_complete)
  
  
  # Convert to matrx
  Y_mat <- Y_df_complete %>%
    dplyr::select(-c(ISTAT_Code_LAU,Anno)) %>%
    st_drop_geometry() %>%
    as.matrix()
  
  # Convert to matrx
  Offset_mat <- Offset_df_complete %>%
    dplyr::select(-c(ISTAT_Code_LAU,Anno)) %>%
    st_drop_geometry() %>%
    as.matrix()
  
  # Convert to matrx
  X_mat <- X_df_complete %>%
    dplyr::select(-c(ISTAT_Code_LAU,Anno)) %>%
    st_drop_geometry() %>%
    as.matrix()
  
  fit <- glm(Y_mat ~ X_mat, offset = Offset_mat, family="poisson")
  
  return(list(fit_optim = fit,
              X_mat_opt = X_mat, Y_mat_opt = Y_mat, Offset_mat_opt = Offset_mat))
}


GLM_GoF <- function(GLModel_list) {
  
  ##### Matrix to store results
  GLModel_GoF <- matrix(NA,ncol = 9, nrow = length(GLModel_list))
  
  ##### Compute the Goodness-of-fit statistics
  for (i in 1:length(GLModel_list)) {
    GLModel_GoF[i,1] <- as.numeric(performance::performance(GLModel_list[[i]])["RMSE"])
    GLModel_GoF[i,2] <- length(GLModel_list[[i]]$y)
    GLModel_GoF[i,3] <- length(GLModel_list[[i]]$coefficients)
    GLModel_GoF[i,4] <- as.numeric(stats::logLik(GLModel_list[[i]]))
    GLModel_GoF[i,5] <- as.numeric(performance::performance(GLModel_list[[i]])["AIC"])
    GLModel_GoF[i,6] <- as.numeric(performance::performance(GLModel_list[[i]])["BIC"])
    GLModel_GoF[i,7] <- as.numeric(performance::performance(GLModel_list[[i]])["R2_Nagelkerke"])
    GLModel_GoF[i,8] <- as.numeric(GLModel_list[[i]]$deviance)
    GLModel_GoF[i,9] <- as.numeric(GLModel_list[[i]]$null.deviance)
  }
  
  ##### Customize table
  colnames(GLModel_GoF) <- c(
    "Residuals standard deviation", "Numer of observations","Number of parameters","Log-likelihood",
    "AIC","BIC","Nagelkerke's Pseudo R2","Deviance","Deviance of the null model"
  )
  
  ##### Output
  GLModel_GoF <- as.data.frame(t(GLModel_GoF))
  return(list(GLModel_GoF = GLModel_GoF))
}
