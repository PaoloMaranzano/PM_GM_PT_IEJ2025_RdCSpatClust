SCR_OptimEstim <- function(Dataset,Grid_models,optim_idx,Specs) {
  
  m <- optim_idx
  
  ##### Define the dataset
  s <- Grid_models[m,2]
  Year <- Grid_models[m,1]
  
  # Response variable
  Y_name <- pull(Specs[which(Specs[,s+1] == "Response"),1])
  Y_df <- Dataset %>%
    filter(Anno == Year) %>%
    dplyr::select(ISTAT_Code_LAU = istat.x, all_of(Y_name))
  Y_df <- Y_df[which(complete.cases(Y_df)),]
  
  # Offset
  Offset_name <- "log_NucleiFamTot"
  Offset_df <- Dataset %>%
    rename(ISTAT_Code_LAU = istat.x) %>%
    filter(Anno == Year,
           ISTAT_Code_LAU %in% Y_df$ISTAT_Code_LAU) %>%
    dplyr::select(ISTAT_Code_LAU, all_of(Offset_name))
  Offset_df <- Offset_df[which(complete.cases(Offset_df)),]
  
  # Covariates
  Xlag_name <- pull(Specs[which(Specs[,s+1] == "Covariate Lagged"),1])
  X_lagged <- Dataset %>%
    rename(ISTAT_Code_LAU = istat.x) %>%
    filter(Anno == Year - 1,
           ISTAT_Code_LAU %in% Y_df$ISTAT_Code_LAU) %>%
    dplyr::select(ISTAT_Code_LAU, all_of(Xlag_name))
  colnames(X_lagged)[-1] <- paste0(colnames(X_lagged)[-1],"_lag1")
  X_name <- pull(Specs[which(Specs[,s+1] == "Covariate"),1])
  X_cont <- Dataset %>%
    rename(ISTAT_Code_LAU = istat.x) %>%
    filter(Anno == Year,
           ISTAT_Code_LAU %in% Y_df$ISTAT_Code_LAU) %>%
    dplyr::select(ISTAT_Code_LAU, all_of(X_name))
  X_df <- full_join(x = X_cont, y = X_lagged, by = "ISTAT_Code_LAU")
  X_df <- X_df[which(complete.cases(X_df)),]
  
  YX_df <- inner_join(x = X_df, y = Y_df, by = "ISTAT_Code_LAU")
  YX_df <- inner_join(x = YX_df, y = Offset_df, by = "ISTAT_Code_LAU")
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
  
  # Shapefile
  shape_red <- shape %>%
    filter(ISTAT_Code_LAU %in% unique(YX_df$ISTAT_Code_LAU))
  
  # Merge with shapefile
  Y_sp <- left_join(x = shape_red, y = Y_df_complete, by = "ISTAT_Code_LAU")
  # Convert to matrx
  Y_mat <- Y_sp %>%
    dplyr::select(-ISTAT_Code_LAU) %>%
    st_drop_geometry() %>%
    as.matrix()
  
  # Merge with shapefile
  Offset_sp <- left_join(x = shape_red, y = Offset_df_complete, by = "ISTAT_Code_LAU")
  # Convert to matrx
  Offset_mat <- Offset_sp %>%
    dplyr::select(-ISTAT_Code_LAU) %>%
    st_drop_geometry() %>%
    as.matrix()
  
  # Merge with shapefile
  X_sp <- left_join(x = shape_red, y = X_df_complete, by = "ISTAT_Code_LAU")
  # Convert to matrx
  X_mat <- X_sp %>%
    dplyr::select(-ISTAT_Code_LAU) %>%
    st_drop_geometry() %>%
    as.matrix()
  
  # Coordinates of the centroids
  centr <- sf::st_centroid(x = Y_sp)
  coords <- sf::st_coordinates(x = centr)
  
  
  
  ##### Spatially clustered regression
  K <- Grid_models[m,4]
  knn <- knn2nb(knearneigh(centr, k=K))
  listw_Knn_dates <- nb2listw(knn)
  W <- as(as_dgRMatrix_listw(listw_Knn_dates), "CsparseMatrix")*K
  
  fit <- SCR(Y = Y_mat,
             X = X_mat,
             Sp = coords,
             offset = as.numeric(Offset_mat),
             W = W,
             G = Grid_models[m,3],
             family = Grid_models[m,5])
  
  return(list(fit_optim = fit, W_opt = W, K_opt = K, centr_opt = centr, coords_opt = coords,
              X_mat_opt = X_mat, Y_mat_opt = Y_mat, Offset_mat_opt = Offset_mat,
              X_sp_opt = X_sp, Y_sp_opt = Y_sp, Offset_sp_opt = Offset_sp,
              shape_red_opt = shape_red))
}
