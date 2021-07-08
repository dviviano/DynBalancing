
#' Estimate dynamic average treatment effects
#'
#' This function computes the ATE from a panel.
#' Each row of the panel must contain a unit observed at a given date.
#' The function estimates the ATE using linear projections for potential outcomes
#' and balancing covariates dynamically. Reference: Viviano, Davide and Bradic Jelena. 2021. Dynamic covariate balancing: estimating treatment effect over time.
#'
#' @param panel_data data.frame format. It must be a panel with each row corresponding to a unit at a given date.
#' @param covariates_names vector of strings with column names of the covariates to include in the regression.
#' @param Time_name string denoting the column name of the panel corresponding to the date.
#' @param unit_name string denoting the column name of the panel corresponding to the individual.
#' @param outcome_name string denoting the column name of the panel corresponding to the outcome variable. Rows with missing outcomes will be removed.
#' @param treatment_name name of the column of the panel corresponding to the treatment assignment. Rows with missing treatment assignments will be removed.
#' @param ds1 counterfactual treatment history of interest of the first potential outcome.
#' @param ds2 counterfactual treatment history of interest of the second potential outcome.
#' @param fixed_effect vector of strings containing the columns to include in the regression as fixed effects.
#' @param pooled boolean indicating whether to run a pooled regression. If pooled = TRUE, the last T periods will not be used, with T equal to the length of the treatment ds1.
#' @param param See documentation.
#' @return summaries: matrix with ATE, variance, critical quantile, expected value of each potential outcome and corresponding critical quantile;
#' @return imbalances_summaries: summaries of imbalances over different periods;
#' @return plots: list of plots. These contain the imbalance plot, ATE and conditional expectation plots;
#' @return all_results: raw results from the optimization;
#' @export
#'
DynBalancing_ATE  <- function(panel_data,
                      covariates_names,
                      Time_name, unit_name, outcome_name,
                      treatment_name, ds1, ds2,
                      fixed_effects = NA,
                      pooled = F, cluster_SE = NA,
                      params = list()){

  params_default$variables_to_plot = covariates_names
  params_default$lags = length(ds1)
  params = exctract_params(params, params_default = params_default)
  alpha = params$alpha
  variables_to_plot = params$variables_to_plot
  lb = params$lb
  ub = params$ub
  params_gurobi = params$params_gurobi
  debias = params$debias
  adaptive_balancing = params$adaptive_balancing
  pooled_all = params$pooled_all
  method = params$method
  conservative_quantile = params$robust_quantile
  open_source = params$open_source
  continuous_treatment = params$continuous_treatment
  instrument_name = params$instrument_name
  final_period = params$final_period
  regularization = params$regularization
  boot = params$boot
  balancing = params$balancing
  comparison_plot  = params$comparison_plot
  with_sparsity_plot = params$with_sparsity_plot
  comparison_plot = params$comparison_plot
  nfolds = params$nfolds
  histogram_plot = params$histogram_plot
  cluster_SE = params$cluster_SE
  fast_adaptive = params$fast_adaptive
  grid_length  = params$grid_length
  n_beta_nonsparse = params$n_beta_nonsparse
  ratio_coefficients = params$ratio_coefficients
  lags = params$lags
  demeaned_fe = params$demeaned_fe
  if(continuous_treatment){
    conservative_quantile = F
    warning('Continuous treatment selected. No balancing performed.')
  }


  if(regularization & continuous_treatment){
    stop('You cannot regularized using a continuous treatment. This has not been implemented yet.')
  }
  if(all(variables_to_plot %in% covariates_names) == F){
    stop('variables_to_plot contains names outside the covariates_names passed.')
  }
  if(open_source == F){
    tryCatch(library(gurobi), error= function(e) stop('Gurobi is not installed. Install Gurobi or select open_source = T to use an open source software.'))
    warning('Selected commercial software by default (Gurobi). You can choose the built-in open source software choosein open_source = T')
    }

  if(lb > ub){
    stop('Lower bound is greater than upper bound')
  }
  if(method %in% c('lasso_plain', 'lasso_subsample') == F){
    stop('Method should be either lasso_plain or lasso_subsample. You can remove the penalization by selecting regularization = F.')
  }

  if(is.na(params_gurobi) == F){
    warning('Params is not the default for the quadratic program. Errors may occur if it is not a list. See help(gurobi) for details.')
  }

  if(alpha > 0.1){
    warning('Size smaller than 0.9 selected.')
  }

  if(balancing %in% c('DCB', 'AIPW', 'IPW', 'IPW_MSM') == F){
    stop('Select for balancing one among DCB, AIPW, IPW, IPW_MSM (default DCB).')
  }

  if(is.null(names(panel_data))){
    stop('Column names are missing in panel data.')
  }

  if(pooled & is.na(fixed_effects[1]) == F){
    ## Adjust time fixed effects when pooling the data to reflect the data shaping
    ## used for pooled regression
    if(Time_name %in% fixed_effects){
      tt <- which(fixed_effects == Time_name)
      fixed_effects[tt] = 'new_Time'
    }
  }

  ## cluster also at the unit level for pooled regression
  if(pooled & is.na(cluster_SE)){
    cluster_SE = unit_name
    ## if cluster is at a different level then regression is automatically adjusted also at the unit level
  }

  Time_name <- which(names(panel_data) == Time_name)
  unit_name <- which(names(panel_data) == unit_name)
  outcome_name <- which(names(panel_data) == outcome_name)
  treatment_name <- which(names(panel_data) == treatment_name)
  if(is.na(final_period)) final_period = max(panel_data[,Time_name], na.rm = T)
  if(is.null(params$initial_period)){
  initial_period = final_period - length(ds1) + 1
  } else{
    initial_period = params$initial_period
    ## If initial period we only pool starting from that period
    pooled_all = F
  }

  if(length(Time_name) == 0){
    stop('Wrong Time_name passed')
  }

  if(length(unit_name) == 0){
    stop('Wrong unit_name passed')
  }

  if(length(outcome_name) == 0){
    stop('Wrong outcome_name passed')
  }

  if(length(treatment_name) == 0){
    stop('Wrong treatment_name passed')
  }

  if(final_period %in% panel_data[, Time_name] == F){
    stop('Wrong final_period passed')
  }


  if(initial_period %in% panel_data[, Time_name] == F){
    stop('You asked to evaluate a treatment history that is too long with respect to the number of periods in your data. Try to change ds1, ds2.')
  }

  if(initial_period >= final_period){
    warning('ds only contains one element. No dynamics will be considered.')
  }



  panel_data <- panel_data[!is.na(panel_data[, treatment_name]),]

  ## Instrument for continuous regressor
  Z = c(NA, NA)
  if(is.na(instrument_name) == F){
    ## First stage
    first_stage = predict(lm(panel_data[, treatment_name] ~. , data =  panel_data[, names(panel_data) %in% c(covariates_names,
                                                                                                             instrument_name,
                                                                                                             fixed_effects)]),
                          new_data = panel_data[, names(panel_data) %in% c(covariates_names, instrument_name,
                                                                                   fixed_effects)] )
    ZZ = rep(NA, nrow(panel_data))
    which_is_NA = apply(panel_data[, names(panel_data) %in% c(covariates_names, instrument_name,
                                                              fixed_effects)], 1, function(x) ifelse(sum(is.na(x)) > 0, 1, 0))
    ZZ[which(which_is_NA == 0)] = first_stage
    all_names = names(panel_data)
    panel_data = cbind(panel_data, treatment_new = ZZ)
    panel_data = as.data.frame(panel_data)
    treatment_name = dim(panel_data)[2]
    names(panel_data) = all_names
  }

  if(pooled){
  panel_data = create_pooled_matrix(panel_data, Time_name, unit_name, outcome_name,
                                  treatment_name,
                                  final_period, initial_period, pooled_all = pooled_all, length_treatment = length(ds1))
  unit_name <- which(names(panel_data) == 'new_name')
  }


  ## Remove units after the final period
  panel_data <- panel_data[panel_data[, Time_name] <= final_period,]
  ## Keep rows of interest for treatment history (note: pooling is keeping also past rows after replicate rows and assign different periods of exposures, see the data_shaping file)
  panel_data <- panel_data[panel_data[, Time_name] >= final_period - length(ds1) + 1,]
  panel_data = clean_data_from_missing_time(panel_data, unit_name, Time_name)
  periods <- panel_data[, Time_name]

  individuals <- panel_data[, unit_name]
  individuals_unique <- unique(individuals)

  if(length(unique(periods)) != length(ds1)) {
    stop('The number of periods is different from the length of ds1.')
  }

  if(length(unique(periods)) != length(ds2)) {
    stop('The number of periods is different from the length of ds2.')
  }

  ## Outcome
  Y_T <- panel_data[panel_data[, Time_name] == final_period, outcome_name]
  ## Handle missing values
  individuals_final_period <- individuals[panel_data[, Time_name] == final_period]
  missing_ind <- which(individuals_unique %in% individuals_final_period == FALSE)
  if(length(missing_ind) > 0){
    panel_data <- panel_data[individuals %in% missing_ind == F, ]
    warning(paste0('Missing outcome values for ', length(missing_ind), ' units. These have been removed.'))
    individuals <- panel_data[, unit_name]
    individuals_unique <- unique(individuals)
    periods <- panel_data[, Time_name]
  }



  ## Create treatment assignments
  D <- create_matrix_of_D(unit_name,  treatment_name, Time_name, panel_data)

  ## Check for missing values on assignments
  missings <- apply(D, 1, function(x) sum(is.na(x)))
  which_missings <- which(missings > 0)
  if(length(which_missings) > 0){
    unit_to_remove <- individuals_unique[which_missings]
    panel_data <- panel_data[as.character(individuals) %in% as.character(unit_to_remove) == F, ]
    warning(paste0('Missing treatment values for ', length(which_missings), ' units. These have been removed.'))
    individuals <- panel_data[, unit_name]
    individuals_unique <- unique(individuals)
    D <- D[-which_missings, ]
    if(length(ds1) == 1) D = matrix(D, ncol = 1)
    Y_T <- Y_T[-which_missings]
    periods <- panel_data[, Time_name]
  }


  ## Construct fixed effects in the regression
  ## Use fe only for the first regression and then project the estimated outcome net of fixed effects
  ## Here save the dimension of dim_fe
  dim_fe = 0
  if(is.na(fixed_effects[1]) == F){
    data_matrix = panel_data
    dim_saved = 0
    for(j in fixed_effects){
      columns_fe = which(names(panel_data)  == j)
      matrix_fe <- model.matrix(~ panel_data[, columns_fe] - 1)
      dim_saved = dim_saved + dim(matrix_fe)[2]
      data_matrix <- cbind(data_matrix, matrix_fe)
    }
    data_matrix = as.data.frame(data_matrix)
    names(data_matrix) = c(names(panel_data), paste0('FE_V', c(1:dim_saved)))
    covariates_names = c(covariates_names,  paste0('FE_V', c(1:dim_saved)))
    panel_data = data_matrix
    if(params$demeaned_fe) dim_fe = dim_saved
  }
 
  #plot_treatments <- create_plot_treatments(D, individuals_unique, time_periods = sort(unique(periods), decreasing = F))

  ## Create vector of idexes for cluster robust standard errors
  indexes_clustering = NA
  if(is.na(cluster_SE) == F){
    column_selected = which(names(panel_data) %in% cluster_SE)
    if(length(column_selected) > 1) {
      column_selected = column_selected[1]
      warning('Only the first entry of cluster_SE will be used. One way clustering only is implemented.')
    }
    my_columns = panel_data[panel_data[, Time_name] == final_period, column_selected]
    unique_elements = unique(my_columns)
    indexes_clustering = rep(NA, length(Y_T))
    k = 1
    for(j in unique_elements){
      indexes_clustering[my_columns == j] = k
      k = k + 1
    }
  }

  ## Create covariates
  ## Note: missing data with covariates are handled internally in DCB
  covariates <- create_matrix_of_covariates(unit_name, covariates_names, Time_name, panel_data, dim_fe)

  res <- compute_ATE(Covariates_t = covariates, Y_T, Ds = D, ds1 = ds1,
                     ds2 = ds2, Time = length(ds1), params = params_gurobi,
                     tolerance_constraint1 = 10**(-8),
                     lb = lb, ub = ub,  method = method, adaptive_balancing = adaptive_balancing,
                     debias = debias, penalization = regularization, balancing = balancing,
                     open_source = open_source, continuous_treatment = continuous_treatment,
                     Z = Z, indexes_clusters =  indexes_clustering, boot = boot, nfolds = nfolds,
                     fast_adaptive = fast_adaptive, grid_length  = grid_length,
                     n_beta_nonsparse = n_beta_nonsparse,
                     ratio_coefficients = ratio_coefficients, lags = lags, dim_fe = dim_fe)




  if(continuous_treatment == F & comparison_plot){
  est1_propensity <- compute_estimator_with_propensity(covariates, Y_T,
                                                       D, ds1, length(unique(periods)),
                                                       method = method, penalization = regularization, lags = lags)
  est2_propensity <- compute_estimator_with_propensity(covariates, Y_T,
                                                       D, ds2, length(unique(periods)),
                                                       method = method, penalization = regularization, lags = lags)


  plot_imbalance1 <- make_imbalance_plot_comparisons(res_estimator= res$estimators[[1]], est1_propensity,
                                         covariates, covariates_names,variables_to_plot, with_sparsity = with_sparsity_plot,
                                         histogram_plot = histogram_plot)
  plot_imbalance2 <- make_imbalance_plot_comparisons(res$estimators[[2]], est2_propensity, covariates, covariates_names, variables_to_plot, with_sparsity =with_sparsity_plot,
                                                     histogram_plot = histogram_plot)

  if(balancing != 'DCB'){
    plot_imbalance1 = NA
    plot_imbalance2 = NA
    warning('Imbalance plot available if selecting balancing as default (DCB).')
  }
  }

  if(continuous_treatment == F & comparison_plot == F){


    plot_imbalance1 <- make_imbalance_plot(res_estimator= res$estimators[[1]],
                                           covariates, covariates_names,variables_to_plot, with_sparsity = with_sparsity_plot,
                                           histogram_plot = histogram_plot)
    plot_imbalance2 <- make_imbalance_plot(res$estimators[[2]], covariates, covariates_names, variables_to_plot, with_sparsity = with_sparsity_plot,
                                           histogram_plot = histogram_plot)

  }
  if(continuous_treatment){
   plot_imbalance1 = plot_imbalance2 = list()
  }

  plot_coefficient_1 <- make_coefficients_plot(res$estimators[[1]], covariates, covariates_names, variables_to_plot, histogram_plot = histogram_plot)
  plot_coefficient_2 <- make_coefficients_plot(res$estimators[[2]], covariates, covariates_names, variables_to_plot, histogram_plot = histogram_plot)

  summary_plot <- make_summary_plot(res, length(ds1), conservative_quantile, alpha)
  quantile1 = ifelse(conservative_quantile, sqrt(qchisq(1 - alpha,  length(ds1))), qnorm(1 - alpha/2))
  quantile2 = ifelse(conservative_quantile, sqrt(qchisq(1 - alpha,  2* length(ds1))), qnorm(1 - alpha/2))
  summaries <- c(res$ATE, sqrt(sum(res$variances)),
                 quantile2,  qnorm(1 - alpha/2), 
                 res$mu_hat[1], sqrt(res$variances[1]),
                 res$mu_hat[2], sqrt(res$variances[2]),
                 quantile1,  qnorm(1 - alpha/2))
  summaries <- as.data.frame(matrix(summaries, nrow = 1))
  names(summaries) <- c('ATE', 'SE_ATE', 'Robust_Quantile_ATE', 'Gaussian_Quantile_ATE', 
                        'Mu1', 'SE_mu1', 'Mu2', 'Variance_mu2',
                        'Robust_Quantile_mu',  'Gaussian_Quantile_ATE')
  return(list(summaries = summaries,
              imbalances_summaries = list(po_1 = plot_imbalance1$data_imbalance, po2 = plot_imbalance2$data_imbalance),
              plots = list(ATE = summary_plot,
                           imbalance1 = plot_imbalance1$plot,
                           imbalance2 = plot_imbalance2$plot,
                           coefficient1 = plot_coefficient_1,
                           coefficient2 = plot_coefficient_2),
              all_results = res ))
}



#' Estimate dynamic average treatment effects for different lags effects
#'
#' This function computes the ATE from a panel for a sequence of treatment histories.
#' Each row of the panel must contain a unit observed at a given date.
#' The function estimates the ATE using linear projections for potential outcomes
#' and balancing covariates dynamically. Reference: Viviano, Davide and Bradic Jelena. 2021. Dynamic covariate balancing: estimating treatment effect over time.
#'
#' @param panel_data data.frame format. It must be a panel with each row corresponding to a unit at a given date.
#' @param covariates_names vector of strings with column names of the covariates to include in the regression.
#' @param Time_name string denoting the column name of the panel corresponding to the date.
#' @param unit_name string denoting the column name of the panel corresponding to the individual.
#' @param outcome_name string denoting the column name of the panel corresponding to the outcome variable. Rows with missing outcomes will be removed.
#' @param treatment_name name of the column of the panel corresponding to the treatment assignment. Rows with missing treatment assignments will be removed.
#' @param ds1 counterfactual treatment history of interest of the first potential outcome.
#' @param ds2 counterfactual treatment history of interest of the second potential outcome.
#' @param histories_length numeric vector containing the number of lags for which the effect must be computed
#' @param fixed_effect vector of strings containing the columns to include in the regression as fixed effects.
#' @param pooled boolean indicating whether to run a pooled regression. If pooled = TRUE, the last T periods will not be used, with T equal to the length of the treatment ds1.
#' @param param See documentation.
#' @return plots: plots of the main results
#' @return all_results:  matrix with ATE, variance, critical quantile, expected value of each potential outcome and corresponding critical quantile;
#' @export

DynBalancing_History <- function(panel_data,
                                       covariates_names,
                                       Time_name, unit_name, outcome_name,
                                       treatment_name, ds1, ds2,
                                       histories_length, fixed_effects = NA,
                                       pooled = F, cluster_SE = NA,
                                       params = list( )){


  params_default$variables_to_plot = covariates_names
  params = exctract_params(params, params_default = params_default)
  alpha = params$alpha
  variables_to_plot = params_default$variables_to_plot
  lb = params$lb
  ub = params$ub
  params_gurobi = params$params_gurobi
  debias = params$debias
  adaptive_balancing = params$adaptive_balancing
  pooled_all = params$pooled_all
  method = params$method
  conservative_quantile = params$robust_quantile
  open_source = params$open_source
  continuous_treatment = params$continuous_treatment
  instrument_name = params$instrument_name
  final_period = params$final_period
  regularization = params$regularization
  numcores = params$numcores
  balancing = params$balancing
  comparison_plot  = params$comparison_plot
  impulse_response = params$impulse_response
  ## Balancing also includes fixed effects
  if(continuous_treatment == F) params$demeaned_fe = F

  if(impulse_response){
    warning('impulse response selected. ds1, ds2 will be ignored.')
  }
  T_all = length(ds1)
  doParallel::registerDoParallel(numcores)

  store_results = foreach::foreach(i = sort(histories_length, decreasing = F), .combine = append)%dopar%{

  t_minus = T_all - i + 1
  dd1 = ds1[t_minus:T_all]
  dd2 = ds2[t_minus:T_all]
  if(impulse_response) {
    dd1 = c(1, rep(0, length(t_minus:T_all) - 1))
    dd2 = c(0, rep(0, length(t_minus:T_all) - 1))
  }
  list(DynBalancing_ATE(panel_data,
                        covariates_names,
                        Time_name, unit_name, outcome_name,
                        treatment_name, dd1, dd2,
                        fixed_effects = fixed_effects,
                        pooled = pooled, cluster_SE = cluster_SE,
                        params = params))
  }

  if(continuous_treatment){
    conservative_quantile = F
    warning('Continuous treatment selected. No balancing performed.')
  }
  ATEs <- sapply(store_results, function(x) x$all_results$ATE)
  Variances1 <- sapply(store_results, function(x) x$all_results$variances[1])
  Variances2 <- sapply(store_results, function(x) x$all_results$variances[2])
  mu1 <- sapply(store_results, function(x) x$all_results$mu_hat[1])
  mu2 <- sapply(store_results, function(x) x$all_results$mu_hat[2])
  matrix_of_results <- cbind(ATEs, Variances1 + Variances2,
                             sapply(sort(histories_length, decreasing = F), function(x) sqrt(qchisq(1 - alpha, 2 * (x)))),
                             mu1, Variances1, mu2, Variances2,
                             sapply(sort(histories_length, decreasing = F), function(x) sqrt(qchisq(1 - alpha,  (x)))),
                             sort(histories_length, decreasing = F)
                             )
  if(conservative_quantile == F){
    matrix_of_results <- cbind(ATEs, sqrt(Variances1 + Variances2),
                               sapply(sort(histories_length, decreasing = F), function(x) qnorm(1- alpha/2)),
                               mu1, sqrt(Variances1), mu2, sqrt(Variances2),
                               sapply(sort(histories_length, decreasing = F), function(x) 1 - alpha/2),
                               sort(histories_length, decreasing = F)
    )

  }
  matrix_of_results <- as.data.frame(matrix_of_results)
  names(matrix_of_results) <- c('ATE', 'SE_ATE', 'Quantile_ATE',
                                'mu1', 'SE1', 'mu2',
                                'SE2', 'Quantile_mu',
                                'Period_length')

  time_all <-  sapply(sort(histories_length, decreasing = F), function(x) (x))
  if(balancing %in% c('DCB', "AIPW")){
  final_results <- compute_ATE_plot_back_time(store_results,
                                              sort(histories_length, decreasing = F),
                                              final_period = final_period, conservative_quantile, alpha)
  } else {
    final_results = NA
  }
  return(list(plots = final_results,
              all_results = matrix_of_results))
  }


#' Estimate dynamic average treatment effects for different years, capturing heterogeneity in time.
#'
#' This function computes the ATE from a panel for a sequence of different years.
#' Each row of the panel must contain a unit observed at a given date.
#' The function estimates the ATE using linear projections for potential outcomes
#' and balancing covariates dynamically. Reference: Viviano, Davide and Bradic Jelena. 2021. Dynamic covariate balancing: estimating treatment effect over time.
#'
#' @param panel_data data.frame format. It must be a panel with each row corresponding to a unit at a given date.
#' @param covariates_names vector of strings with column names of the covariates to include in the regression.
#' @param Time_name string denoting the column name of the panel corresponding to the date.
#' @param unit_name string denoting the column name of the panel corresponding to the individual.
#' @param outcome_name string denoting the column name of the panel corresponding to the outcome variable. Rows with missing outcomes will be removed.
#' @param treatment_name name of the column of the panel corresponding to the treatment assignment. Rows with missing treatment assignments will be removed.
#' @param ds1 counterfactual treatment history of interest of the first potential outcome.
#' @param ds2 counterfactual treatment history of interest of the second potential outcome.
#' @param final_periods numeric vector containing the years for which the effect needs to be computed.
#' @param fixed_effect vector of strings containing the columns to include in the regression as fixed effects.
#' @param pooled boolean indicating whether to run a pooled regression. If pooled = TRUE, the first T periods will not be used, with T equal to the length of the treatment ds1.
#' @param param See documentation.
#' @return plots: plots of the main results
#' @return all_results:  matrix with ATE, variance, critical quantile, expected value of each potential outcome and corresponding critical quantile;
#' @export

DynBalncing_Het_ATE <- function(panel_data,
                                         covariates_names,
                                         Time_name, unit_name, outcome_name,
                                         treatment_name, ds1, ds2,
                                         final_periods,  fixed_effects = NA,
                                         pooled = F, cluster_SE = NA,
                                         params = list( )){

  params_default$variables_to_plot = covariates_names
  params = exctract_params(params, params_default = params_default)
  alpha = params$alpha
  variables_to_plot = params_default$variables_to_plot
  lb = params$lb
  ub = params$ub
  params_gurobi = params$params_gurobi
  debias = params$debias
  adaptive_balancing = params$adaptive_balancing
  pooled_all = params$pooled_all
  method = params$method
  conservative_quantile = params$robust_quantile
  open_source = params$open_source
  continuous_treatment = params$continuous_treatment
  instrument_name = params$instrument_name
  final_period = params$final_period
  regularization = params$regularization
  numcores = params$numcores
  balancing = params$balancing
  comparison_plot  = params$comparison_plot
  doParallel::registerDoParallel(numcores)

  store_results = foreach::foreach(t = sort(final_periods, decreasing = F ), .combine = append)%dopar%{
    params$final_period = t
    list(DynBalancing_ATE(panel_data,
                          covariates_names,
                          Time_name, unit_name, outcome_name,
                          treatment_name, ds1, ds2,
                          fixed_effects = fixed_effects,
                          pooled = pooled, cluster_SE = cluster_SE,
                          params = params))
  }

  final_results <- compute_ATE_plot_forward_time(store_results, sort(final_periods, decreasing = F ),
                                                 length(ds1))

  ATEs <- sapply(store_results, function(x) x$all_results$ATE)
  Variances1 <- sapply(store_results, function(x) x$all_results$variances[1])
  Variances2 <- sapply(store_results, function(x) x$all_results$variances[2])
  mu1 <- sapply(store_results, function(x) x$all_results$mu_hat[1])
  mu2 <- sapply(store_results, function(x) x$all_results$mu_hat[2])
  matrix_of_results <- cbind(ATEs, sqrt(Variances1 + Variances2),
                             sapply(rep(length(ds1), length(final_period)),
                                    function(x) sqrt(qchisq(1 - alpha, 2 * x))),
                             mu1, sqrt(Variances1), mu2, sqrt(Variances2),
                             sapply(rep(length(ds1), length(final_period)), function(x) sqrt(qchisq(1 - alpha, x))),
                             final_periods
  )
  matrix_of_results <- as.data.frame(matrix_of_results)
  names(matrix_of_results) <- c('ATE', 'SE_ATE', 'Conservative_Quantile_ATE',
                                'mu1', 'SE1', 'SE2',
                                'Variance2', 'Conservative_quantile_mu',
                                'Period_length')

  return(list(plots = final_results,
              all_results = matrix_of_results))

}
