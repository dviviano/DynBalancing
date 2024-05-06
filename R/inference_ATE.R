
## Compute variance

compute_variance <- function(gammas, predictions, not_nas, Y_T){

  Time <- dim(predictions)[2]
  epsilon_hat <- Y_T[not_nas[[Time]]] - predictions[not_nas[[Time]],Time]
  ## Weights are weights/sum(weights) (stabilized) for IPW
  gammas_times_epsilon_squared <- (gammas[not_nas[[Time]],Time]**2)%*%(epsilon_hat**2)
  if(Time == 1) return( gammas_times_epsilon_squared)
  Var <- rep(0, Time - 1)
  k <- 1
  for(x in 2:Time){
    ind1 <- not_nas[[x]]
    ind2 <- not_nas[[x-1]]
    ind_int <- intersect(ind1, ind2)
    pp1 <- predictions[ind_int, x]
    pp2 <- predictions[ind_int, x-1]
    gam <- gammas[ind_int,x-1]
    Var[k] <- (gam**2)%*%((pp1 - pp2)**2)
    k <- k + 1
    }
  Var <- sum(Var) + gammas_times_epsilon_squared
  return(Var)
}


## Variance cluster robust

help_clustering = function(x, index_T, gg_T, epsilon_hat){
  sum(gg_T[index_T == x] * epsilon_hat[index_T == x])**2
}
compute_variance_cluster_robust <- function(gammas, predictions, not_nas, Y_T, indexes_clusters){

  indexes_clusters_unique = unique(indexes_clusters)
  Time <- dim(predictions)[2]
  epsilon_hat <- Y_T[not_nas[[Time]]] - predictions[not_nas[[Time]],Time]
  gg_T = gammas[not_nas[[Time]],Time]
  index_T = indexes_clusters[not_nas[[Time]]]
  gammas_times_epsilon_squared = sum(sapply(indexes_clusters_unique, function(x) help_clustering(x, index_T, gg_T, epsilon_hat )))
  if(Time == 1) return(gammas_times_epsilon_squared)
  Var <- rep(0, Time - 1)
  k <- 1
  for(x in 2:Time){
    ind1 <- not_nas[[x]]
    ind2 <- not_nas[[x-1]]
    ind_int <- intersect(ind1, ind2)
    index_T = indexes_clusters[ind_int]
    pp1 <- predictions[ind_int, x]
    pp2 <- predictions[ind_int, x-1]
    gam <- gammas[ind_int,x-1]
    Var[k] <- sum(sapply(indexes_clusters_unique, function(x) help_clustering(x, index_T, gam, pp1 - pp2 )))
    k <- k + 1
  }
  Var <- sum(Var) + gammas_times_epsilon_squared
  return(Var)
}



## Covariates_t: list where each element is the matrix of covariates and outcomes of the past periods
## Y_T: outcome at endline period
## Ds: matrix of treatment assignments, each column corresponds to a different period
## ds: history of treatment of interest
## Time: how many time periods
## params: params to be passed to Gurobi (NA fix default parameters)
## tolerance_constraint1: lower bounds on gammas not to be negative
## lb: lower bound on the constants to be used in the balancing equation
## ub: upper bound on the constants to be used in the balancing equation
## Methods: lasso_plain, lasso_subsample
## balancing: either DCB, IPW, AIPW or hdCBPS
## open_source: use an open source software?
## continuous treatment: use a contd treatment, Z instrument for continuous treatment
## indexes_clusters = indexes for cluster robust variance estimation
## boot: bootstap samples for contd treat
## nfolds: nfolds for cv
## fast_adaptive: whether use a fast adaptive choice of the tuning parameters
# grid_length: length size for the choice of the tuning parameters
## n_beta_nonsparse: small number to indicate whether an estimated coefficient can be considered equal to zero
## ratio_coefficients: if we are in high-dimensions, and less than 1 - ratio_coefficients (e.g., 2/3rd) of the coefficients are sparse,
##                     impose stricter constants for balancing on the largest 1/3rd (idea: it imposes more stringent balancing on most relevant variables)
## (see Algorithm C.1 in the Appendix )
## dim_fe > 0 only if demeaned = T

compute_ATE <- function(Covariates_t, Y_T, Ds, ds1,ds2, Time, params = NA, tolerance_constraint1 = 10**(-8),
                        lb = 0.0001, ub = 10,  method = 'lasso_subsample', adaptive_balancing = T, debias = F,
                        penalization= T, balancing = 'DCB', open_source = F, continuous_treatment = F,
                        Z = NA, indexes_clusters = NA, boot = 50, nfolds = 10, fast_adaptive = F,
                        grid_length = grid_length,
                        n_beta_nonsparse = 10**(-4),
                        ratio_coefficients = 1/3, lags = Time, dim_fe = 0){

  if(continuous_treatment){

    estimator1 <-  compute_estimator_cont(Covariates_t, Y_T, Ds, ds = ds1, Time, params, tolerance_constraint1,
                                     lb, ub,  method, adaptive_balancing, debias = debias,
                                     penalization = penalization, open_source = open_source,
                                     Z = Z,
                                     indexes_clusters = indexes_clusters, boot = boot, nfolds = nfolds, lags = lags,
                                     dim_fe = dim_fe)

    estimator2 <- compute_estimator_cont(Covariates_t, Y_T, Ds, ds2, Time, params, tolerance_constraint1,
                                    lb, ub,  method, adaptive_balancing, debias = debias,
                                    penalization = penalization, open_source = open_source,
                                    Z = Z,
                                    indexes_clusters = indexes_clusters, boot = boot, nfolds = nfolds, lags = lags,
                                    dim_fe = dim_fe)

  }
  if(balancing == 'DCB' & continuous_treatment == F){
    if(fast_adaptive){
  estimator1 <-  compute_estimator_fast_adaptive(Covariates_t, Y_T, Ds, ds = ds1, Time, params, tolerance_constraint1,
                                   lb, ub,  method, adaptive_balancing, debias = debias,
                                   penalization = penalization, open_source = open_source, nfolds = nfolds,
                                   length_grid = grid_length, n_beta_nonsparse = n_beta_nonsparse,
                                   ratio_coefficients = ratio_coefficients, lags = lags)

  estimator2 <- compute_estimator_fast_adaptive(Covariates_t, Y_T, Ds, ds2, Time, params, tolerance_constraint1,
                                  lb, ub,  method, adaptive_balancing, debias = debias,
                                  penalization = penalization, open_source = open_source, nfolds = nfolds,
                                  length_grid = grid_length, n_beta_nonsparse = n_beta_nonsparse,
                                  ratio_coefficients = ratio_coefficients, lags = lags)
    } else {
      estimator1 <-  compute_estimator(Covariates_t, Y_T, Ds, ds = ds1, Time, params, tolerance_constraint1,
                                                     lb, ub,  method, adaptive_balancing, debias = debias,
                                                     penalization = penalization, open_source = open_source, nfolds = nfolds,
                                       length_grid = grid_length, n_beta_nonsparse = n_beta_nonsparse,
                                       ratio_coefficients = ratio_coefficients, lags = lags)

      estimator2 <- compute_estimator(Covariates_t, Y_T, Ds, ds2, Time, params, tolerance_constraint1,
                                                    lb, ub,  method, adaptive_balancing, debias = debias,
                                                    penalization = penalization, open_source = open_source, nfolds = nfolds,
                                      length_grid = grid_length, n_beta_nonsparse = n_beta_nonsparse,
                                      ratio_coefficients = ratio_coefficients, lags = lags)
    }
  }
  if (balancing %in% c('AIPW', 'IPW', 'IPW_MSM') & continuous_treatment == F) {

    estimator1 <-  compute_estimator_with_propensity(Covariates_t, Y_T, Ds, ds1, Time, method,
                                     penalization = penalization, type = balancing, nfolds = nfolds, lags = lags)
    estimator2 <- compute_estimator_with_propensity(Covariates_t, Y_T, Ds, ds2, Time, method,
                                                    penalization = penalization, type = balancing, nfolds = nfolds, lags = lags)

    if(balancing != 'AIPW') warning('Balancing selected has not DR adjustment. Size distortion may occur.')
  }
  bias1 <- ifelse(is.na(estimator1$bias)[1], 0, estimator1$bias)
  bias2 <- ifelse(is.na(estimator2$bias)[1], 0, estimator2$bias)
  ATE <- estimator1$mu_hat - estimator2$mu_hat - bias1 + bias2
  if(continuous_treatment){
    Variance1 = estimator1$variance
    Variance2 = estimator2$variance
    return(list(ATE = ATE, mu_hat = c(estimator1$mu_hat, estimator2$mu_hat), variances = c(Variance1, Variance2),
                estimators = list(estimator1, estimator2)))

  }

  if(is.na(indexes_clusters[1])){
  Variance1 <- compute_variance(gammas = estimator1$gammas, predictions = estimator1$predictions,
                                not_nas = estimator1$not_nas , Y_T =  Y_T)
  Variance2 <- compute_variance(estimator2$gammas, estimator2$predictions, estimator2$not_nas, Y_T)
  } else {
    Variance1 <- compute_variance_cluster_robust (gammas = estimator1$gammas, predictions = estimator1$predictions,
                                  not_nas = estimator1$not_nas , Y_T =  Y_T, indexes_clusters)
    Variance2 <- compute_variance_cluster_robust (estimator2$gammas, estimator2$predictions, estimator2$not_nas, Y_T,
                                                  indexes_clusters)

  }
  return(list(ATE = ATE, mu_hat = c(estimator1$mu_hat, estimator2$mu_hat), variances = c(Variance1, Variance2),
              estimators = list(estimator1, estimator2)))
}



