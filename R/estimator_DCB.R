


## Function to compute mu_hat(d_1:T)

## Covariates_t: list where each element is the matrix of covariates and outcomes of the past periods used in the regression (Covariates_[[t]] corresponds to H_{i,t})
## Y_T: outcome at endline period
## Ds: matrix of treatment assignments, each column corresponds to a different period
## ds: history of treatment of interest
## Time: how many time periods
## params: params to be passed to Gurobi (NA fix default parameters)
## tolerance_constraint1: lower bounds on gammas not to be negative
## lb: lower bound on the constants to be used in the balancing equation
## ub: upper bound on the constants to be used in the balancing equation
## Methods: lasso_plain, lasso_subsample
## Adaptive balancing: impose weaker moment conditions on inactive variables
## open_source: use estimation of gamma with quad.prog instead of gurobi?
## continuous treatment: use a contd treatment, Z instrument for continuous treatment
## return: final estimator and matrix of gamma_hat
## n_beta_nonsparse: constant to indicate whether a certain estimated coefficient can be considered equal to zero
## ratio_coefficients: it prioritize the largest ratio coefficients variables for balancing if most of the variables are non zero
## lags: number of carryovers to include in the regressions

## Faster version is compute_estimator_faster
compute_estimator <- function(Covariates_t, Y_T, Ds, ds, Time, params = NA, tolerance_constraint1 = 10**(-8),
                              lb = 0.0001, ub = 10,  method = 'lasso_subsample',
                              adaptive_balancing = T, debias = F,
                              penalization = T, open_source = F, nfolds = 10,
                              length_grid = 1000, n_beta_nonsparse = 10**(-4),
                              ratio_coefficients = 1/3, lags = Time){


  component_mu <- rep(NA, Time)
  n <- dim(Ds)[1]
  coefs = compute_coefficients(Time, Y_T, Ds, Covariates_t, ds, method, penalization, nfolds = nfolds, lags = lags)
  pred_t <- coefs[[2]]
  coef_t <- coefs[[1]]
  Covariates_t_nonna = coefs[[3]]
  not_nas = coefs[[4]]
  sequence_bounds <- seq(from = ub/3, to = ub, length = 3)
  K_first <- seq(from = lb, to = sequence_bounds[1], length = floor(length_grid**(1/3)))
  K_second <- seq(from = sequence_bounds[1], to = sequence_bounds[2], length = floor(length_grid**(1/3)))
  K_third <- seq(from = sequence_bounds[2], to = sequence_bounds[3], length = floor(length_grid**(1/3)))
  all_Ks <- list(K_first, K_second, K_third)


  ## Estimation of weights in the first period
  ## Set smaller constraints on non-active variables in high dimensions
  with_beta = F
  beta_hat <- NA
  if(adaptive_balancing){

    beta_hat <- coef_t[[1]]
    with_beta = T
  }
  ## Loop over each possible K_1 and pick the feasible solutions with smallest constant
  for (JJ in all_Ks){
    for(k1 in JJ){
      K2 = k1
      if(adaptive_balancing) K2 <- JJ
      for(k2 in K2){
    if(open_source){
      gammas_first <- tryCatch(compute_gamma1_os(X1 = Covariates_t_nonna[[1]],
                                                 D1 = Ds[not_nas[[1]],1], d_1 = ds[1],  K_1 = k1, tolerance_constraint1 = tolerance_constraint1,
                                                 with_beta = with_beta, beta_hat = beta_hat, K_2 = k2,
                                                 n_beta_nonsparse = n_beta_nonsparse, ratio_coefficients = ratio_coefficients), error = function(e){NA})
    } else {

      gammas_first <- tryCatch(compute_gamma1(Covariates_t_nonna[[1]], Ds[not_nas[[1]],1], d_1 = ds[1], params =params, K_1 = k1, tolerance_constraint1 = tolerance_constraint1,
                                              with_beta = with_beta, beta_hat = beta_hat, K_2 = k2,
                                              n_beta_nonsparse = n_beta_nonsparse, ratio_coefficients = ratio_coefficients), error = function(e){NA})

    }
    suppressWarnings(if(is.na(gammas_first) == F){break})
      }
      suppressWarnings(if(is.na(gammas_first) == F){break})
    }
    suppressWarnings(if(is.na(gammas_first) == F){break})
  }

  suppressWarnings(if(is.na(gammas_first)){  stop('Error: infeasible problem for period 1. Try to increase the ub parameter.') })
  previous_component <- 1/length(not_nas[[1]])
  component_mu[1] <- (gammas_first - previous_component)%*%pred_t[[1]]
  ## Sequential estimation
  keep_gammas <- matrix(0, nrow = n, ncol = Time) ## Assign zero weights to missing values
  keep_gammas[not_nas[[1]],1] <- gammas_first
  if(Time > 1){
    for(t in 2:Time){

      ## Repeat for subsequent periods
      with_beta = F
      beta_hat <- NA
      if(adaptive_balancing){
        beta_hat <- coef_t[[t]]
        with_beta = T
      }
      ## Loop over each possible K_1 and pick the feasible solutions with smallest constant
      for (JJ in all_Ks){
        for(k1 in JJ){
          K2 = k1
          if(adaptive_balancing) K2 <- JJ
          for(k2 in K2){
        if(open_source){
          gammas_tj <- tryCatch(compute_gammat_os(gammat_before = keep_gammas[not_nas[[t]],t-1],
                                                  XX = Covariates_t_nonna[[t]],  d_t = ds[1:t],
                                                  D = Ds[not_nas[[t]],1:t],  K_1 = k1, tolerance_constraint1 = tolerance_constraint1,
                                                  with_beta = with_beta, beta_hat = beta_hat, K_2 = k2,
                                                  n_beta_nonsparse = n_beta_nonsparse, ratio_coefficients = ratio_coefficients), error = function(e){NA})
        } else{
          gammas_tj <- tryCatch(compute_gammat(keep_gammas[not_nas[[t]],t-1], Covariates_t_nonna[[t]],  d_t = ds[1:t], Ds[not_nas[[t]],1:t], params =params, K_1 = k1, tolerance_constraint1 = tolerance_constraint1,
                                               with_beta = with_beta, beta_hat = beta_hat, K_2 = k2,
                                               n_beta_nonsparse = n_beta_nonsparse, ratio_coefficients = ratio_coefficients), error = function(e){NA})

        }
        suppressWarnings(if(is.na(gammas_tj) == F){break})
        }
        suppressWarnings(if(is.na(gammas_tj) == F){break})
      }
      suppressWarnings(if(is.na(gammas_tj) == F){break})
    }


      suppressWarnings(if(is.na(gammas_tj)){  stop(paste0('Error: infeasible problem for period ', t,  '. Try to increase the ub parameter.')) })
      keep_gammas[not_nas[[t]],t] <- gammas_tj
      previous_component <- keep_gammas[not_nas[[t]],t-1]/sum(keep_gammas[not_nas[[t]],t-1])
      component_mu[t] <- (keep_gammas[not_nas[[t]],t] - keep_gammas[not_nas[[t]],t-1])%*%pred_t[[t]]

    }
  }
  final_predictions <- matrix(0, nrow = n, ncol = Time)
  for(t in 1:Time){
    final_predictions[not_nas[[t]],t] <- pred_t[[t]]
  }
  mu_hat <- keep_gammas[not_nas[[Time]],Time]%*%Y_T[not_nas[[Time]]] - sum(component_mu)

  final_bias <- NA

  ## Consider a debias component in the final expression
  if(debias){
    coef_B <- coef_t
    for(i in 1:20){
      indexes <- sample(c(1:n), replace = T, size = n)
      coef_bb <- compute_coefficients_boot(Time, Y_T[indexes], Ds[indexes,], Covariates_t, ds, method, indexes, penalization, nfolds = nfolds, lags = lags)[[1]]
      for(t in 1:Time){
        ## Add ups coefficients
        coef_B[[t]] <- as.vector(coef_B[[t]]) + as.vector(coef_bb[[t]])
      }
    }

    ## Take average over coefficients and debias
    bias <- rep(NA, Time)
    coef_forbias <- coef_t[[1]] - coef_B[[1]]/20
    coef_forbias <- as.vector(coef_forbias)
    coef_forbias <- coef_forbias[-1]

    bias[1] <- coef_forbias%*%as.vector((keep_gammas[not_nas[[1]],1] - 1/n)%*%Covariates_t_nonna[[1]])
    for(t in 2:Time){
      coef_forbias <- coef_t[[t]] - coef_B[[t]]/20
      coef_forbias <- as.vector(coef_forbias)
      coef_forbias <- coef_forbias[-1]

      bias[t] <- coef_forbias%*%as.vector((keep_gammas[not_nas[[t]],t] - keep_gammas[not_nas[[t]],t-1])%*%Covariates_t_nonna[[t]])
    }
    final_bias <- sum(bias)
  }

  return(list(mu_hat = mu_hat , gammas = keep_gammas, predictions = final_predictions, not_nas = not_nas, coef = coef_t, bias = final_bias))
}

