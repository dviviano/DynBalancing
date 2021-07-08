
## MAin function for coefficients

compute_coefficients <- function(Time, Y_T, Ds, Covariates_t, ds, method, penalization, continuous_treatment,
                                 vector_for_predictions = NA, nfolds = 10, lags = Time, dim_fe = 0){


  pred_t <- list()
  coef_t <- list()
  Covariates_t_nonna <- list()
  ## Deal with missing values
  not_nas <- list()
  if(method == 'lasso_subsample'){


    subsample <- sapply(c(1:Time),  function(x) apply(as.matrix(Ds[,1:x],ncol = x), 1, function(y) all(y == ds[1:x])))
    nas_Y <- which(is.na(Y_T))
    predictions <- Y_T
    for(t in rev(1:Time)){
      XX_t <- Covariates_t[[t]]
      predictions[nas_Y] <- NA
      ## Remove missing values
      all_matrix <- cbind(XX_t, Y_final = predictions)
      not_nas[[t]] <- which(apply(all_matrix, 1, function(x) sum(is.na(x)) == 0))
      all_matrix <- na.omit(all_matrix)
      Covariates_t_nonna[[t]] <- as.matrix(all_matrix[, -dim(all_matrix)[2]])
      ## Estimate on the sub-sample incurring the treatment history of interest
      if(penalization){
        reg <- glmnet::cv.glmnet(y = all_matrix[subsample[ not_nas[[t]],t], dim(all_matrix)[2]],
                                 x = as.matrix(all_matrix[subsample[ not_nas[[t]],t], -dim(all_matrix)[2]]), nfolds = nfolds)
      } else {
        reg <- glmnet::glmnet(y = all_matrix[subsample[ not_nas[[t]],t], dim(all_matrix)[2]], x = as.matrix(all_matrix[subsample[ not_nas[[t]],t], -dim(all_matrix)[2]]),
                      lambda = exp(-8), alpha = 0)
      }
      predictions <- rep(NA, length(Y_T))
      ## Remove fixed effects from the predictions if dim_fe > 0 (active only for continuous treatments)
      if(dim_fe > 0 & t == Time) XX_t[, c((dim(XX_t)[2] - dim_fe + 1):dim(XX_t)[2])] = 0
      predictions[which(apply(XX_t, 1, function(x) sum(is.na(x)) == 0))] <-  predict(reg, newx = as.matrix(na.omit(XX_t)))
      coef_t[[t]] <- coef(reg)
      pred_t[[t]] <-   predict(reg, newx = as.matrix(all_matrix[, -dim(all_matrix)[2]]) )
    }
  }  else if(method == 'lasso_plain'){

    predictions <- Y_T
    nas_Y <- which(is.na(Y_T))
    model_effect = list()
    for(t in rev(1:Time)){
      XX_t <- Covariates_t[[t]]
      predictions[nas_Y] <- NA
      all_matrix <- cbind(XX_t, Y_final = predictions)

      ## Deal with missing values
      not_nas[[t]] <- which(apply(all_matrix, 1, function(x) sum(is.na(x)) == 0))
      all_matrix <- na.omit(all_matrix)
      Covariates_t_nonna[[t]] <- as.matrix(all_matrix[, -dim(all_matrix)[2]])
      if(penalization){
        coefficients_not_to_penalize = min(t, lags)
        reg <- glmnet::cv.glmnet(y = all_matrix[, dim(all_matrix)[2]], x = cbind(all_matrix[, -dim(all_matrix)[2]], Ds[not_nas[[t]],1:t]), nfolds = nfolds
                                 , penalty.factor = c(rep(1, dim(XX_t)[2] + t - coefficients_not_to_penalize),
                                                      rep(0, coefficients_not_to_penalize)))
      } else {
        reg <- glmnet::glmnet(y = all_matrix[, dim(all_matrix)[2]], x = cbind(all_matrix[, -dim(all_matrix)[2]], Ds[not_nas[[t]],1:t]), lambda = exp(-8),
                      alpha = 0)
      }
      ## Remove fixed effects from the predictions if dim_fe > 0 (active only for continuous treatments)
      if(dim_fe > 0 & t == Time) XX_t[, c((dim(XX_t)[2] - dim_fe + 1):dim(XX_t)[2])] = 0
      if(t > 1) {
        matrix_for_predictions <- cbind(XX_t, Ds[, 1:(t-1)], ds[t])
      } else {
        dd = rep(ds[1], dim(XX_t)[1])
        ## THis is used for a continuous treatment only
        if(all(sapply(vector_for_predictions, is.na)) == F){
          dd = vector_for_predictions
        }
        matrix_for_predictions <- cbind(XX_t, dd)
      }
      predictions <- rep(NA, length(predictions))
      predictions[which(apply(matrix_for_predictions, 1, function(x) sum(is.na(x)) == 0))] <- predict(reg, newx = na.omit(matrix_for_predictions) )
      ## last_coef saved for continuous treatment use only
      model_effect[[t]] = coef(reg)[length(coef(reg))]
      coef_t[[t]] <- coef(reg)[-c((dim(XX_t)[2] + 2):length(coef(reg)))]
      pred_t[[t]] <-   predict(reg, newx = cbind(all_matrix[, -dim(all_matrix)[2]], t(matrix(ds[1:t], nrow = t, ncol = dim(all_matrix)[1])) ))
    }
  }
  return(list(coef_t, pred_t, Covariates_t_nonna, not_nas, model_effect))
}

