############# Compute Gammas


## Gamma in first Period

## X1: Pass X in first period
## D1: Treatment assignment in first period
## d_1:treatment of interest in first period
## params: for the Quadratic Program formulation in Gurobi
## K_1: constant in front of the balancing constraint
## tolerance constraint: minimum tolerance to avoid negative weights (weights greater or qual than zero)
## with_beta: boolena: whether covariates with non-zero coefficieints must be have stricter balancing constraints (in such case beta_hat argument must be the vector of coefficients)
## beta_hat: vector of coefficients obtained from the estimated linear model (only if with_beta = T)
## K2 : constant for those parameters with smaller (zero) coefficient
## n_beta_nonsparse: small number indicating the tuning under which a (scale-free) parameter is considered almost zero (e.g., exp(-4))
## ratio_coefficients: if there are too many non sparse coefficients (and we are in high-dim, i.e., length(coef) >= 100)
##                     and more than 1 - ratio_coefficients are non sparse impose stronger balancing on the largest ratio_coefficients of the coefficients
## Return: vector of gammas in the first period and coefficients

compute_gamma1 <- function(X1, D1, d_1, params = NA, K_1, tolerance_constraint1 = 10**(-8),
                           with_beta = F, beta_hat = NA, K_2 = 0, n_beta_nonsparse = 10**(-4),
                           ratio_coefficients = 1/3){

  ## Set parameters for Gurobi

  if(is.na(params)){
    params <- list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit =2000, Threads = 1, Heuristics=0,Cuts=0)
  }
  ## History for ones
  ## Rescale coefficients to make them comparable
  beta_hat[-1] = beta_hat[-1] * apply(X1, 2, function(x) sqrt(var(x, na.rm= T)))

  X1_1 <- X1[D1 == d_1,]
  Xbar_1 <- apply(X1,2 ,mean)
  p <- dim(X1_1)[2]
  n <- dim(X1_1)[1]

  Q <- matrix(0, nrow = n, ncol = n)
  diag(Q) <- 1 #nu_hat_1**2
  model <- list()
  ## If beta_hat are passed imposed constraints only on the positive coefficients of beta_hat



  constraint_matrix <- matrix(0, nrow = 2*p + 1, ncol= n )

  ## Constraints over each variable
  for(k in 1:p){
    constraint_matrix[k,] <- X1_1[,k]
  }
  for(k in (p + 1):(2*p)){
    constraint_matrix[k,] <- X1_1[,k - p]
  }
  ## Constraint on sum of gammas being one
  constraint_matrix[2*p + 1,] <- rep(1, n)


  b_vector <- c((Xbar_1) + K_1*sqrt(log(p)/(n)), (Xbar_1) - K_1*sqrt(log(p)/(n)),
                1)

  if(with_beta){
    non_zerocoef1 <- which(abs(beta_hat) > n_beta_nonsparse) - 1 ## Remove the intercept
    non_zerocoef1 <- non_zerocoef1[-1]
    ## If too few coefficients are sparse and we are in high-dim give priority to the largest ones for balancing
    if(length(beta_hat) >= 90 & sum(beta_hat == 0) < (1  - ratio_coefficients) * length(beta_hat)){
      non_zerocoef1 = sort(abs(beta_hat[-1]), decreasing = TRUE, index.return = TRUE)$ix[1:floor(ratio_coefficients * length(beta_hat))]
    }
    b_vector <- c((Xbar_1) + K_2*sqrt(log(p)/(n)), (Xbar_1) - K_2*sqrt(log(p)/(n)),
                  1)
    b_vector[c(non_zerocoef1, non_zerocoef1 + p)] <- c((Xbar_1[non_zerocoef1]) + K_1*sqrt(log(p)/((n))), (Xbar_1[non_zerocoef1]) - K_1*sqrt(log(p)/((n))))

  }




  model <- list()

  model$sense <- c(rep('<=', p), rep('>=', p), '=')

  ## Constant vector
  c <- rep(0, n)
  model$A <- constraint_matrix
  model$obj <- c
  model$rhs <- b_vector
  model$lb <- rep(0 + tolerance_constraint1, n)
  model$up <- rep(log(n) * (n**(-2/3)), n)
  model$vtype <- rep('C', n)
  model$modelsense <- 'min'
  model$Q <- Q
  results <- gurobi::gurobi(model = model, params = params)
  gamma_1 <- rep(0, dim(X1)[1])
  gamma_1[D1 == d_1] <- results$x
  return(gamma_1)

}

## Gamma in the second periood

## gammat_minus1: vector of gammas computed in period t - 1
## XX: all covariates, and past outcomes used in the regression
## D: matrix, each column corresponds to the treatment assignment in corresponding period period
## d_t : treatment history of interest
## params: for the Quadratic Program formulation in Gurobi
## K_1: constant in front of the balancing constraint
## tolerance constraint: minimum tolerance to avoid negative weights (weights greater or qual than zero)
## with_beta: boolena: whether covariates with non-zero coefficieints must be have stricter balancing constraints (in such case beta_hat argument must be the vector of coefficients)
## beta_hat: vector of coefficients obtained from the estimated linear model (only if with_beta = T)
## K2 : constant for those parameters with smaller (zero) coefficient
## n_beta_nonsparse: select n variables with the largest coefficient for more stringent balancing if too many variables are non-zero

## Return: matrix of gammas where each column correspond to each period and last column correspond to period t (estimated gamma)

compute_gammat <- function(gammat_before, XX, d_t, D, params = NA, K_1,  tolerance_constraint1 = 10**(-8),
                           with_beta = F, beta_hat = NA, K_2 = 0, n_beta_nonsparse = 10,
                           ratio_coefficients = 1/3){

  ## Deal with one against multiple periods
  if(is.null(dim(gammat_before))){
    gammat_minus1 <- gammat_before
    } else{
  gammat_minus1 <- gammat_before[, dim(gammat_before)[2]] }

  if(is.na(params)){
    params <- list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit =2000, Threads = 1, Heuristics=0,Cuts=0)
  }
  subsample <- apply(D, 1, function(x) all(x == d_t))
  nn <- dim(XX)[1]
  ## Rescale coefficients to make them comparable
  beta_hat[-1] = beta_hat[-1] * apply(XX, 2, function(x) sqrt(var(x, na.rm= T)))
  XX_1 <- XX[subsample,]
  Xbar_1 <- apply(XX, 2, function(x) gammat_minus1%*%x)
  Xbar_1 <- Xbar_1/sum(gammat_minus1) ## Rescale in case of missing values
  p <- dim(XX_1)[2]
  n <- dim(XX_1)[1]

  Q <- matrix(0, nrow = n, ncol = n)
  diag(Q) <- 1

  constraint_matrix <- matrix(0, nrow = 2*p + 1, ncol= n )

  ## Constraints over each variable
  for(k in 1:p){
    constraint_matrix[k,] <- XX_1[,k]
  }
  for(k in (p + 1):(2*p)){
    constraint_matrix[k,] <- XX_1[,k - p]
  }
  ## Constraint on sum of gammas being one
  constraint_matrix[2*p + 1,] <- rep(1, n)

  b_vector <- c((Xbar_1) + K_1*sqrt(log(p)/sqrt(n)), (Xbar_1) - K_1*sqrt(log(p)/sqrt(n)),
                1)

  if(with_beta){
    non_zerocoef2 <- which(abs(beta_hat) > n_beta_nonsparse) - 1 ## Remove the intercept
    non_zerocoef2 <- non_zerocoef2[-1]
    ## If too few coefficients are sparse and we are in high-dim give priority to the largest ones for balancing
    if(length(beta_hat) >= 90 & sum(beta_hat == 0) < (1 - ratio_coefficients) * length(beta_hat)){
      non_zerocoef2 = sort(abs(beta_hat[-1]), decreasing = TRUE, index.return = TRUE)$ix[1:floor(ratio_coefficients * length(beta_hat))]
    }
    b_vector <- c((Xbar_1) + K_2*sqrt(log(p)/(sqrt(n))), (Xbar_1) - K_2*sqrt(log(p)/(sqrt(n))),
                  1)
    b_vector[c(non_zerocoef2, non_zerocoef2 + p)] <- c((Xbar_1[non_zerocoef2]) + K_1*sqrt(log(p)/(sqrt(n))), (Xbar_1[non_zerocoef2]) - K_1*sqrt(log(p)/(sqrt(n))))

  }



  ## Constant vector
  c <- rep(0, n)
  model <- list()
  model$A <- constraint_matrix
  model$obj <- c
  model$sense <- c(rep('<=', p), rep('>=', p), '=')
  model$rhs <- b_vector
  model$lb <- rep(0 + tolerance_constraint1, n)
  model$up <- rep(log(n) * (n**(-2/3)), n)
  model$vtype <- rep('C', n)
  model$modelsense <- 'min'
  model$Q <- Q


  results2 <- gurobi::gurobi(model = model, params = params)
  gamma_t <- rep(0, dim(XX)[1])
  gamma_t[subsample] <- results2$x

  return(gamma_t)

}


