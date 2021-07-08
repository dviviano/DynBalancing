
## Gamma in first Period (open source function)

## X1: Pass X in first period
## D1: Treatment assignment in first period
## d_1:treatment of interest in first period
## params: for the Quadratic Program formulation in Gurobi
## K_1: constant in front of the balancing constraint
## tolerance constraint: minimum tolerance to avoid negative weights (weights greater or qual than zero)
## with_beta: boolena: whether covariates with non-zero coefficieints must be have stricter balancing constraints (in such case beta_hat argument must be the vector of coefficients)
## beta_hat: vector of coefficients obtained from the estimated linear model (only if with_beta = T)
## K2 : constant for those parameters with smaller (zero) coefficient
## n_beta_nonsparse: small number indicating the tuning under which a (scale-free)  parameter is considered almost zero (e.g., exp(-4))
## ratio_coefficients: if there are too many non sparse coefficients (and we are in high-dim, i.e., length(coef) >= 100)
##                     and more than 1 - ratio_coefficients are non sparse impose stronger balancing on the largest 1/3rd of the coefficients

## Return: vector of gammas in the first period and coefficients

compute_gamma1_os <- function(X1, D1, d_1, K_1, tolerance_constraint1 = 10**(-8),
                           with_beta = F, beta_hat = NA, K_2 = 0, n_beta_nonsparse = 10**(-4),
                           ratio_coefficients = 1/3){


  ## History for ones
  ## Rescale coefficients to make them comparable
  beta_hat[-1] = beta_hat[-1] * apply(X1, 2, function(x) sqrt(var(x, na.rm= T)))
  X1_1 <- X1[D1 == d_1,]
  Xbar_1 <- apply(X1,2 ,mean)
  p <- dim(X1_1)[2]
  n <- dim(X1_1)[1]

  Q <- matrix(0, nrow = n, ncol = n)
  diag(Q) <- 1

  constraint_matrix <- matrix(0, nrow = 2*p, ncol= n )

  ## Constraints over each variable
  for(k in 1:p){
    constraint_matrix[k,] <- -X1_1[,k]
  }
  for(k in (p + 1):(2*p)){
    constraint_matrix[k,] <- X1_1[,k - p]
  }
  ## Constraint on sum of gammas being one
  constraint_matrix <- rbind(rep(1, n), constraint_matrix)


  b_vector <- c(1, -(Xbar_1) - K_1*sqrt(log(p)/sqrt(n)), (Xbar_1) - K_1*sqrt(log(p)/sqrt(n)))

  if(with_beta){
    non_zerocoef1 <- which(abs(beta_hat) > n_beta_nonsparse) - 1 ## Remove the intercept
    non_zerocoef1 <- non_zerocoef1[-1]
    if(length(beta_hat) >= 90 & sum(beta_hat == 0) < (1 - ratio_coefficients) * length(beta_hat)){
      non_zerocoef1 = sort(abs(beta_hat[-1]), decreasing = TRUE, index.return = TRUE)$ix[1:floor(ratio_coefficients * length(beta_hat))]
    }
    b_vector <- c(1, -(Xbar_1) - K_2*sqrt(log(p)/(n)), (Xbar_1) - K_2*sqrt(log(p)/(n)))
    if(length(non_zerocoef1) > 0){
    b_vector[c(non_zerocoef1 + 1, non_zerocoef1 + p + 1)] <- c((-Xbar_1[non_zerocoef1]) - K_1*sqrt(log(p)/((n))), (Xbar_1[non_zerocoef1]) - K_1*sqrt(log(p)/((n))))
    }
  }



  lower_bound = 0
  b_vector <- c(b_vector, rep(0, n), rep(- log(n) * n**(-2/3), n))
  constraint_matrix <- rbind(constraint_matrix, diag(1, n), diag(-1, n))
  ## Constant vector
  c <- rep(0, n)
  my_sol <- quadprog::solve.QP(Q, c, t(constraint_matrix), b_vector, meq = 1)
  gamma_1 <- rep(0, dim(X1)[1])
  gamma_1[D1 == d_1] <- my_sol$solution
  return(gamma_1)

}



## Gamma in the second period (open source function)

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


compute_gammat_os <- function(gammat_before, XX, d_t, D,  K_1,  tolerance_constraint1 = 10**(-8),
                           with_beta = F, beta_hat = NA, K_2 = 0, n_beta_nonsparse = 10**(-4),
                           ratio_coefficients  = 1/3){

  ## Deal with one against multiple periods
  if(is.null(dim(gammat_before))){
    gammat_minus1 <- gammat_before
  } else{
    gammat_minus1 <- gammat_before[, dim(gammat_before)[2]] }

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

  constraint_matrix <- matrix(0, nrow = 2*p, ncol= n )

  ## Constraints over each variable
  for(k in 1:p){
    constraint_matrix[k,] <- -XX_1[,k]
  }
  for(k in (p + 1):(2*p)){
    constraint_matrix[k,] <- XX_1[,k - p]
  }
  ## Constraint on sum of gammas being one
  constraint_matrix <- rbind(rep(1, n), constraint_matrix)
  constraint_matrix <- rbind(constraint_matrix, diag(1, n), diag(-1, n))
  b_vector <- c(1, -(Xbar_1) - K_1*sqrt(log(p)/sqrt(n)), (Xbar_1) - K_1*sqrt(log(p)/sqrt(n)))

  if(with_beta){
    non_zerocoef2 <- which(abs(beta_hat) > n_beta_nonsparse) - 1 ## Remove the intercept
    non_zerocoef2 <- non_zerocoef2[-1]
    if(length(beta_hat) >= 90 & sum(beta_hat == 0) < (1 - ratio_coefficients) * length(beta_hat)){
      non_zerocoef2 = sort(abs(beta_hat[-1]), decreasing = TRUE, index.return = TRUE)$ix[1:floor(ratio_coefficients * length(beta_hat))]
    }
    b_vector <- c(1, -(Xbar_1) - K_2*sqrt(log(p)/(sqrt(n))), (Xbar_1) - K_2*sqrt(log(p)/(sqrt(n))))
    if(length(non_zerocoef2) > 0){
    b_vector[c(1 + non_zerocoef2, non_zerocoef2 + p + 1)] <- c((-Xbar_1[non_zerocoef2]) - K_1*sqrt(log(p)/(sqrt(n))), (Xbar_1[non_zerocoef2]) - K_1*sqrt(log(p)/(sqrt(n))))
    }
  }




  b_vector <- c(b_vector, rep(0, n), rep(- log(n) * n**(-2/3), n))
  ## Constant vector
  c <- rep(0, n)
  my_sol <- quadprog::solve.QP(Q, c, t(constraint_matrix), b_vector, meq = 1)
  gamma_t <- rep(0, dim(XX)[1])
  gamma_t[subsample] <- my_sol$solution

  return(gamma_t)

}
