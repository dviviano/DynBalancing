### Plot function: return the imbalance plot for IPW weights and our weights
## Input: results of DCB, IPW results, covariates (list) and names of covariates
## with_sparsity: if T it also reports which coefficients are zero

## Make plots with comparisons among methods
make_imbalance_plot_comparisons <- function(res_estimator, IPW_results, covariates, covariates_names,  variables_to_plot,
                                with_sparsity = F, histogram_plot){

  all_time <- length(covariates)
  plot_var = which(covariates_names %in% variables_to_plot )

  imbalance_t1 <- apply(covariates[[1]][, plot_var], 2, function(x) sum((res_estimator$gammas[,1] - 1/(length(x) - sum(is.na(x))))*x/sd(x,na.rm = T), na.rm = T))
  imbalance_t1_IPW <- apply(covariates[[1]][, plot_var], 2, function(x) sum((IPW_results$gammas[,1] - 1/(length(x) -sum(is.na(x))))*x/sd(x,na.rm = T), na.rm = T))


  types <- sapply(res_estimator$coef[[1]][c(1 + plot_var)], function(x) ifelse(x == 0, 'DCB - s', 'DCB - ns'))

  data_to_plot <- cbind(c(log(abs(imbalance_t1) + 1), log(abs(imbalance_t1_IPW) + 1)),
                        rep(1, 2 * dim(covariates[[1]][, plot_var])[2]),
                        c(rep('DCB', dim(covariates[[1]][, plot_var])[2]),
                          rep('IPW', dim(covariates[[1]][, plot_var])[2])
                        ),
                        c(rep(variables_to_plot, 2)
                        ),
                        c(types, rep('IPW', dim(covariates[[1]][, plot_var])[2]))
  )
  if(all_time > 1){
  for(t in 2:all_time){

    imbalance_t <- apply(covariates[[t]][, plot_var], 2, function(x) sum((res_estimator$gammas[,t] - res_estimator$gammas[,t - 1])*x/sd(x,na.rm = T), na.rm = T))
    imbalance_t_IPW <- apply(covariates[[t]][, plot_var], 2, function(x) sum((IPW_results$gammas[,t] - IPW_results$gammas[,t - 1])*x/sd(x,na.rm = T), na.rm = T))
    types <- sapply(res_estimator$coef[[t]][c(plot_var + 1)], function(x) ifelse(x == 0, 'DCB - s', 'DCB - ns'))


    data_to_plot <- rbind(data_to_plot, cbind(c(log(abs(imbalance_t) + 1), log(abs(imbalance_t_IPW) + 1)),
                                              c(rep(t, 2 * dim(covariates[[t]][, plot_var])[2])),
                                              c(rep('DCB', dim(covariates[[t]][, plot_var])[2]),
                                                rep('IPW', dim(covariates[[t]][, plot_var])[2]))
                                              ,
                                              c(rep(variables_to_plot, 2)
                                              ),
                                              c(types, rep('IPW', dim(covariates[[t]][, plot_var])[2]))))

  }
  }
  data_to_plot <- as.data.frame(data_to_plot)
  data_to_plot[,1] <- as.numeric(as.character(data_to_plot[,1]))
  names(data_to_plot) <- c('LogImbalance', 'Period', 'Method', 'Covariates', 'Type')
  data_to_plot[,3] <- as.character(data_to_plot[,3])
  data_to_plot[,4] <- as.character(data_to_plot[,4])
  data_to_plot[,2] <- factor(as.numeric(as.character(data_to_plot[,2])),
                             levels = sort(unique(as.numeric(as.character(data_to_plot[,2]))), decreasing = F))
  minimum <- min(data_to_plot[,1])
  maximum <- max(data_to_plot[,1])
  if(length(variables_to_plot) < 10 & histogram_plot == F){
  plot1 <- ggplot2::ggplot(data_to_plot, ggplot2::aes(LogImbalance, Covariates)) +
    ggplot2::geom_point(ggplot2::aes(shape  = Method, col = Method), size = 4) +
   # geom_segment(ggplot2::aes(x = minimu,
  #                   y = Covariates, xend = LogImbalance,
  #                   yend = Covariates), size = 2, color = 'red', data = data_to_plot[data_to_plot$Method != 'IPW',]) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~Period,  scales = "free_x") +
    ggplot2::theme(axis.title.x=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
          axis.title.y=ggplot2::element_text(size = 0),
          plot.title = ggplot2::element_text(size=22),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 20),
          legend.position="top")

  # +
  #scale_shape_manual(values=c(13, 2))
  if(with_sparsity){
    plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(LogImbalance, Covariates), shape = 13,
                                data = data_to_plot[data_to_plot$Type == 'DCB - s', ], size = 3.5)
  }
  } else {
    plot1 <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = LogImbalance, color = Method, fill = Method)) +
      ggplot2::geom_histogram(position="identity", alpha=0.5) +
      # geom_segment(ggplot2::aes(x = minimu,
      #                   y = Covariates, xend = LogImbalance,
      #                   yend = Covariates), size = 2, color = 'red', data = data_to_plot[data_to_plot$Method != 'IPW',]) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~Period, scales = "free_x") +
      ggplot2::theme(axis.title.x=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
                     axis.title.y=ggplot2::element_text(size = 0),
                     plot.title = ggplot2::element_text(size=22),
                     axis.text.x = ggplot2::element_text(size = 20),
                     axis.text.y = ggplot2::element_text(size = 20),
                     legend.position="top")
}
  return(list(plot = plot1, data_imbalance = data_to_plot))
}

## Make plots only for the method considered
make_imbalance_plot <- function(res_estimator, covariates, covariates_names,  variables_to_plot,
                                            with_sparsity = F, histogram_plot){

  all_time <- length(covariates)
  plot_var = which(covariates_names %in% variables_to_plot )

  imbalance_t1 <- apply(covariates[[1]][, plot_var], 2, function(x) sum((res_estimator$gammas[,1] - 1/(length(x) - sum(is.na(x))))*x/sd(x,na.rm = T), na.rm = T))


  types <- sapply(res_estimator$coef[[1]][c(1 + plot_var)], function(x) ifelse(x == 0, 'DCB - s', 'DCB - ns'))

  data_to_plot <- cbind(log(abs(imbalance_t1) + 1),
                        rep(1, dim(covariates[[1]][, plot_var])[2]),
                        c(rep('DCB', dim(covariates[[1]][, plot_var])[2])
                        ), variables_to_plot,
                        types)

  if(all_time > 1){
    for(t in 2:all_time){

      imbalance_t <- apply(covariates[[t]][, plot_var], 2, function(x) sum((res_estimator$gammas[,t] - res_estimator$gammas[,t - 1])*x/sd(x,na.rm = T), na.rm = T))
      types <- sapply(res_estimator$coef[[t]][c(plot_var + 1)], function(x) ifelse(x == 0, 'DCB - s', 'DCB - ns'))


      data_to_plot <- rbind(data_to_plot,cbind(log(abs(imbalance_t) + 1),
                                               rep(t, dim(covariates[[t]][, plot_var])[2]),
                                               c(rep('DCB', dim(covariates[[t]][, plot_var])[2])
                                               ), variables_to_plot,
                                               types))

    }
  }
  data_to_plot <- as.data.frame(data_to_plot)
  data_to_plot[,1] <- as.numeric(as.character(data_to_plot[,1]))
  names(data_to_plot) <- c('LogImbalance', 'Period', 'Method', 'Covariates', 'Type')
  data_to_plot[,3] <- as.character(data_to_plot[,3])
  data_to_plot[,4] <- as.character(data_to_plot[,4])
  data_to_plot[,2] <- factor(as.numeric(as.character(data_to_plot[,2])), levels = sort(unique(as.numeric(as.character(data_to_plot[,2]))), decreasing = F))
  minimum <- min(data_to_plot[,1])
  maximum <- max(data_to_plot[,1])
  if(length(variables_to_plot) < 10 & histogram_plot == F){
  plot1 <- ggplot2::ggplot(data_to_plot, ggplot2::aes(LogImbalance, Covariates)) +
    ggplot2::geom_point( size = 4, shape = 1) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~Period, scales = "free_x") +
    ggplot2::theme(axis.title.x=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
          axis.title.y=ggplot2::element_text(size = 0),
          plot.title = ggplot2::element_text(size=22),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 20),
          legend.position="top")
  # +
  #scale_shape_manual(values=c(13, 2))
  if(with_sparsity){
    plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(LogImbalance, Covariates), shape = 13,
                                data = data_to_plot[data_to_plot$Type == 'DCB - s', ], size = 3.5)
  }
  } else {
    plot1 <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = LogImbalance)) +
      ggplot2::geom_histogram(position="identity", alpha=0.5) +
      # geom_segment(ggplot2::aes(x = minimu,
      #                   y = Covariates, xend = LogImbalance,
      #                   yend = Covariates), size = 2, color = 'red', data = data_to_plot[data_to_plot$Method != 'IPW',]) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~Period, scales = "free_x") +
      ggplot2::theme(axis.title.x=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
                     axis.title.y=ggplot2::element_text(size = 0),
                     plot.title = ggplot2::element_text(size=22),
                     axis.text.x = ggplot2::element_text(size = 20),
                     axis.text.y = ggplot2::element_text(size = 20),
                     legend.position="top")
  }
  return(list(plot = plot1, data_imbalance = data_to_plot))
}



## Return the coefficient values
## Same inputs as before

make_coefficients_plot <- function(res_estimator, covariates, covariates_names, variables_to_plot, histogram_plot){

  ## Coefficients plot
  plot_var = which(covariates_names %in% variables_to_plot )
  coef1 <- res_estimator$coef[[1]]
  coef1 = coef1[-1]
  coef1 = coef1[plot_var]
  data_to_plot <- cbind(c(coef1),
                        c(rep(1, dim(covariates[[1]][, plot_var])[2])),
                        c(rep(variables_to_plot, 1)
                        ))
  if(length(covariates) > 1){
  for(t in 2:length(covariates)){
    coef1 <- res_estimator$coef[[t]]
    coef1 = coef1[-1]
    coef1 = coef1[plot_var]
    data_to_plot <- rbind(data_to_plot, cbind(c(coef1),
                                              c(rep(t, dim(covariates[[t]][, plot_var])[2])),
                                              c(rep(variables_to_plot, 1)
                                              )))

  }
  }
  data_to_plot <- as.data.frame(data_to_plot)
  data_to_plot[,1] <- as.numeric(as.character(data_to_plot[,1]))
  data_to_plot[,3] <- as.character(data_to_plot[,3])
  names(data_to_plot) <- c('Coefficient', 'Period', 'Covariates')
  data_to_plot[,2] <- factor(as.numeric(as.character(data_to_plot[,2])), levels = sort(unique(as.numeric(as.character(data_to_plot[,2]))), decreasing = F))
  if(length(variables_to_plot) < 10 & histogram_plot == F){
  plot1 <- ggplot2::ggplot(data_to_plot, ggplot2::aes(Coefficient, Covariates)) +
    ggplot2::geom_point( size = 4, shape = 'triangle') +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~Period) +
    ggplot2::theme(axis.title=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
          plot.title = ggplot2::element_text(size=22),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 20))
  } else {
    plot1 <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = Coefficient)) +
      ggplot2::geom_histogram(position="identity", alpha=0.5) +
      # geom_segment(ggplot2::aes(x = minimu,
      #                   y = Covariates, xend = LogImbalance,
      #                   yend = Covariates), size = 2, color = 'red', data = data_to_plot[data_to_plot$Method != 'IPW',]) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~Period) +
      ggplot2::theme(axis.title.x=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
                     axis.title.y=ggplot2::element_text(size = 0),
                     plot.title = ggplot2::element_text(size=22),
                     axis.text.x = ggplot2::element_text(size = 20),
                     axis.text.y = ggplot2::element_text(size = 20),
                     legend.position="top")
    }

  return(plot1)
}

## It plots the ATE and the mu_hat with se
## Inputs: results of DCB and the length of the period
## dark gay is Gaussian standard errors and light gray is adjusted standard errors

make_summary_plot <- function(results, num_periods, conservative_quantile, alpha){

  quantile1 = ifelse(conservative_quantile, sqrt(qchisq(1 - alpha,  num_periods)), qnorm(1 - alpha/2))
  quantile2 = ifelse(conservative_quantile, sqrt(qchisq(1 - alpha,  2* num_periods)), qnorm(1 - alpha/2))
  ATE <- results$ATE
  mu1 <- results$mu_hat[1]
  mu2 <- results$mu_hat[2]
  var1 <- results$variances[1]
  var2 <- results$variances[2]

  data_to_plot <- cbind(c(ATE, mu1, mu2), c('ATE', 'mu1', 'mu2'))
  data_to_plot <- as.data.frame(data_to_plot)
  data_to_plot[,1] <- as.numeric(as.character(data_to_plot[,1]))
  names(data_to_plot) <- c('Results', 'Type')
  plot1 <- ggplot2::ggplot(data_to_plot[data_to_plot$Type != 'ATE', ], ggplot2::aes(Results, Type)) +
    ggplot2::geom_segment(ggplot2::aes(x = mu1 - quantile1  * sqrt(var1),
                     y = 'mu1', xend = mu1 + quantile1 * sqrt(var1),
                     yend = 'mu1'), size = 2, color = 'lightgray') +
    ggplot2::geom_segment(ggplot2::aes(x = mu1 - 1.64 * sqrt(var1), y = 'mu1', xend = mu1 + 1.64 * sqrt(var1),
                     yend = 'mu1'), size = 2, color = 'darkgray') +
    ggplot2::geom_segment(ggplot2::aes(x = mu2 - quantile1 * sqrt(var2),
                     y = 'mu2', xend = mu2 + quantile2 * sqrt(var2),
                     yend = 'mu2'), size = 2, color = 'lightgray') +
    ggplot2::geom_segment(ggplot2::aes(x = mu2 - 1.64 * sqrt(var2), y = 'mu2', xend = mu2 + 1.64 * sqrt(var2),
                     yend = 'mu2'), size = 2, color = 'darkgray') +
    ggplot2::geom_point( size = 4, shape = 'triangle') +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title=ggplot2::element_text(size = 0), legend.text=ggplot2::element_text(size = 25),
          plot.title = ggplot2::element_text(size=22),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 20))


  plot2 <- ggplot2::ggplot(data_to_plot[data_to_plot$Type == 'ATE', ], ggplot2::aes(Results, Type)) +
    ggplot2::geom_segment(ggplot2::aes(x = ATE - quantile2 * sqrt(var1 + var2),
                     y = 'ATE', xend = ATE + quantile2 * sqrt(var1 + var2),
                     yend = 'ATE'), size = 2, color = 'lightgray') +
    ggplot2::geom_segment(ggplot2::aes(x = ATE - 1.64 * sqrt(var1 + var2), y = 'ATE', xend = ATE + 1.64 * sqrt(var1 + var2),
                     yend = 'ATE'), size = 2, color = 'darkgray') +
    ggplot2::geom_point( size = 4, shape = 'triangle') +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title=ggplot2::element_text(size = 0), legend.text=ggplot2::element_text(size = 25),
          plot.title = ggplot2::element_text(size=22),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 20))
  final_plot <- list(mus = plot1, ATE = plot2)
  return(final_plot)
}

## Plot the ATE as we coonsider the same final period (Y_T), but we allow for
## longer and longer carry-over effects
## Inputs: store results: list of DCB results over each period
## periods used to go back in time
## the total time length

compute_ATE_plot_back_time <- function(store_results, periods, final_period, conservative_quantile, alpha){

  ATEs <- sapply(store_results, function(x) x$all_results$ATE )
  mu1 <- sapply(store_results, function(x) x$all_results$mu_hat[1] )
  mu2 <- sapply(store_results, function(x) x$all_results$mu_hat[2] )
  Var1 <- sapply(store_results, function(x) x$all_results$variances[1] )
  Var2 <- sapply(store_results, function(x) x$all_results$variances[2] )

  ATEs1 <- cbind(ATEs, Var1, Var2, periods )
  ATEs1 <- as.data.frame(ATEs1)
  names(ATEs1) <- c('ATE', 'Var1', 'Var2', 'Periods')
  ATEs1[,1] <- as.numeric(as.character(ATEs1[,1]))
  ATEs1[,2] <- as.numeric(as.character(ATEs1[,2]))
  ATEs1[,3] <- as.numeric(as.character(ATEs1[,3]))
  ATEs1[,4] <- factor(as.character(ATEs1[,4]), levels = sort(as.numeric(ATEs1[,4]), decreasing = F))

  if(conservative_quantile){
  plot1 <- ggplot2::ggplot(ATEs1, ggplot2::aes(Periods, ATE)) +
    ggplot2::geom_segment(ggplot2::aes(y = ATE - sqrt(qchisq(1 - alpha,  2 * as.numeric(as.character(Periods)))) * sqrt(Var1 + Var2),
                     x = Periods, yend = ATE + sqrt(qchisq(1 - alpha, 2 * as.numeric(as.character(Periods)))) * sqrt(Var1 + Var2),
                     xend = Periods), size = 2, color = 'lightgray') +
    ggplot2::geom_segment(ggplot2::aes(y = ATE - qnorm(1 - alpha/2) * sqrt(Var1 + Var2),
                     x = Periods, yend = ATE +  qnorm(1 - alpha/2) * sqrt(Var1 + Var2),
                     xend = Periods), size = 2, color = 'darkgray') +
    ggplot2::geom_point( size = 4, shape = 'triangle') +
    ggplot2::theme_bw()  +
    ggplot2::theme(axis.title=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
          plot.title = ggplot2::element_text(size=22),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 20))
  } else {
    plot1 <- ggplot2::ggplot(ATEs1, ggplot2::aes(Periods, ATE)) +
      ggplot2::geom_segment(ggplot2::aes(y = ATE -  qnorm(1 - alpha/2) * sqrt(Var1 + Var2),
                       x = Periods, yend = ATE +  qnorm(1 - alpha/2) * sqrt(Var1 + Var2),
                       xend = Periods), size = 2, color = 'darkgray') +
      ggplot2::geom_point( size = 4, shape = 'triangle') +
      ggplot2::theme_bw()  +
      ggplot2::theme(axis.title=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
            plot.title = ggplot2::element_text(size=22),
            axis.text.x = ggplot2::element_text(size = 20),
            axis.text.y = ggplot2::element_text(size = 20))

  }
  mus1 <- cbind(mu1, Var1, periods)
  mus1 <- as.data.frame(mus1)
  names(mus1) <- c('Mu1', 'Var1', 'Periods')
  mus1[,1] <- as.numeric(as.character(mus1[,1]))
  mus1[,2] <- as.numeric(as.character(mus1[,2]))
  mus1[,3] <- factor(as.character(mus1[,3]), levels = sort(as.numeric(mus1[,3]), decreasing = F))

  if(conservative_quantile){
  plot2 <- ggplot2::ggplot(mus1, ggplot2::aes(Periods, Mu1)) +
    ggplot2::geom_segment(ggplot2::aes(y = Mu1 - sqrt(qchisq(1 - alpha,  as.numeric(as.character(Periods)))) * sqrt(Var1),
                     x = Periods, yend = Mu1 + sqrt(qchisq(1 - alpha, as.numeric(as.character(Periods)))) * sqrt(Var1),
                     xend = Periods), size = 2, color = 'lightgray') +
    ggplot2::geom_segment(ggplot2::aes(y = Mu1 - qnorm(1 - alpha/2) * sqrt(Var1),
                     x = Periods, yend = Mu1 + qnorm(1 - alpha/2) * sqrt(Var1),
                     xend = Periods), size = 2, color = 'darkgray') +
    ggplot2::geom_point( size = 4, shape = 'triangle') +
    ggplot2::theme_bw()  +
    ggplot2::theme(axis.title=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
          plot.title = ggplot2::element_text(size=22),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 20))
  } else {
    plot2 <- ggplot2::ggplot(mus1, ggplot2::aes(Periods, Mu1)) +
      ggplot2::geom_segment(ggplot2::aes(y = Mu1 - qnorm(1- alpha/2) * sqrt(Var1),
                       x = Periods, yend = Mu1 + qnorm(1 - alpha/2) * sqrt(Var1),
                       xend = Periods), size = 2, color = 'darkgray') +
      ggplot2::geom_point( size = 4, shape = 'triangle') +
      ggplot2::theme_bw()  +
      ggplot2::theme(axis.title=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
            plot.title = ggplot2::element_text(size=22),
            axis.text.x = ggplot2::element_text(size = 20),
            axis.text.y = ggplot2::element_text(size = 20))
  }
  mus2 <- cbind(mu2, Var2, periods)
  mus2 <- as.data.frame(mus2)
  names(mus2) <- c('Mu2', 'Var2', 'Periods')
  mus2[,1] <- as.numeric(as.character(mus2[,1]))
  mus2[,2] <- as.numeric(as.character(mus2[,2]))
  mus2[,3] <- factor(as.character(mus2[,3]), levels = sort(as.numeric(mus2[,3]), decreasing = F))

  if(conservative_quantile){
  plot3 <- ggplot2::ggplot(mus2, ggplot2::aes(Periods, Mu2)) +
    ggplot2::geom_segment(ggplot2::aes(y = Mu2 - sqrt(qchisq(0.9,  as.numeric(as.character(Periods)))) * sqrt(Var2),
                     x = Periods, yend = Mu2 + sqrt(qchisq(0.9, as.numeric(as.character(Periods)))) * sqrt(Var2),
                     xend = Periods), size = 2, color = 'lightgray') +
    ggplot2::geom_segment(ggplot2::aes(y = Mu2 - 1.64 * sqrt(Var2),
                     x = Periods, yend = Mu2 + 1.64 * sqrt(Var2),
                     xend = Periods), size = 2, color = 'darkgray') +
    ggplot2::geom_point( size = 4, shape = 'triangle') +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
          plot.title = ggplot2::element_text(size=22),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 20))
  } else {
    plot3 <- ggplot2::ggplot(mus2, ggplot2::aes(Periods, Mu2)) +
      ggplot2::geom_segment(ggplot2::aes(y = Mu2 - qnorm(1 - alpha/2) * sqrt(Var2),
                       x = Periods, yend = Mu2 + qnorm(1 - alpha/2) * sqrt(Var2),
                       xend = Periods), size = 2, color = 'darkgray') +
      ggplot2::geom_point( size = 4, shape = 'triangle') +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
            plot.title = ggplot2::element_text(size=22),
            axis.text.x = ggplot2::element_text(size = 20),
            axis.text.y = ggplot2::element_text(size = 20))
  }
    return(list(ATE = plot1, mu1 = plot2, mu2 = plot3))
}

## Plot the effects forward in time: we fix the same carry-over length
## but vary the final period of interest (Y_T)
## Input a list of results of DCB, the periods and the time length

compute_ATE_plot_forward_time <- function(store_results, periods, Time){

  ATEs <- sapply(store_results, function(x) x$all_results$ATE )
  mu1 <- sapply(store_results, function(x) x$all_results$mu_hat[1] )
  mu2 <- sapply(store_results, function(x) x$all_results$mu_hat[2] )
  Var1 <- sapply(store_results, function(x) x$all_results$variances[1] )
  Var2 <- sapply(store_results, function(x) x$all_results$variances[2] )

  ATEs1 <- cbind(ATEs, Var1, Var2, periods)
  ATEs1 <- as.data.frame(ATEs1)
  names(ATEs1) <- c('ATE', 'Var1', 'Var2', 'Periods')
  ATEs1[,1] <- as.numeric(as.character(ATEs1[,1]))
  ATEs1[,2] <- as.numeric(as.character(ATEs1[,2]))
  ATEs1[,3] <- as.numeric(as.character(ATEs1[,3]))
  ATEs1[,4] <- factor(ATEs1[,4], levels = sort(unique(ATEs1[,4]), decreasing =  F))



  plot1 <- ggplot2::ggplot(ATEs1, ggplot2::aes(Periods, ATE)) +
    ggplot2::geom_segment(ggplot2::aes(y = ATE - sqrt(qchisq(0.9,  2 * Time)) * sqrt(Var1 + Var2),
                     x = Periods, yend = ATE + sqrt(qchisq(0.9, 2 * Time)) * sqrt(Var1 + Var2),
                     xend = Periods), size = 2, color = 'lightgray') +
    ggplot2::geom_segment(ggplot2::aes(y = ATE - 1.64 * sqrt(Var1 + Var2),
                     x = Periods, yend = ATE + 1.64 * sqrt(Var1 + Var2),
                     xend = Periods), size = 2, color = 'darkgray') +
    ggplot2::geom_point( size = 4, shape = 'triangle') +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
          plot.title = ggplot2::element_text(size=22),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 20))

  mus1 <- cbind(mu1, Var1, periods)
  mus1 <- as.data.frame(mus1)
  names(mus1) <- c('Mu1', 'Var1', 'Periods')
  mus1[,1] <- as.numeric(as.character(mus1[,1]))
  mus1[,2] <- as.numeric(as.character(mus1[,2]))
  mus1[,3] <- factor(mus1[,3], levels = sort(unique(mus1[,3]), decreasing =  F))


  plot2 <- ggplot2::ggplot(mus1, ggplot2::aes(Periods, Mu1)) +
    ggplot2::geom_segment(ggplot2::aes(y = Mu1 - sqrt(qchisq(0.9,  Time)) * sqrt(Var1),
                     x = Periods, yend = Mu1 + sqrt(qchisq(0.9, Time)) * sqrt(Var1),
                     xend = Periods), size = 2, color = 'lightgray') +
    ggplot2::geom_segment(ggplot2::aes(y = Mu1 - 1.64 * sqrt(Var1),
                     x = Periods, yend = Mu1 + 1.64 * sqrt(Var1),
                     xend = Periods), size = 2, color = 'darkgray') +
    ggplot2::geom_point( size = 4, shape = 'triangle') +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
          plot.title = ggplot2::element_text(size=22),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 20))

  mus2 <- cbind(mu2, Var2, periods)
  mus2 <- as.data.frame(mus2)
  names(mus2) <- c('Mu2', 'Var2', 'Periods')
  mus2[,1] <- as.numeric(as.character(mus2[,1]))
  mus2[,2] <- as.numeric(as.character(mus2[,2]))
  mus2[,3] <- factor(mus2[,3], levels = sort(unique(mus2[,3]), decreasing =  F))


  plot3 <- ggplot2::ggplot(mus2, ggplot2::aes(Periods, Mu2)) +
    ggplot2::geom_segment(ggplot2::aes(y = Mu2 - sqrt(qchisq(0.9,  Time)) * sqrt(Var2),
                     x = Periods, yend = Mu2 + sqrt(qchisq(0.9, Time)) * sqrt(Var2),
                     xend = Periods), size = 2, color = 'lightgray') +
    ggplot2::geom_segment(ggplot2::aes(y = Mu2 - 1.64 * sqrt(Var2),
                     x = Periods, yend = Mu2 + 1.64 * sqrt(Var2),
                     xend = Periods), size = 2, color = 'darkgray') +
    ggplot2::geom_point( size = 4, shape = 'triangle') +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title=ggplot2::element_text(size = 25), legend.text=ggplot2::element_text(size = 25),
          plot.title = ggplot2::element_text(size=22),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 20))

  return(list(ATE = plot1, mu1 = plot2, mu2 = plot3))
}

create_plot_treatments_helper <- function(D, individuals_unique, time_periods, axis.lab.gap = c(0,0)){

  df <- data.frame(treatment = c(D), unit = rep(individuals_unique, dim(D)[2]), Y= 1:length(D),
                   Time = rep(time_periods, each = length(individuals_unique) ))
  my_plot = panelView::panelView(Y ~ treatment, data = df,  index = c('unit','Time'), main = '', axis.lab.gap = axis.lab.gap)
  return(my_plot)
}

#' Plot treatment status of individuals over time
#'
#' @param panel_data data.frame format. It must be a panel with each row corresponding to a unit at a given date.
#' @param covariates_names vector of strings with column names of the covariates to include in the regression.
#' @param Time_name string denoting the column name of the panel corresponding to the date.
#' @param unit_name string denoting the column name of the panel corresponding to the individual.
#' @param outcome_name string denoting the column name of the panel corresponding to the outcome variable. Rows with missing outcomes will be removed.
#' @param treatment_name name of the column of the panel corresponding to the treatment assignment. Rows with missing treatment assignments will be removed.
#' @param final_period last period considered
#' @param initial_period first period considered
#' @param axis.lab.gap gap for x axis
#' @return plots

plot_treatment_status <- function(panel_data,
                                  Time_name, unit_name, outcome_name,
                                  treatment_name,
                                  final_period, initial_period, axis.lab.gap = c(0,0)){

  if(is.null(names(panel_data))){
    stop('Column names are missing in panel data.')
  }
  Time_name <- which(names(panel_data) == Time_name)
  unit_name <- which(names(panel_data) == unit_name)
  outcome_name <- which(names(panel_data) == outcome_name)
  treatment_name <- which(names(panel_data) == treatment_name)

  if(length(Time_name) == 0){
    stop('Wrong Time_name passed')
  }

  if(length(unit_name) == 0){
    stop('Wrong unit_name passed')
  }

  if(length(treatment_name) == 0){
    stop('Wrong treatment_name passed')
  }

  if(final_period %in% panel_data[, Time_name] == F){
    stop('Wrong final_period passed')
  }


  if(initial_period %in% panel_data[, Time_name] == F){
    stop('Wrong initial_period passed')
  }

  if(initial_period >= final_period){
    stop('initial_period is >= final_period. initial_period must be strictly smaller than final_period')
  }



  ## Remove units after the final period
  panel_data <- panel_data[panel_data[, Time_name] <= final_period,]
  ## Remove units before the initial period
  panel_data <- panel_data[panel_data[, Time_name] >= initial_period,]

  periods <- panel_data[, Time_name]
  individuals <- panel_data[, unit_name]
  individuals_unique <- unique(individuals)


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
    periods <- panel_data[, Time_name]
  }
  plot_treatments <- create_plot_treatments_helper(D, individuals_unique, time_periods = sort(unique(periods), decreasing = F),
                                                   axis.lab.gap)
  return(plot_treatments)
}
