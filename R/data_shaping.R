## Clean data : remove units whose periods of interest

help_clean_data_missing_time = function(i, times, time_1, time_T){
  if(length(times) == length(time_1:time_T)) return(i)
  return(NA)
}

clean_data_from_missing_time <- function(panel_data, unit_name, time_name){
  units = unique(panel_data[, unit_name])
  time_1 = min(panel_data[, time_name])
  time_T = max(panel_data[, time_name])
  units = sapply(units, function(x) help_clean_data_missing_time(x, panel_data[panel_data[,unit_name] == x, time_name], time_1, time_T))
  return(panel_data[panel_data[, unit_name] %in% na.omit(units),])
}



## Inputs: column name for units
##         column name for treatment assignment
##         column name for period
##         panel_data
## Output: matrix of treatment assignments for each unit

create_matrix_of_D <- function(unit_name,
                               treatment_name,
                               Time_name,
                               panel_data){

  units <- panel_data[, unit_name]
  periods <- panel_data[, Time_name]
  D <- matrix(NA, nrow = length(unique(units)), ncol = length(unique(periods)))
  k = 1
  for(t in sort(unique(periods), decreasing = F) ){
    D_this_year <- rep(NA, length(unique(units)))
    D_this_year[unique(units) %in% panel_data[panel_data[, Time_name] == t, unit_name]  ] <- panel_data[panel_data[, Time_name] == t, treatment_name]
    D[,k] <-   D_this_year
    k = k + 1
  }
  return(D)
}

## inputs: column name corresponding to the unit
##         covariates names used
##         column name of the period
##         panel_data
## Output: list of covariates, each element is a matrix corresponding to covariates in a certain period

create_matrix_of_covariates <- function(unit_name, covariates_names, Time_name, panel_data, dim_fe){

  units <- panel_data[, unit_name]
  unique_units <- unique(units)
  n <- length(unique_units)
  periods <- panel_data[, Time_name]
  covariates <- list()
  final_period = max(unique(periods))
  k = 1
  if(dim_fe == 0){
  for(t in sort(unique(periods), decreasing = F) ){
    covariates_today <- matrix(NA, nrow = n, ncol = length(covariates_names))
    covariates_today[unique_units %in% panel_data[periods == t, unit_name], ] <- t(apply(as.matrix(panel_data[periods == t, names(panel_data) %in% covariates_names ]), 1, as.numeric))
    covariates[[k]] <- covariates_today
    k = k + 1
  }
  } else {
    ## Include fe only in the final regression if dim_fe > 0 (i.e., params$demeaned_fe = T)
    ## You then remove the fixed effects once predicting the previous regressions
    for(t in sort(unique(periods), decreasing = F) ){
      keep_covariates = covariates_names[-c((length(covariates_names) - dim_fe + 1):length(covariates_names))]
      if(t == final_period) keep_covariates = covariates_names
      covariates_today <- matrix(NA, nrow = n, ncol = length(keep_covariates))
      covariates_today[unique_units %in% panel_data[periods == t, unit_name], ] <- t(apply(as.matrix(panel_data[periods == t, names(panel_data) %in% keep_covariates ]), 1, as.numeric))
      covariates[[k]] <- covariates_today
      k = k + 1
    }
  }
  return(covariates)
}


## create a sample with pooled data
## the new data has a matrix year, corresponding to the period of interest
## and new_Time, corresponding to the original entry of the column
## the name of each element is the name of the country lagged

create_pooled_matrix = function(panel_data, Time_name, unit_name, outcome_name,
                                treatment_name,
                                final_period, initial_period, pooled_all = F,
                                length_treatment){

  original_name = names(panel_data)
  num_periods = final_period - initial_period
  if(pooled_all) num_periods = final_period - min(panel_data[,Time_name]) - length_treatment - 1
  individuals = unique(panel_data[, unit_name])
  k = 1
  for(i in individuals){

    columns = panel_data[panel_data[,unit_name] == i, ]
    flag = T
    k = 1
    for(j in 1:num_periods){
      keep_periods_final = final_period - j
      keep_periods_init = final_period - j - length_treatment + 1
      new_element = columns[columns[, Time_name] %in% c(keep_periods_init:keep_periods_final) ,]
      if(dim(new_element)[1]  == length(c(keep_periods_init:keep_periods_final))){
        flag = F
        new_element =   cbind(new_element, new_Time = new_element[, Time_name],
                              new_name = paste0(i, 'l',j, 'unit', k))

        new_element[sort(new_element[,Time_name], decreasing = F, index.return = T)$ix, Time_name] = c(keep_periods_init:keep_periods_final) + j
        if(j == 1) init = new_element
        if(j > 1) init <- rbind(new_element, init)
        k = k +1
      }
    }

    if(i == individuals[1]){
      acc = init
    } else if (flag == F) {
      acc = rbind(acc, init)
    }
  }
  panel_data_new = cbind(panel_data, new_Time = panel_data[, Time_name], new_name = panel_data[, unit_name])
  final_data = rbind(panel_data_new, acc)
  final_data = as.data.frame(final_data)
  names(final_data) = c(original_name, 'new_Time', 'new_name')
  return(final_data)
}


