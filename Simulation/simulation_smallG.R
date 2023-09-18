## Simulation to compare Auction and FW with exact solution for small graphs ##

library(graphmetrics)
# Setting of parameters
n = rep(rep(c(4,8,11), each = 3))

m = rep(rep(c(4,8,11), times = 3))
C = c(0.1,0.4,0.8) 

iterations <- 100
#sessionInfo()
result <- list()

###### simulation ######
if(TRUE){
  # for every scenario (independent (1) or dependent (2))
  for(scenarioN in c(1,2)){ 
    # for every metric type  
    scenario_list <- list()
    for(type in c(  "OSPA2")){#"OSPA1", "TT")){
      # setting data.frame with all parameter combinations
      if(scenarioN == 1){
        number_params <-length(n) * length(C)
        results_df <- data.frame(setting = 1:number_params, n = rep( n , each= length(C)), m = rep( m ,each=length(C)), C = rep(C, times = length(n)), res_auction = rep(0,number_params),  res_FW = rep(0,number_params), res_FW_start = rep(0,number_params))
        times_df   <- data.frame(setting = 1:number_params, n = rep( n , each= length(C)), m = rep( m ,each=length(C)), C = rep(C, times = length(n)), time_exact = rep(0,number_params), time_auction = rep(0,number_params),  time_FW = rep(0,number_params), time_FW_start = rep(0,number_params))
      }
      else{
        number_params <- length(unique(n)) * length(C)
        results_df <- data.frame(setting = 1:number_params, n = rep( unique(n) , each= length(C)),  C = rep(C, times = length(unique(n))), res_auction = rep(0,number_params), res_auction = rep(0,number_params),  res_FW = rep(0,number_params), res_FW_start = rep(0,number_params))
        times_df   <- data.frame(setting = 1:number_params, n = rep( unique(n) , each= length(C)),  C = rep(C, times = length(unique(n))), time_exact = rep(0,number_params), time_auction = rep(0,number_params), time_auction = rep(0,number_params),  time_FW = rep(0,number_params), time_FW_start = rep(0,number_params))
      }
      
      type_list <- list(type = type, mean = list(result_df = results_df, times_df = times_df), "5%_quant" = list(result_df = results_df, times_df = times_df),
                        "95%_quant" = list(result_df = results_df, times_df = times_df) , data = list())
      
      # for every parameter setting (filling df)
      for (i in 1:number_params) {
        cat("parameter setting number", i)
        type_list$data <- append(type_list$data, list(list(parameter_setting = i, results_iteration = data.frame(number = 1:iterations, res_exact = rep(0,iterations), res_auction = rep(0,iterations),  res_auction = rep(0,iterations),  res_FW = rep(0,iterations), res_FW_start = rep(0,iterations)),
                                                           time_iteration = data.frame(number = 1:iterations, time_exact = rep(0,iterations), time_auction = rep(0,iterations),  time_auction = rep(0,iterations),  time_FW = rep(0,iterations), time_FW_start = rep(0,iterations)))))
        # iterations-many iterations for each parameter setting
        for (j in 1:iterations) {
          cat("Scenario",scenarioN, "parameter setting number", i, "iteration", j)
          if(scenarioN == 1){
            g1 <- rspatER(n=type_list$mean$result_df$n[i] , p=0.3)
            g2 <- rspatER(n=type_list$mean$result_df$m[i], p=0.4)
            
          }
          else if(scenarioN ==2){
            g1 <- rspatER(n=type_list$mean$result_df$n[i], p=0.4)
            g2 <- rperturb(g1, scatter=0.1, flip=0.3)
          }
          type_list$data[[i]]$time_iteration$time_exact[j] <- system.time(type_list$data[[i]]$results_iteration$res_exact[j] <- gmspat(g1, g2, method = cplex_match, type=type, CV = type_list$mean$result_df$C[i], CE = type_list$mean$result_df$C[i], vpen = type_list$mean$result_df$C[i])$dist)[3]
          type_list$data[[i]]$time_iteration$time_auction[j] <- system.time(type_list$data[[i]]$results_iteration$res_auction[j] <- gmspat(g1, g2, method = graph_auction, compensate=TRUE, type=type, CV = type_list$mean$result_df$C[i], CE = type_list$mean$result_df$C[i], vpen = type_list$mean$result_df$C[i], lang = "Cpp", stop_at=3, maxiter=100, eps = 0.01 * type_list$mean$result_df$C[i], verbose=0)$dist)[3]
          type_list$data[[i]]$time_iteration$time_FW[j] <- system.time(type_list$data[[i]]$results_iteration$res_FW[j] <- gmspat(g1, g2, method = FW_match, start_help = FALSE, type=type, CV = type_list$mean$result_df$C[i], CE = type_list$mean$result_df$C[i], vpen = type_list$mean$result_df$C[i])$dist)[3]
        } 
        
        ##### compute mean time and mean deviation for current parameter setting #####
        
        #### mean relative deviation, i.e. (res_approx - res_exact )/res_exact #### 
        # with 5% and 95% quantile 
        
        # vector of iterations-many exact results for current parameter setting
        res_exact <- type_list$data[[i]]$results_iteration$res_exact 
        
        #auction
        type_list$mean$result_df$res_auction[i] <- mean((type_list$data[[i]]$results_iteration$res_auction - res_exact)/res_exact)
        type_list$'5%_quant'$result_df$res_auction[i] <- sort((type_list$data[[i]]$results_iteration$res_auction - res_exact)/res_exact)[iterations*0.05]
        type_list$'95%_quant'$result_df$res_auction[i] <- sort((type_list$data[[i]]$results_iteration$res_auction - res_exact)/res_exact)[iterations*0.95]
        
        #FW
        type_list$mean$result_df$res_FW[i] <- mean((type_list$data[[i]]$results_iteration$res_FW - res_exact)/res_exact)
        type_list$'5%_quant'$result_df$res_FW[i] <- sort((type_list$data[[i]]$results_iteration$res_FW - res_exact)/res_exact)[iterations*0.05]
        type_list$'95%_quant'$result_df$res_FW[i] <- sort((type_list$data[[i]]$results_iteration$res_FW - res_exact)/res_exact)[iterations*0.95]
        
        #### mean time #### 
        # with 5% and 95% quantile 
        
        #exact
        type_list$mean$times_df$time_exact[i] <- mean(type_list$data[[i]]$time_iteration$time_exact )
        type_list$'5%_quant'$times_df$time_exact[i] <- sort(type_list$data[[i]]$time_iteration$time_exact)[iterations*0.05]
        type_list$'95%_quant'$times_df$time_exact[i] <- sort(type_list$data[[i]]$time_iteration$time_exact)[iterations*0.95]
        
        #auction
        type_list$mean$times_df$time_auction[i] <- mean(type_list$data[[i]]$time_iteration$time_auction)
        type_list$'5%_quant'$times_df$time_auction[i] <- sort(type_list$data[[i]]$time_iteration$time_auction)[iterations*0.05]
        type_list$'95%_quant'$times_df$time_auction[i] <- sort(type_list$data[[i]]$time_iteration$time_auction)[iterations*0.95]
        
        #FW
        type_list$mean$times_df$time_FW[i] <- mean(type_list$data[[i]]$time_iteration$time_FW)
        type_list$'5%_quant'$times_df$time_FW[i] <- sort(type_list$data[[i]]$time_iteration$time_FW)[iterations*0.05]
        type_list$'95%_quant'$times_df$time_FW[i] <- sort(type_list$data[[i]]$time_iteration$time_FW)[iterations*0.95]
        
        
        #### save ####
        type_list <- type_list
        #save(type_list, file = "type_list_data_small.RData")
      }
      scenario_list <- append(scenario_list, list(type_list))
      #save(scenario_list, file = "scenario_data_small.RData")
    }
    result <- append(result, list(scenario_list))
    save(result, file = "simulation_data_small.RData")
  }
  
}
