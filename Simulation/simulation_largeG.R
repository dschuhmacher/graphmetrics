### Simulation for large graphs ###
library(graphmetrics)
# Setting of parameters

#n = rep(rep(c(20,30,50,100), each = 4))
n = c(20,30,50,100)
#m = rep(rep(c(20,30,50,100), times = 4))
m = c(20,30,50,100)
C = c(0.1,0.4) 


iterations <- 100

###### simulation ######
if(TRUE){
  result_large <- list()
  # for every scenario
  for(scenarioN in c(2,1)){ 
    # for every metric type 
    scenario_list_large <- list()
    for(type in c(  "OSPA2")){#"OSPA1",, "TT")){
      # setting data.frame with all parameter combinations
      if(scenarioN == 1){
        number_params <-length(n) * length(C)
        results_df <- data.frame(setting = 1:number_params, n = rep( n , each= length(C)), m = rep( m ,each=length(C)), C = rep(C, times = length(n)), res_FW = rep(0,number_params))
        times_df   <- data.frame(setting = 1:number_params, n = rep( n , each= length(C)), m = rep( m ,each=length(C)), C = rep(C, times = length(n)), time_auction = rep(0,number_params),  time_FW = rep(0,number_params))
      }
      else{
        number_params <- length(unique(n)) * length(C)
        results_df <- data.frame(setting = 1:number_params, n = rep( unique(n) , each= length(C)),  C = rep(C, times = length(unique(n))),res_FW = rep(0,number_params))
        times_df   <- data.frame(setting = 1:number_params, n = rep( unique(n) , each= length(C)),  C = rep(C, times = length(unique(n))), time_auction = rep(0,number_params),  time_FW = rep(0,number_params))
      }
      
      type_list_large <- list(type = type, mean = list(result_df = results_df, times_df = times_df), "5%_quant" = list(result_df = results_df, times_df = times_df),
                              "95%_quant" = list(result_df = results_df, times_df = times_df) , data = list())
      
      # for every parameter setting (filling df)
      for (i in 1:number_params) {
        cat("parameter setting number", i)
        type_list_large$data <- append(type_list_large$data, list(list(parameter_setting = i, results_iteration = data.frame(number = 1:iterations, res_auction = rep(0,iterations),  res_FW = rep(0,iterations)),
                                                                       time_iteration = data.frame(number = 1:iterations, time_auction = rep(0,iterations),  time_FW = rep(0,iterations)))))
        # iterations-many iterations for each parameter setting
        for (j in 1:iterations) {
          cat("Scenario",scenarioN, "parameter setting number", i, "iteration", j)
          if(scenarioN == 1){
            #g1 <- igraph::sample_grg(nodes = type_list_large$mean$result_df$n[i], radius = 3/type_list_large$mean$result_df$n[i], coords = TRUE) 
            g1 <- rspatER(n=type_list_large$mean$result_df$n[i] , p= 3/type_list_large$mean$result_df$n[i])
            #g2 <- igraph::sample_grg(nodes = type_list_large$mean$result_df$n[i], radius = 4/type_list_large$mean$result_df$n[i], coords = TRUE) 
            g2 <- rspatER(n=type_list_large$mean$result_df$m[i], p= 4/type_list_large$mean$result_df$m[i])
            c <- 4
          }
          else if(scenarioN ==2){
            g1 <- rspatER(n=type_list_large$mean$result_df$n[i], p= 4/type_list_large$mean$result_df$n[i])
            g2 <- rperturb(g1, scatter=1/type_list_large$mean$result_df$n[i], flip=0.3)
            c <- 2
          }
          type_list_large$data[[i]]$time_iteration$time_auction[j] <- system.time(type_list_large$data[[i]]$results_iteration$res_auction[j] <- gmspat(g1, g2, method = graph_auction, compensate=TRUE, type=type, CV = type_list_large$mean$result_df$C[i], CE = type_list_large$mean$result_df$C[i], vpen = type_list_large$mean$result_df$C[i], lang = "Cpp", stop_at=15, maxiter=10000, eps =0.01 * c * type_list_large$mean$result_df$C[i] * (type_list_large$mean$result_df$n[i]/100) , verbose=0)$dist)[3] # * #1/(type_list_large$mean$result_df$n[i]+1 ), verbose=0)$dist)[3]#0.01* 2 * type_list_large$mean$result_df$C[i] , verbose=0)$dist)[3]#* 2 
          cat("\n")
          type_list_large$data[[i]]$time_iteration$time_FW[j] <- system.time(type_list_large$data[[i]]$results_iteration$res_FW[j] <- gmspat(g1, g2, method = FW_match, start_help = FALSE, type=type, CV = type_list_large$mean$result_df$C[i], CE = type_list_large$mean$result_df$C[i], vpen = type_list_large$mean$result_df$C[i])$dist)[3]
          #type_list_large$data[[i]]$time_iteration$time_FW_start[j] <- system.time(type_list_large$data[[i]]$results_iteration$res_FW_start[j] <- gmspat(g1, g2, method = FW_match, start_help = TRUE, type=type, CV = type_list_large$mean$result_df$C[i], CE = type_list_large$mean$result_df$C[i], vpen = type_list_large$mean$result_df$C[i])$dist)[3]
        } 
        
        ##### compute mean time and mean deviation for current parameter setting #####
        
        #### mean relative deviation from auction result, i.e. (res_approx - res_auction )/res_auction #### 
        # with 5% and 95% quantile 
        
        #auction 
        type_list_large$mean$result_df$res_auction[i] <- mean(type_list_large$data[[i]]$results_iteration$res_auction )
        type_list_large$'5%_quant'$result_df$res_auction[i] <- sort((type_list_large$data[[i]]$results_iteration$res_auction ))[iterations*0.05]
        type_list_large$'95%_quant'$result_df$res_auction[i] <- sort((type_list_large$data[[i]]$results_iteration$res_auction ))[iterations*0.95]
        
        #auction - FW
        type_list_large$mean$result_df$res_FW[i] <- mean((type_list_large$data[[i]]$results_iteration$res_FW - type_list_large$data[[i]]$results_iteration$res_auction)/type_list_large$data[[i]]$results_iteration$res_auction)
        type_list_large$'5%_quant'$result_df$res_FW[i] <- sort((type_list_large$data[[i]]$results_iteration$res_FW - type_list_large$data[[i]]$results_iteration$res_auction)/type_list_large$data[[i]]$results_iteration$res_auction)[iterations*0.05]
        type_list_large$'95%_quant'$result_df$res_FW[i] <- sort((type_list_large$data[[i]]$results_iteration$res_FW - type_list_large$data[[i]]$results_iteration$res_auction)/type_list_large$data[[i]]$results_iteration$res_auction)[iterations*0.95]
        
        #auction - FW_start
        #type_list_large$mean$result_df$res_FW_start[i] <- mean((type_list_large$data[[i]]$results_iteration$res_FW_start - type_list_large$data[[i]]$results_iteration$res_auction)/type_list_large$data[[i]]$results_iteration$res_auction)
        #type_list_large$'5%_quant'$result_df$res_FW_start[i] <- sort((type_list_large$data[[i]]$results_iteration$res_FW_start - type_list_large$data[[i]]$results_iteration$res_auction)/type_list_large$data[[i]]$results_iteration$res_auction)[iterations*0.05]
        #type_list_large$'95%_quant'$result_df$res_FW_start[i] <- sort((type_list_large$data[[i]]$results_iteration$res_FW_start - type_list_large$data[[i]]$results_iteration$res_auction)/type_list_large$data[[i]]$results_iteration$res_auction)[iterations*0.95]
        
        #### mean time #### 
        # with 5% and 95% quantile 
        #auction
        type_list_large$mean$times_df$time_auction[i] <- mean(type_list_large$data[[i]]$time_iteration$time_auction)
        type_list_large$'5%_quant'$times_df$time_auction[i] <- sort(type_list_large$data[[i]]$time_iteration$time_auction)[iterations*0.05]
        type_list_large$'95%_quant'$times_df$time_auction[i] <- sort(type_list_large$data[[i]]$time_iteration$time_auction)[iterations*0.95]
        
        #FW
        type_list_large$mean$times_df$time_FW[i] <- mean(type_list_large$data[[i]]$time_iteration$time_FW)
        type_list_large$'5%_quant'$times_df$time_FW[i] <- sort(type_list_large$data[[i]]$time_iteration$time_FW)[iterations*0.05]
        type_list_large$'95%_quant'$times_df$time_FW[i] <- sort(type_list_large$data[[i]]$time_iteration$time_FW)[iterations*0.95]
        
        #FW_start
        #type_list_large$mean$times_df$time_FW_start[i] <- mean(type_list_large$data[[i]]$time_iteration$time_FW_start)
        #type_list_large$'5%_quant'$times_df$time_FW_start[i] <- sort(type_list_large$data[[i]]$time_iteration$time_FW_start)[iterations*0.05]
        #type_list_large$'95%_quant'$times_df$time_FW_start[i] <- sort(type_list_large$data[[i]]$time_iteration$time_FW_start)[iterations*0.95]
        
        
        #### save ####
        #save(type_list_large, file = "type_list_data_large.RData")
      }
      scenario_list_large <- append(scenario_list_large, list(type_list_large))
      save(scenario_list_large, file = "scenario_data_large.RData")
    }
    result_large <- append(result_large, list(scenario_list_large))
    save(result_large, file = "simulation_data_large.RData")
  }
  
}