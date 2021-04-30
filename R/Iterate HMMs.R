
run.HMMs = function(list1, elements) {
  hmm.res<- list()
  
  for (j in elements) {
    start.time<- Sys.time()
    
    # Empty list for order selection
    k.models<- list()
    
    
    ## K = 2
    
    allm<- list()
    niter<- 30
    stateNames <- c("Encamped","Exploratory")
    whichzero <- which(list1[[j]]$step == 0)
    propzero <- length(whichzero)/nrow(list1[[j]])
    zeromass0 <- c(propzero, 0)        #for zero distances by state
    
    for (i in 1:niter) {
      print(paste("Simulation", j))
      print(paste("K=2"))
      print(paste("Iteration", i))
      
      # Step length mean
      stepMean0 <- runif(2,
                         min = c(0.1, 2),
                         max = c(1, 15))
      # Step length standard deviation
      stepSD0 <- runif(2,
                       min = c(0.1, 3),
                       max = c(2, 8))
      # Turning angle mean
      angleMean0 <- c(pi, 0)
      # Turning angle concentration
      angleCon0 <- c(0.8, 0.8)
      # Fit model
      if(propzero > 0) {  #don't include zero mass if no 0s present
        stepPar0 <- c(stepMean0, stepSD0, zeromass0)
      } else {
        stepPar0 <- c(stepMean0, stepSD0)
      }
      anglePar0 <- c(angleMean0, angleCon0)
      allm[[i]] <- fitHMM(data = list1[[j]], nbStates = 2, 
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          dist = list(step = "gamma", angle = "wrpcauchy"),
                          formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                          estAngleMean = list(angle=TRUE),
                          stateNames = stateNames,
                          optMethod = "Nelder-Mead")
    }
    
    # Extract likelihoods of fitted models
    allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))
    
    # Index of best fitting model (smallest negative log-likelihood)
    whichbest <- which.min(allnllk)
    
    # Best fitting model
    k.models[[1]] <- allm[[whichbest]]
    
    
    
    
    ## K = 3
    
    allm<- list()
    niter<- 30
    stateNames <- c("Encamped","ARS","Transit")
    whichzero <- which(list1[[j]]$step == 0)
    propzero <- length(whichzero)/nrow(list1[[j]])
    zeromass0 <- c(propzero, 0, 0)        #for zero distances by state
    
    for (i in 1:niter) {
      print(paste("Simulation", j))
      print(paste("K=3"))
      print(paste("Iteration", i))
      
      # Step length mean
      stepMean0 <- runif(3,
                         min = c(0.1, 1, 5),
                         max = c(0.5, 3, 15))
      # Step length standard deviation
      stepSD0 <- runif(3,
                       min = c(0.1, 1, 2),
                       max = c(1, 2, 5))
      # Turning angle mean
      angleMean0 <- c(pi, 0, 0)
      # Turning angle concentration
      angleCon0 <- runif(3,
                         min = c(0.7, 0.01, 0.7),
                         max = c(0.99, 0.1, 0.99))
      # Fit model
      if(propzero > 0) {  #don't include zero mass if no 0s present
        stepPar0 <- c(stepMean0, stepSD0, zeromass0)
      } else {
        stepPar0 <- c(stepMean0, stepSD0)
      }
      anglePar0 <- c(angleMean0, angleCon0)
      allm[[i]] <- fitHMM(data = list1[[j]], nbStates = 3, 
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          dist = list(step = "gamma", angle = "wrpcauchy"),
                          formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                          estAngleMean = list(angle=TRUE),
                          stateNames = stateNames,
                          optMethod = "Nelder-Mead")  
    }
    
    # Extract likelihoods of fitted models
    allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))
    
    # Index of best fitting model (smallest negative log-likelihood)
    whichbest <- which.min(allnllk)
    
    # Best fitting model
    k.models[[2]] <- allm[[whichbest]]
    
    
    
    
    
    ## K = 4
    
    allm<- list()
    niter<- 30
    stateNames <- c("Encamped","ARS","Exploratory","Transit")
    whichzero <- which(list1[[j]]$step == 0)
    propzero <- length(whichzero)/nrow(list1[[j]])
    zeromass0 <- c(propzero, propzero, 0, 0)        #for zero distances by state
    
    for (i in 1:niter) {
      print(paste("Simulation", j))
      print(paste("K=4"))
      print(paste("Iteration", i))
      
      # Step length mean
      stepMean0 <- runif(4,
                         min = c(0.1, 1, 5, 8),
                         max = c(1, 3, 10, 15))
      # Step length standard deviation
      stepSD0 <- runif(4,
                       min = c(0.1, 1, 2, 3),
                       max = c(1, 2, 4, 8))
      # Turning angle mean
      angleMean0 <- c(pi, pi, 0, 0)
      # Turning angle concentration
      angleCon0 <- runif(4,
                         min = c(0.5, 0.01, 0.5, 0.5),
                         max = c(0.99, 0.5, 0.9, 0.99))
      # Fit model
      if(propzero > 0) {  #don't include zero mass if no 0s present
        stepPar0 <- c(stepMean0, stepSD0, zeromass0)
      } else {
        stepPar0 <- c(stepMean0, stepSD0)
      }
      anglePar0 <- c(angleMean0, angleCon0)
      allm[[i]] <- fitHMM(data = list1[[j]], nbStates = 4, 
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          dist = list(step = "gamma", angle = "wrpcauchy"),
                          formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                          estAngleMean = list(angle=TRUE),
                          stateNames = stateNames,
                          optMethod = "Nelder-Mead")  
    }
    
    # Extract likelihoods of fitted models
    allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))
    
    # Index of best fitting model (smallest negative log-likelihood)
    whichbest <- which.min(allnllk)
    
    # Best fitting model
    k.models[[3]] <- allm[[whichbest]]
    
    
    
    end.time<- Sys.time()
    elapsed.time<- difftime(end.time, start.time, units = "min")
    
    hmm.res[[j]]<- list(models = k.models, elapsed.time = elapsed.time)
  }
  
  hmm.res
}
