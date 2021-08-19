### CRW model

CRW.sim=function(nsim, ntseg, nstep, SL.params, TA.params, Z0) {  
  #nsim=number of simulations
  #ntseg=number of time segments
  #nstep=number of steps per time segment; if single number, applied to all segments; if vector of length 2, these values used as min and max to be drawn from runif() per segment
  #nbehav=number of behavioral states ==> derive from nrow SL.params
  #SL.params=df of shape and scale params
  #TA.params=df of mean TA and concen. param
  #Z0=initial location
  
  #uses gamma and wrapped cauchy distribs
  #behaviors params must be in order
  #for simulating w/ 3 behavioral states
  
  track.list<- list()  #store all tracks of all durations
  brkpt.list<- list()  #store all true breakpoints
  
  for (k in 1:length(ntseg)) {
    
    behav<- list()  #to store behaviors for each time segment
    behav.all<- list()  #to store behaviors for each observation
    sim.list<- list()  #to store all sim tracks of same ntseg
    sim.brkpts<- list()  #to store all sim breakpoints of same ntseg
    
    #randomly choose n segments of each behavior to be 'pure' instead of mixed
    npure<- ifelse(ntseg[k] <= 10, 0,
                   ifelse(ntseg[k] > 10 & ntseg[k] <=50, 3,
                          ifelse(ntseg[k] > 50 & ntseg[k] <=100, 6,
                                 9)))
    
    for (j in 1:nsim) {
      
      #create vector of dominant behaviors per time segment (50 segments)
      behav[[j]]<- rmultinom(ntseg[k], 1, c(0.5, 0.30, 0.20)) %>% 
        apply(., 2, function(x) which(x == 1)) #change from matrix to vector
      # table(behav[[k]])/ntseg[k] #check freq
      
      #assign positions for pure time segments (if any)
      behav1.pure<- sample(which(behav[[j]]==1), npure, replace = FALSE)
      behav2.pure<- sample(which(behav[[j]]==2), npure, replace = FALSE)
      behav3.pure<- sample(which(behav[[j]]==3), npure, replace = FALSE)
      
      pure<- c(behav1.pure, behav2.pure, behav3.pure) %>% sort()
      
      #define number of steps per segment
      if (length(nstep) == 1) {
        nstep1<- nstep
      } else if (length(nstep) == 2) {
        nstep1<- round(runif(ntseg[k], nstep[1], nstep[2]))
      } else {
        stop("nstep must be length of 1 or 2")
      }
      
      #create vector of behaviors within each time segment
      behav.full<- list()
      for (i in 1:length(behav[[j]])) {
        nstep.i<- ifelse(length(nstep1) == 1, nstep1, nstep1[i])  #define nstep per tseg
        
        if (i %in% pure) {
          behav.full[[i]]<- rep(behav[[j]][i], nstep.i)
        }else if (behav[[j]][i] == 1) {
          behav.full[[i]]<- rmultinom(nstep.i, 1, c(0.8, 0.1, 0.1))  %>% 
            apply(., 2, function(x) which(x == 1))
        } else if (behav[[j]][i] == 2) {
          behav.full[[i]]<- rmultinom(nstep.i, 1, c(0.1, 0.8, 0.1))  %>% 
            apply(., 2, function(x) which(x == 1))
        } else if (behav[[j]][i] == 3) {
          behav.full[[i]]<- rmultinom(nstep.i, 1, c(0.1, 0.1, 0.8))  %>% 
            apply(., 2, function(x) which(x == 1))
        }
      }
      
      behav.all[[j]]<- unlist(behav.full)
      
      
      
      #create vector of step lengths
      SL<- list()
      for (i in 1:length(behav.all[[j]])) {
        if (behav.all[[j]][i] == 1) {
          SL[[i]]<- rgamma(1, shape = SL.params[1,1], scale = SL.params[1,2])  #Rest
        } else if (behav.all[[j]][i] == 2) {
          SL[[i]]<- rgamma(1, shape = SL.params[2,1], scale = SL.params[2,2])  #Exploratory
        } else {
          SL[[i]]<- rgamma(1, shape = SL.params[3,1], scale = SL.params[3,2])  #Transit
        }
      }
      SL<- unlist(SL)
      
      
      #create vector of turning angles
      TA<- list()
      for (i in 1:length(behav.all[[j]])) {
        if (behav.all[[j]][i] == 1) {
          TA[[i]]<- rwrappedcauchy(1, mu=circular(TA.params[1,1]), rho=TA.params[1,2]) %>%
            ifelse(. > pi, .-(2*pi), .)  #Rest
        } else if (behav.all[[j]][i] == 2) {
          TA[[i]]<- rwrappedcauchy(1, mu=circular(TA.params[2,1]), rho=TA.params[2,2]) %>%
            ifelse(. > pi, .-(2*pi), .)  #Exploratory
        } else {
          TA[[i]]<- rwrappedcauchy(1, mu=circular(TA.params[3,1]), rho=TA.params[3,2]) %>%
            ifelse(. > pi, .-(2*pi), .)  #Transit
        }
      }
      TA<- unlist(TA)
      
      
      # cumulative angle
      Phi <- cumsum(TA)
      
      # step length components
      dX <- SL*cos(Phi)
      dY <- SL*sin(Phi)
      
      # actual X-Y values
      X <- c(Z0[1], Z0[1] + cumsum(dX))
      Y <- c(Z0[2], Z0[2] + cumsum(dY))
      
      #define coarse-scale behav
      if (length(nstep) == 2) {
        behav_coarse1<- rep(behav[[j]], nstep1)
      } else {
        behav_coarse1<- rep(behav[[j]], each=nstep1)
      }
      
      
      track<- data.frame(id = as.character(paste0(j,"_",k)),
                         x = X,
                         y = Y,
                         SL = c(NA, SL),
                         TA = c(NA, TA),
                         behav_fine = as.factor(c(NA, behav.all[[j]])),
                         behav_coarse = as.factor(c(NA, behav_coarse1)),
                         track_length = length(behav.all[[j]]))
      
      sim.list[[j]]<- track
      
      #define breakpoints based on how nstep was set
      if (length(nstep) == 1) {
        sim.brkpts[[j]]<- which(diff(behav[[j]]) != 0) * nstep1
      } else {
        sim.brkpts[[j]]<- cumsum(nstep1)[which(diff(behav[[j]]) != 0)]
      }
      
    }
    names(sim.list)<- sim.list %>% modify_depth(1, ~unique(.$id)) %>% unlist()
    names(sim.brkpts)<- sim.list %>% modify_depth(1, ~unique(.$id)) %>% unlist()
    
    track.list[[k]]<- sim.list
    brkpt.list[[k]]<- sim.brkpts
  }
  names(track.list)<- ntseg
  names(brkpt.list)<- ntseg
  
  list(tracks = track.list, brkpts = brkpt.list)
}

#----------------------------------------------


### CRW model using weird distributions

CRW.sim2=function(nsim, ntseg, nstep, SL.params, TA.params, Z0) {  
  #nsim=number of simulations
  #ntseg=number of time segments
  #nstep=number of steps per time segment; if single number, applied to all segments; if vector of length 2, these values used as min and max to be drawn from runif() per segment
  #nbehav=number of behavioral states ==> derive from nrow SL.params
  #SL.params=df of shape and scale params
  #TA.params=df of mean TA and concen. param
  #Z0=initial location
  
  #uses truncated normal for SL and truncated normal, uniform, and beta for TA
  #behaviors params must be in order
  #for simulating w/ 3 behavioral states
  
  track.list<- list()  #store all tracks of all durations
  brkpt.list<- list()  #store all true breakpoints
  
  for (k in 1:length(ntseg)) {
    
    behav<- list()  #to store behaviors for each time segment
    behav.all<- list()  #to store behaviors for each observation
    sim.list<- list()  #to store all sim tracks of same ntseg
    sim.brkpts<- list()  #to store all sim breakpoints of same ntseg
    
    #randomly choose n segments of each behavior to be 'pure' instead of mixed
    npure<- ifelse(ntseg[k] <= 10, 0,
                   ifelse(ntseg[k] > 10 & ntseg[k] <=50, 3,
                          ifelse(ntseg[k] > 50 & ntseg[k] <=100, 6,
                                 9)))
    
    for (j in 1:nsim) {
      
      #create vector of dominant behaviors per time segment (50 segments)
      behav[[j]]<- rmultinom(ntseg[k], 1, c(0.5, 0.30, 0.20)) %>% 
        apply(., 2, function(x) which(x == 1)) #change from matrix to vector
      # table(behav[[k]])/ntseg[k] #check freq
      
      #assign positions for pure time segments (if any)
      behav1.pure<- sample(which(behav[[j]]==1), npure, replace = FALSE)
      behav2.pure<- sample(which(behav[[j]]==2), npure, replace = FALSE)
      behav3.pure<- sample(which(behav[[j]]==3), npure, replace = FALSE)
      
      pure<- c(behav1.pure, behav2.pure, behav3.pure) %>% sort()
      
      #define number of steps per segment
      if (length(nstep) == 1) {
        nstep1<- nstep
      } else if (length(nstep) == 2) {
        nstep1<- round(runif(ntseg[k], nstep[1], nstep[2]))
      } else {
        stop("nstep must be length of 1 or 2")
      }
      
      
      #create vector of behaviors within each time segment
      behav.full<- list()
      for (i in 1:length(behav[[j]])) {
        nstep.i<- ifelse(length(nstep1) == 1, nstep1, nstep1[i])  #define nstep per tseg
        
        if (i %in% pure) {
          behav.full[[i]]<- rep(behav[[j]][i], nstep.i)
        }else if (behav[[j]][i] == 1) {
          behav.full[[i]]<- rmultinom(nstep.i, 1, c(0.8, 0.1, 0.1))  %>% 
            apply(., 2, function(x) which(x == 1))
        } else if (behav[[j]][i] == 2) {
          behav.full[[i]]<- rmultinom(nstep.i, 1, c(0.1, 0.8, 0.1))  %>% 
            apply(., 2, function(x) which(x == 1))
        } else if (behav[[j]][i] == 3) {
          behav.full[[i]]<- rmultinom(nstep.i, 1, c(0.1, 0.1, 0.8))  %>% 
            apply(., 2, function(x) which(x == 1))
        }
      }
      
      behav.all[[j]]<- unlist(behav.full)
      
      
      
      #create vector of step lengths
      SL<- list()
      for (i in 1:length(behav.all[[j]])) {
        if (behav.all[[j]][i] == 1) {
          SL[[i]]<- rtnorm(1, lo = SL.params[1,1], hi = SL.params[1,2], mu = SL.params[1,3],
                           sig = SL.params[1,4])  #Rest
        } else if (behav.all[[j]][i] == 2) {
          SL[[i]]<- rtnorm(1, lo = SL.params[2,1], hi = SL.params[2,2], mu = SL.params[2,3],
                           sig = SL.params[2,4])  #Exploratory
        } else {
          SL[[i]]<- rtnorm(1, lo = SL.params[3,1], hi = SL.params[3,2], mu = SL.params[3,3],
                           sig = SL.params[3,4])  #Transit
        }
      }
      SL<- unlist(SL)
      
      
      #create vector of turning angles
      TA<- list()
      for (i in 1:length(behav.all[[j]])) {
        if (behav.all[[j]][i] == 1) {
          TA[[i]]<- rbeta(1, shape1 = TA.params[1,1], shape2 = TA.params[1,2])*2*pi-pi  #Rest
        } else if (behav.all[[j]][i] == 2) {
          TA[[i]]<- runif(1, min = TA.params[2,1], max = TA.params[2,2])  #Exploratory
        } else {
          TA[[i]]<- rtnorm(1, lo=-pi, hi=pi, mu=TA.params[3,1], sig=TA.params[3,2])  #Transit
          
        }
      }
      TA<- unlist(TA)
      
      
      # cumulative angle
      Phi <- cumsum(TA)
      
      # step length components
      dX <- SL*cos(Phi)
      dY <- SL*sin(Phi)
      
      # actual X-Y values
      X <- c(Z0[1], Z0[1] + cumsum(dX))
      Y <- c(Z0[2], Z0[2] + cumsum(dY))
      
      #define coarse-scale behav
      if (length(nstep) == 2) {
        behav_coarse1<- rep(behav[[j]], nstep1)
      } else {
        behav_coarse1<- rep(behav[[j]], each=nstep1)
      }
      
      track<- data.frame(id = as.character(paste0(j,"_",k)),
                         x = X,
                         y = Y,
                         SL = c(NA,SL),
                         TA = c(NA, TA),
                         behav_fine = as.factor(c(NA, behav.all[[j]])),
                         behav_coarse = as.factor(c(NA, behav_coarse1)),
                         track_length = length(behav.all[[j]]))
      
      sim.list[[j]]<- track
      
      #define breakpoints based on how nstep was set
      if (length(nstep) == 1) {
        sim.brkpts[[j]]<- which(diff(behav[[j]]) != 0) * nstep1
      } else {
        sim.brkpts[[j]]<- cumsum(nstep1)[which(diff(behav[[j]]) != 0)]
      }
      
    }
    names(sim.list)<- sim.list %>% modify_depth(1, ~unique(.$id)) %>% unlist()
    names(sim.brkpts)<- sim.list %>% modify_depth(1, ~unique(.$id)) %>% unlist()
    
    track.list[[k]]<- sim.list
    brkpt.list[[k]]<- sim.brkpts
  }
  names(track.list)<- ntseg
  names(brkpt.list)<- ntseg
  
  list(tracks = track.list, brkpts = brkpt.list)
}

#----------------------------------------------


### HMM assuming CRW using common distributions; modified from Leos-Barajas & Michelot (2018) arXiv pre-print (https://arxiv.org/pdf/1806.10639.pdf)

HMM.sim = function(nsim, nobs, SL.params, TA.params, Z0) {  
  #nsim=number of simulations
  #nobs=number of observations per simulation (vector)
  #SL.params=df of shape and rate params
  #TA.params=df of mean TA and concen. param
  #Z0=coords of initial location
  
  #uses gamma and wrapped cauchy distribs
  #behaviors params must be in order
  #for simulating w/ 3 behavioral states
  
  
  track.list<- list()  #store all tracks of all durations
  max.cnt<- table(nobs) %>% max()
  sim.ID<- rep(1:max.cnt, nsim/max.cnt)  #running 5 different simulations per track length
  
  for (j in 1:nsim) {
      
    # Number of states
    N <- 3
    # transition probabilities
    Gamma <- matrix(c(0.9, 0.05, 0.05,
                      0.05, 0.9, 0.05,
                      0.05, 0.05, 0.9),
                    nrow = 3, ncol = 3)
    # initial distribution set to the stationary distribution 
    delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))
    
    # state-dependent distribution params
    #step length (gamma)
    gamma.shape <- SL.params[,1]
    gamma.rate <- SL.params[,2]
    
    #turning angle (wrapped Cauchy)
    WC.mean <- TA.params[,1]
    WC.concen <- TA.params[,2]
    
    nobs1 <- nobs[j]
    S <- rep(NA, nobs1) 
    y <- matrix(NA, nobs1, 2)  #2 cols; for step length and turning angle
    
    # initialize state and observation
    S[1] <- sample(1:N, size=1, prob=delta)  #latent state
    y[1,1] <- rgamma(n=1, shape=gamma.shape[S[1]], rate=gamma.rate[S[1]])  #SL
    y[1,2] <- rwrappedcauchy(n=1, mu=circular(WC.mean[S[1]]), rho=WC.concen[S[1]]) %>%  #TA
      ifelse(. > pi, .-(2*pi), .)  #deals with values between pi and 2pi
    
    # simulate state and observation processes forward
    for(t in 2:nobs1) {
      S[t] <- sample(1:N, size=1, prob=Gamma[S[t-1],])  #latent state
      y[t,1] <- rgamma(n=1, shape=gamma.shape[S[t]], rate=gamma.rate[S[t]])  #SL
      y[t,2] <- rwrappedcauchy(n=1, mu=circular(WC.mean[S[t]]), rho=WC.concen[S[t]]) %>%  #TA
        ifelse(. > pi, .-(2*pi), .)  #deals with values between pi and 2pi
    }
    
    # cumulative angle
      Phi <- cumsum(y[,2])
      
      # step length components
      dX <- y[,1]*cos(Phi)
      dY <- y[,1]*sin(Phi)
      
      # actual X-Y values
      X <- c(Z0[1], Z0[1] + cumsum(dX))
      Y <- c(Z0[2], Z0[2] + cumsum(dY))
      
      track<- data.frame(id = as.character(paste0(sim.ID[j],"_",nobs1)),
                         time1 = 1:(nobs1 + 1),
                         x = X,
                         y = Y,
                         SL = c(NA, y[,1]),
                         TA = c(NA, y[,2]),
                         state = c(NA, S),
                         track_length = nobs1)
      
      track.list[[j]]<- track
  }
  
    names(track.list)<- track.list %>% modify_depth(1, ~unique(.$id)) %>% unlist()
    
  
  return(track.list)
}
