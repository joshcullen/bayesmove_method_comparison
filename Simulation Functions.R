### Correlated random walk (CRW) model

CRW.sim=function(nsim, ntseg, nstep, SL.params, TA.params, Z0) {  
  #nsim=number of simulations
  #ntseg=number of time segments
  #nstep=number of steps per time segment
  #nbehav=number of behavioral states ==> derive from nrow SL.params
  #SL.params=df of shape and scale params for step lengths
  #TA.params=df of mean turning angle and concentration param
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
      
      
      #create vector of behaviors within each time segment
      behav.full<- list()
      for (i in 1:length(behav[[j]])) {
        if (i %in% pure) {
          behav.full[[i]]<- rep(behav[[j]][i], nstep)
        }else if (behav[[j]][i] == 1) {
          behav.full[[i]]<- rmultinom(nstep, 1, c(0.8, 0.1, 0.1))  %>% 
            apply(., 2, function(x) which(x == 1))
        } else if (behav[[j]][i] == 2) {
          behav.full[[i]]<- rmultinom(nstep, 1, c(0.1, 0.8, 0.1))  %>% 
            apply(., 2, function(x) which(x == 1))
        } else if (behav[[j]][i] == 3) {
          behav.full[[i]]<- rmultinom(nstep, 1, c(0.1, 0.1, 0.8))  %>% 
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
      track<- data.frame(id = as.character(paste0(j,"_",k)), x = X, y = Y, SL = c(NA,SL),
                         TA = c(NA, TA),
                         behav_fine = as.factor(c(NA, behav.all[[j]])),
                         behav_coarse = c(NA, rep(behav[[j]], each=nstep)) %>% factor(),
                         track_length = length(behav.all[[j]]))
      
      sim.list[[j]]<- track
      
      sim.brkpts[[j]]<- which(diff(behav[[j]]) != 0) * nstep
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

