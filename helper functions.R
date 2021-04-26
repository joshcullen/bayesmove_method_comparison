
brkpt.accuracy=function(model.brkpts, true.brkpts, acc.tol, dup.tol, miss.tol) {
  
  model.brkpts<- as.numeric(model.brkpts) %>% na.omit() %>% as.vector()
  true.brkpts<- as.numeric(true.brkpts) %>% na.omit() %>% as.vector()
  all.brkpts<- data.frame(brks = c(true.brkpts, model.brkpts),
                          type = rep(c("True","Model"),c(length(true.brkpts),
                                                         length(model.brkpts))))
  
  accuracy<- matrix(NA, length(model.brkpts), 1)
  for (i in 1:length(model.brkpts)) {  #assign brkpts as accurate or inaccurate
    
    tmp<- c(model.brkpts[i] - (acc.tol:0), model.brkpts[i] + (1:acc.tol)) %in% true.brkpts %>%
      sum()
    
    if (tmp == 0) {
      accuracy[i]<- "Inaccurate"
    } else {
      accuracy[i]<- "Accurate"
    }
  }
  
  if (sum(abs(diff(model.brkpts)) <= dup.tol) >= 0) {  #identify duplicate brkpts
    ind<- which(abs(diff(model.brkpts)) <= dup.tol)
    ind<- sort(c(ind, ind+1))
  }
  
  ind.acc<- ind[which(accuracy[ind] == "Accurate")]
  ind.inacc<- ind[which(accuracy[ind] == "Inaccurate")]
  accuracy[ind.acc]<- "Accurate Duplicate"
  accuracy[ind.inacc]<- "Inaccurate Duplicate"
  accuracy<- c(rep("True",length(true.brkpts)), accuracy)
  
  
  #identify missing breakpoints from model
  status<- matrix(NA,length(true.brkpts),1)
  for (i in 1:length(true.brkpts)) {
    
    tmp<- c(true.brkpts[i] - (miss.tol:0), true.brkpts[i] + (1:miss.tol)) %in% model.brkpts %>%
      sum()
    
    if (tmp == 0) {
      status[i]<- "Missing"
    } else {
      status[i]<- "Present"
    }
  }
  
  miss.ind<- which(status =="Missing")
   
  
  all.brkpts$acc<- accuracy
  all.brkpts[miss.ind, "acc"]<- "Missing"
  
  all.brkpts
}

#--------------------------------------

extract.behav.props=function(params, lims, behav.names){  
  #only for gamma and wrapped Cauchy distributions
  #only defined for step lengths ('dist') and turning angles ('rel.angle')
  #params must be alist of data frames storing 2 cols for the params of the gamma and wrapped cauchy distributions
  #lims must be a list of numeric vectors
  
  #number of bins for both params are already defined as 5 (SL) and 8 (TA)
  #order of states is set as encamped, ARS, transit
  
  #Extract bin proportions for step lengths
  props.SL<- list() #create list to store results per behavior
  
  SL<- params[[1]]  #extract SL params from argument; must be first
  dist.bin.lims<- lims[[1]]  #extract limits for SL from argument; must be first
  
  for (j in 1:nrow(SL)) {
    shape1=SL[j,1]
    rate1=SL[j,2]
    
    bins.estim=rep(NA,(length(dist.bin.lims) - 1))
    
    for (i in 2:length(dist.bin.lims)) {
      if (i-1 == 1) {
        bins.estim[i-1]=pgamma(dist.bin.lims[i],shape=shape1,rate=rate1)
      } else {
        bins.estim[i-1]=pgamma(dist.bin.lims[i],shape=shape1,rate=rate1)-
          pgamma(dist.bin.lims[i-1],shape=shape1,rate=rate1)
      }
    }
    
    props.SL[[j]]<- bins.estim
  }
  
  names(props.SL)<- behav.names
  props.SL1<- props.SL %>% 
    bind_rows() %>% 
    pivot_longer(cols = names(.), names_to = "behav", values_to = "prop") %>% 
    arrange(behav) %>% 
    mutate(var = "Step Length", bin = rep(1:5, length(behav.names)))
  
  
  
  #Extract bin proportions for turning angles
  props.TA<- list()  #create list to store results per behavior
  
  TA<- params[[2]]  #extract TA params from argument; must be second
  angle.bin.lims<- lims[[2]]  #extract limits for TA from argument; must be second
  
  for (j in 1:nrow(TA)) {
    mu1=circular(TA[j,1],type='angles',units='radians')
    rho1=TA[j,2]
    
    bins.estim=rep(NA,(length(angle.bin.lims) - 1))
    
    for (i in 2:length(angle.bin.lims)) {
      pwrappedcauchy=integrate(dwrappedcauchy,
                               angle.bin.lims[i-1],
                               angle.bin.lims[i],
                               mu=mu1,
                               rho=rho1)
      bins.estim[i-1]=pwrappedcauchy$value
    }
    
    props.TA[[j]]<- bins.estim
  }
  
  names(props.TA)<- behav.names
  props.TA1<- props.TA %>% 
    bind_rows() %>% 
    pivot_longer(cols = names(.), names_to = "behav", values_to = "prop") %>% 
    arrange(behav) %>% 
    mutate(var = "Turning Angle", bin = rep(1:8, length(behav.names)))
  
  
  
  #Combine proportions for SL and TA
  props<- rbind(props.SL1, props.TA1)
  props$behav<- factor(props$behav, levels = behav.names)
  props<- props %>% 
    arrange(behav, var)
  
  props
}

#--------------------------------------
extract.behav.props_weird=function(params, lims, behav.names){  
  #only for truncated normal and mixed distributions (beta, uniform, truncated normal);
  #step lengths ('dist') and turning angles ('rel.angle') params must be a list
  #of data frames storing 2 cols for the params of each distribution; lims must be a list
  #of numeric vectors
  
  #number of bins for both params are already defined as 5 (SL) and 8 (TA)
  #order of states is set as encamped, ARS, transit
  
  #Extract bin proportions for step lengths
  props.SL<- list() #create list to store results per behavior
  
  SL<- params[[1]]  #extract SL params from argument; must be first
  dist.bin.lims<- lims[[1]]  #extract limits for SL from argument; must be first
  
  for (j in 1:nrow(SL)) {
    mean1=SL[j,1]
    sd1=SL[j,2]
    
    bins.estim=rep(NA,(length(dist.bin.lims) - 1))
    
    for (i in 2:length(dist.bin.lims)) {
      if (i-1 == 1) {
        bins.estim[i-1]=pnorm(dist.bin.lims[i],mean=mean1,sd=sd1) - 
          pnorm(0,mean=mean1,sd=sd1)
      } else {
        bins.estim[i-1]=pnorm(dist.bin.lims[i],mean=mean1,sd=sd1)-
          pnorm(dist.bin.lims[i-1],mean=mean1,sd=sd1)
      }
    }
    
    props.SL[[j]]<- bins.estim / sum(bins.estim)
  }
  
  names(props.SL)<- behav.names
  props.SL1<- props.SL %>% 
    bind_rows() %>% 
    pivot_longer(cols = names(.), names_to = "behav", values_to = "prop") %>% 
    arrange(behav) %>% 
    mutate(var = "Step Length", bin = rep(1:5, length(behav.names)))
  
  
  
  #Extract bin proportions for turning angles
  props.TA<- list()  #create list to store results per behavior
  
  TA<- params[[2]]  #extract TA params from argument; must be second
  angle.bin.lims<- lims[[2]]  #extract limits for TA from argument; must be second
  
  
  #Encamped state (beta distrib)
  a=TA[1,1]
  b=TA[1,2]
  
  bins.estim=rep(NA,(length(angle.bin.lims) - 1))
  
  for (i in 2:length(angle.bin.lims)) {
    if (i-1 == 1) {
      bins.estim[i-1]=pbeta((angle.bin.lims[i]+pi)/(2*pi),shape1=a,shape2=b) - 
        pbeta(0,shape1=a,shape2=b)
    } else {
      bins.estim[i-1]=pbeta((angle.bin.lims[i]+pi)/(2*pi),shape1=a,shape2=b)-
        pbeta((angle.bin.lims[i-1]+pi)/(2*pi),shape1=a,shape2=b)
    }
  }
  
  props.TA[[1]]<- bins.estim
  
  
  
  #ARS state (uniform distrib)
  lo=TA[2,1]
  up=TA[2,2]
  
  bins.estim=rep(NA,(length(angle.bin.lims) - 1))
  
  for (i in 2:length(angle.bin.lims)) {
    if (i-1 == 1) {
      bins.estim[i-1]=punif(angle.bin.lims[i],min=lo,max=up) - 
        punif(-pi,min=lo,max=up)
    } else {
      bins.estim[i-1]=punif(angle.bin.lims[i],min=lo,max=up)-
        punif(angle.bin.lims[i-1],min=lo,max=up)
    }
  }
  
  props.TA[[2]]<- bins.estim
  
  
  
  #Transit state (truncated normal distrib)
  mean1=TA[3,1]
  sd1=TA[3,2]
  
  bins.estim=rep(NA,(length(angle.bin.lims) - 1))
  
  for (i in 2:length(angle.bin.lims)) {
    if (i-1 == 1) {
      bins.estim[i-1]=pnorm(angle.bin.lims[i],mean=mean1,sd=sd1) - 
        pnorm(-pi,mean=mean1,sd=sd1)
    } else {
      bins.estim[i-1]=pnorm(angle.bin.lims[i],mean=mean1,sd=sd1)-
        pnorm(angle.bin.lims[i-1],mean=mean1,sd=sd1)
    }
  }
  
  props.TA[[3]]<- bins.estim
  
  
  
  
  names(props.TA)<- behav.names
  props.TA1<- props.TA %>% 
    bind_rows() %>% 
    pivot_longer(cols = names(.), names_to = "behav", values_to = "prop") %>% 
    arrange(behav) %>% 
    mutate(var = "Turning Angle", bin = rep(1:8, length(behav.names)))
  
  
  
  #Combine proportions for SL and TA
  props<- rbind(props.SL1, props.TA1)
  props$behav<- factor(props$behav, levels = behav.names)
  props<- props %>% 
    arrange(behav, var)
  
  props
}


#--------------------------------------

extract.behav.props_segclust=function(params, lims, behav.names){  
  #only for segclust2d results that use a normal distribution
  #only defined for step lengths ('dist') and turning angles ('rel.angle')
  #params must be a list of data frames storing 4 cols for the params of the normal distributions (mu.SL, mu.TA, sd.SL, sd.TA)
  #lims must be a list of numeric vectors
  
  #number of bins for both params are already defined as 5 (SL) and 8 (TA)
  #order of states is set as encamped, ARS, transit
  
  #Extract bin proportions for step lengths
  props.SL<- list() #create list to store results per behavior
  
  # SL<- params[[1]]  #extract SL params from argument; must be first
  dist.bin.lims<- lims[[1]]  #extract limits for SL from argument; must be first
  
  for (j in 1:nrow(params)) {
    mu1=params[j,1]
    sd1=params[j,3]
    
    bins.estim=rep(NA,(length(dist.bin.lims) - 1))
    
    for (i in 2:length(dist.bin.lims)) {
      if (i-1 == 1) {
        bins.estim[i-1]=pnorm(dist.bin.lims[i], mean=mu1, sd=sd1)
      } else {
        bins.estim[i-1]=pnorm(dist.bin.lims[i], mean=mu1,sd=sd1)-
          pnorm(dist.bin.lims[i-1], mean=mu1, sd=sd1)
      }
    }
    
    props.SL[[j]]<- bins.estim
  }
  
  names(props.SL)<- behav.names
  props.SL1<- props.SL %>% 
    bind_rows() %>% 
    pivot_longer(cols = names(.), names_to = "behav", values_to = "prop") %>% 
    arrange(behav) %>% 
    mutate(var = "Step Length", bin = rep(1:5, length(behav.names)))
  
  
  
  #Extract bin proportions for turning angles and create as bilaterally symmetrical
  props.TA<- list()  #create list to store results per behavior
  
  angle.bin.lims<- lims[[2]]  #extract limits for TA from argument; must be second
  angle.bin.lims2<- angle.bin.lims[angle.bin.lims >= 0]  #only positive values
  
  for (j in 1:nrow(params)) {
    mu1=params[j,2]
    sd1=params[j,4]
    
    bins.estim=rep(NA,(length(angle.bin.lims2) - 1))
    
    for (i in 2:length(angle.bin.lims2)) {
      if (i-1 == 1) {
        bins.estim[i-1]=pnorm(angle.bin.lims2[i], mean=mu1, sd=sd1)
      } else {
        bins.estim[i-1]=pnorm(angle.bin.lims2[i], mean=mu1,sd=sd1)-
          pnorm(angle.bin.lims2[i-1], mean=mu1, sd=sd1)
      }
    }
    
    props.TA[[j]]<- c(rev(bins.estim), bins.estim)  #create symmetrical distribution
  }
  
  names(props.TA)<- behav.names
  props.TA1<- props.TA %>% 
    bind_rows() %>% 
    pivot_longer(cols = names(.), names_to = "behav", values_to = "prop") %>% 
    arrange(behav) %>% 
    mutate(var = "Turning Angle", bin = rep(1:8, length(behav.names)))
  
  
  
  #Combine proportions for SL and TA
  props<- rbind(props.SL1, props.TA1)
  props$behav<- factor(props$behav, levels = behav.names)
  props<- props %>% 
    arrange(behav, var)
  
  props
}


#--------------------------------------

rtnorm=function(n,lo,hi,mu,sig, log = FALSE){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if (log == TRUE) {
    mu<- exp(mu)
    sig<- exp(sig)
  }
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}

#--------------------------------------

#density of truncated normal distribution
dtnorm=function(x,mean1,sd1,lo,hi){
  denom=pnorm(hi,mean=mean1,sd=sd1)-pnorm(lo,mean=mean1,sd=sd1)
  dnorm(x,mean=mean1,sd=sd1)/denom
}

#--------------------------------------
assign.time.seg_hmm=function(dat, brkpts){  #add tseg assignment to each obs
  tmp=which(unique(dat$id) == brkpts$id)
  breakpt<- brkpts[tmp,-1] %>% purrr::discard(is.na) %>% as.numeric(.[1,])
  breakpt<- breakpt + 1  #only for analysis of HMMs with simulated data
  breakpt1=c(0,breakpt,Inf)
  n=length(breakpt1)
  res=matrix(NA,nrow(dat),1)  #accounts for NA state at beginning of each ID
  for (i in 2:n){
    ind=which(breakpt1[i-1] < dat$time1 & dat$time1 <= breakpt1[i])
    res[ind,]=i-1
  }
  dat$tseg<- as.vector(res)
  dat
}


#--------------------------------------
assign.time.seg_segclust = function(dat){  #add tseg assignment to each obs
  breakpt<- which(diff(dat$state, lag = 1) != 0)
  breakpt1=c(0,breakpt,Inf)
  n=length(breakpt1)
  res=matrix(NA,nrow(dat),1)  #accounts for NA state at beginning of each ID
  for (i in 2:n){
    ind=which(breakpt1[i-1] < dat$time1 & dat$time1 <= breakpt1[i])
    res[ind,]=i-1
  }
  dat$tseg<- as.vector(res)
  dat
}
