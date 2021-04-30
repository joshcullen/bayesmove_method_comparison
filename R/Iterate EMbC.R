
#for mapping EMbC model over a list

run.embc = function(list1) {
  
  set.seed(1)
  
  start.time<- Sys.time()
  
  # run EMbC model
  mod<- embc(as.matrix(list1), info=-1)
  
  end.time<- Sys.time()
  elapsed.time<- difftime(end.time, start.time, units = "min")
  
  
  p()
  
  list(model = mod, elapsed.time = elapsed.time)
  
}

#--------------------------------------

extract.embc.params = function(mod) {
  
  #extract means of velocity and turning angle
  mu<- map(mod@P, ~pluck(., "M")) %>% 
    set_names(c('LL','LH','HL','HH')) %>% 
    bind_rows(.id = "state") %>% 
    t() %>% 
    data.frame() %>% 
    rename(SL.mu = X1, absTA.mu = X2)
  
  #extract SDs of velocity and turning angle
  sig<- map(mod@P, ~pluck(., "S") %>% 
              diag()) %>% 
    set_names(c('LL','LH','HL','HH')) %>% 
    bind_rows(.id = "state") %>% 
    t() %>% 
    data.frame() %>% 
    rename(SL.sd = X1, absTA.sd = X2)
  
  #combine into single DF
  tmp<- cbind(mu, sig)
  
  tmp
}
