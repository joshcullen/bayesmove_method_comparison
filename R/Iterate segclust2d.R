
#for mapping segclust function over a list

run.segclust2d_map = function(list1) {
  
  set.seed(123)
  
  # print(paste("Running track", names(list1)[i]))
  start.time<- Sys.time()
  
  Kmax1<- 1.5 * (list1$track_length[1] / 100)
  
  # to account for different track lengths while maintaining same 'lmin
  if (nrow(list1) <= 1000) {
    tmp<- segclust(list1, Kmax = Kmax1, lmin=50, ncluster = c(2,3,4),
                   seg.var = c("SL","abs.TA"), scale.variable = FALSE, subsample = FALSE)
  } else {
    tmp<- segclust(list1, Kmax = Kmax1, lmin=50, ncluster = c(2,3,4),
                   seg.var = c("SL","abs.TA"), scale.variable = FALSE, subsample = TRUE)
  }
  
  end.time<- Sys.time()
  elapsed.time<- difftime(end.time, start.time, units = "min")
  p()
  
  list(models = tmp, elapsed.time = elapsed.time)
  
}

#------------------------------------

# for extracting the mean and sd for the proper model w/ highesr BIC by ncluster
extract.segclust2d.behav.params = function(list1, ncluster) {
  nseg <- list1$Kopt.BIC[ncluster]
  dat.out<- pluck(list1$outputs[[paste(ncluster, "class -", nseg, "segments")]], "states")
  
  dat.out
}
