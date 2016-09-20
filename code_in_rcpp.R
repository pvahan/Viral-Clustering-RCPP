library(RANN)
library(Rcpp)


setwd("~/Desktop/Virus Clustering/Final Code")
sourceCpp("Spread.cpp")
sourceCpp("Suppress.cpp")



#### Sample Smallest ####
RandSmallest = function(y, n) {
  samp = NULL
  nSort = length(unique(y))
  sortMinFreq = as.numeric(names(sort(table(y))))
  samplength = sapply(1:nSort, function(x) (which(y == sortMinFreq[x])))
  for(sampIter in 1:length(samplength)) {
    if(length(samplength[[sampIter]]) > 1) {
      samp = c(samp, sample(samplength[[sampIter]]))
    } else {
      samp = c(samp, samplength[[sampIter]])
    }
  }
  return(samp)
}


#### Main Function ####
ViralCl = function(data, l = 3, maxK = 50, itr = TRUE, plt = F, maxIter = 1000, 
                   s = 30, eta = 1.2, KTrue = 1) {
  #### Initialization ####
  n = nrow(data); p = ncol(data); K = n; samp = sample(1:n)
  y = NULL; y[[1]] = 1:n; gamma = 1; u = l
  i = 1; delta  = 0; t=0; m = floor(log2(n)); data = as.matrix(data)
  
  knndata = nn2(data, k = m+1)
  minDistNames = t(knndata$nn.idx[, -1])
  sampMat = matrix(sample(1:m, maxIter*n, replace = T), maxIter, n)
  for(jj in 1:n) sampMat[, jj] = minDistNames[sampMat[, jj], jj]
  
  while(K[i] > min(4*maxK, n)) {
    #### Spread Virus ####
    zz = rep(0, n)
    zz[as.numeric(names(sort(table(y[[i]]))))] = sort(table(y[[i]]))
    y[[i+1]] = Spread(y[[i]], samp, sampMat[i,], zz, K[i], min(4*maxK, n), n)
    K[i+1] = length(unique(y[[i+1]]))
    samp = RandSmallest(y[[i+1]], n)
    i = i+1
  }
  
  
  gamma[i+1] = 1; t[i] = n; t[i+1] = n
  while(gamma[i+1] > 10^-6) {
    
    if(u > 0) {
      
      #### Spread Virus ####
      samp = RandSmallest(y[[i]], n)
      zz = rep(0, n)
      zz[as.numeric(names(sort(table(y[[i]]))))] = sort(table(y[[i]]))
      y[[i+1]] = Spread(y[[i]], samp, sampMat[i,], zz, K[i], 1, n)
      u = u - 1
      
    } else {
      #### Suppress Virus ####
      y[[i+1]] = Suppress(data, y[[i]])
      u = l 
    }
    i = i+1
    K[i] = length(unique(y[[i]]))
    
    
    if (i > 1) delta[i] = mean(y[[i]] != y[[i-1]]) # Calculate Delta
    if (i%%s == 0 &&  gamma[i-s+1] < gamma[i]) {   # Update t 
      t[i+1] = round(t[i]/eta)
    } else {
      t[i+1] =  t[i]
    }
    
    if (delta[i] > K[i]/t[i+1]) {                  # Update gamma
      gamma[i+1] = gamma[i]*(1+delta[i])
    } else {
      gamma[i+1] = gamma[i]/2
    }
    if(itr) print(c(K[i], gamma[i+1], i))         # print some statistics
    if (K[i] <= KTrue) break
    
  }
  
  
  # perform Suppress Virus until convergence
  while(sum(y[[i-1]] != y[[i]]) != 0) {
    y[[i+1]] = Suppress(data, y[[i]])
    i = i+1
    K[i] = length(unique(y[[i]]))
  }
  
  return(y[[i]])
}


outvc = ViralCl(iris[, -5])
