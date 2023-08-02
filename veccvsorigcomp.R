source("aux_functions.R")

replicates = 25

p = 300; nsamp = 10 #we assume that the data is from R^p
l_f_matrix = rbind(c(0.2, 0.5), c(0.2, 0.5)); 
sig_f_matrix = rbind(c(0.2, 0.5), c(0.5, 0.2)) 
x_grid = seq(-6, 6, length.out = p)
mvecchia = c(floor(p/20), floor(p/12), floor(p/10))

clusterdata_full = array(dim = c(nrow(sig_f_matrix), replicates, 2 * nsamp))
clustertime = array(dim = c(nrow(sig_f_matrix), replicates))

clusterdata_veccs = array(dim = c(nrow(sig_f_matrix), replicates, 
                                  length(mvecchia), 2 * nsamp))
clustertime_veccs = array(dim = c(nrow(sig_f_matrix), replicates, 
                                  length(mvecchia)))
for(i in 1:nrow(sig_f_matrix)){
  
  #producing original data
  l_f = l_f_matrix[i, ]; sig_f = sig_f_matrix[i, ]
  
  X = x_grid / (max(x_grid) - min(x_grid))
  cov_matrix_1 = sig_f[1]^2 * exp( - (fields::rdist(X))^2/ l_f[1]^2) + diag(10^-10, length(x_grid))
  cov_matrix_2 = sig_f[2]^2 * exp( - (fields::rdist(X))^2/ l_f[2]^2) + diag(10^-10, length(x_grid))
  cov1.chol = t(chol(cov_matrix_1))
  cov2.chol = t(chol(cov_matrix_2))
  
  for(j in 1:replicates){
    
    set.seed(1234 * i + 100 * j)
    normal_samples = matrix(rnorm(2*nsamp*p), nrow = p)
    f1_samps = cov1.chol %*% normal_samples[ , 1:nsamp]
    f2_samps = cov2.chol %*% normal_samples[ , (nsamp + 1):(2*nsamp)]
    sample_collection = Y = t(cbind(f1_samps, f2_samps))
    
    proc.time() -> start.time
    cl_data = gmm.fromscratch.v2(X = X, 
                                 Y = sample_collection, 
                                 k = 2, itern_em = 100)
    end.time <- proc.time()
    
    clusterdata_full[i, j, ] = cl_data$cluster
    clustertime[i, j] = as.numeric((end.time - start.time)[3])
    
    for(k in 1:length(mvecchia)){
      
      m = mvecchia[k]
      orgordering = 1:p; locs = 1:p/p 
      vecchialocspecs = getLocsVecchia(locs = matrix(locs), 
                                       ordering = "maxmin", m = m, seed = 1234)
      
      start.time <- proc.time()
      cl_data_vecchia = gmm.fromscratch.v2.vecchia(X = X[vecchialocspecs$MatchOrigtoPerms], Y = sample_collection[ , vecchialocspecs$MatchOrigtoPerms], 
                                                   k = 2, itern_em = 100, vecchia_obj = vecchialocspecs$vecchia_obj)
      end.time <- proc.time()
      clustertime_veccs[i, j, k] = as.numeric((end.time - start.time)[3])
      clusterdata_veccs[i, j, k, ] = cl_data_vecchia$cluster
      cat(paste("\ti = ", i, ", j = ", j, ", k =", k, ".\n"))
    }
  }
}

save(clustertime_veccs, clusterdata_veccs, clusterdata_full, clustertime,
     file = "origvsvecchia.Rdata")

misclasfull = array(dim = c(2, replicates))
for(i in 1:2){
  for(j in 1:replicates){
    
    misclasfull[i, j] = ratemisclass(clusterA = clusterdata_full[i, j, ], 
                                     clusterB = rep(1:2, each = nsamp))
  }
}

misclasvecchia = array(dim = c(2, replicates, length(mvecchia)))

for(i in 1:2){
  
  for(j in 1:replicates){
    
    for(k in 1:length(mvecchia)){
      
      misclasvecchia[i, j, k] = ratemisclass(clusterA = clusterdata_veccs[i, j, k, ],
                                            clusterB = rep(1:2, each = nsamp))
    }
  }
}

groups = as.factor(rep(1:4, each = replicates))

y = c(misclasfull[1, ],c( misclasvecchia[1, , ]))

pdf("boxplot_case1.pdf")
boxplot(y ~ groups, main = "Original vs Vecchia")
dev.off()

y = c(misclasfull[2, ],c( misclasvecchia[2, , ]))

pdf("boxplot_case2.pdf")
boxplot(y ~ groups, main = "Original vs Vecchia")
dev.off()