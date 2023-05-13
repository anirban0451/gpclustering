source("aux_functions.R")


##creating two potential GP structures
p = 100; nsamp = 10 #we assume that the data is from R^p
l_f = c(0.2, 0.5); sig_f = c(0.5, 0.2) 
x_grid = seq(-6, 6, length.out = p)

#f1_mean = cos(x_grid)
#f2_mean = 15 + 2 * sin(x_grid)
#f3_mean = -(10 + cos(x_grid/2) + sin(x_grid))


#defining covariance matrix and creating samples
X = x_grid / (max(x_grid) - min(x_grid))
cov_matrix_1 = sig_f[1]^2 * exp( - (fields::rdist(X))^2/ l_f[1]^2) + diag(10^-10, length(x_grid))
cov_matrix_2 = sig_f[2]^2 * exp( - (fields::rdist(X))^2/ l_f[2]^2) + diag(10^-10, length(x_grid))
cov1.chol = t(chol(cov_matrix_1))
cov2.chol = t(chol(cov_matrix_2))

set.seed(1234)
normal_samples = matrix(rnorm(2*nsamp*p), nrow = p)
f1_samps = cov1.chol %*% normal_samples[ , 1:nsamp]
f2_samps = cov2.chol %*% normal_samples[ , (nsamp + 1):(2*nsamp)]
sample_collection = Y = t(cbind(f1_samps, f2_samps))



##checking plots to ensure the samples are drawn properly-----
matplot(x = x_grid, y = t(sample_collection), 
        type = "l", lty = rep(1:2, each =nsamp), 
        col = rep(c("red", "black"), each = nsamp))

##gpclustering
proc.time() -> start.time
cl_data = gmm.fromscratch.v2(X = X, Y = sample_collection, k = 2, itern_em = 300)
end.time <- proc.time()

end.time - start.time

##vecchia setup--------------------------------------
seed = 1234
set.seed(seed); m = p
#ordered_indices = sample.int(p, size = p) 
#needed for vecchia, in general no fixed ordering
orgordering = 1:p; locs = 1:p/p 
#original ordering, which is basically natural ordering
#MatchOrigtoPerms = match(ordered_indices, orgordering)
#MatchPermstoOrig = order(ordered_indices) 
#writing ordered_indices[MatchPermstoOrig] gives original ordering
vecchialocspecs = getLocsVecchia(locs = matrix(locs), ordering = "maxmin", m = m, seed = seed)

#vecchia_obj = vecchia_specify_modified(locs = matrix(ordered_indices), ordering = "none", 
#                                       conditioning = "mra", m = m)

start.time <- proc.time()
cl_data = gmm.fromscratch.v2(X = X[vecchialocspecs$MatchOrigtoPerms], Y = sample_collection[ , vecchialocspecs$MatchOrigtoPerms], 
                             k = 2, itern_em = 300)
cl_data_vecchia = gmm.fromscratch.v2.vecchia(X = X[vecchialocspecs$MatchOrigtoPerms], Y = sample_collection[ , vecchialocspecs$MatchOrigtoPerms], 
                                     k = 2, itern_em = 300, vecchia_obj = vecchialocspecs$vecchia_obj)
end.time <- proc.time()
end.time - start.time

clustering = gmm.fromscratch.vecchia(X = sample_collection[ , MatchOrigtoPerms], k = 3,
                                     vecchia_obj = vecchia_obj)

mapped_cluster = clustering$cluster

mapped_cluster2 = mapped_cluster
mapped_cluster2[mapped_cluster == 3] = 1
mapped_cluster2[mapped_cluster == 1] = 2
mapped_cluster2[mapped_cluster == 2] = 3

cluster_misspecification_index = which(mapped_cluster2 != rep(1:3, each = 50))
matplot(x = x_grid, y = cbind(f1_mean, f2_mean, f3_mean, 
                        sample_collection[cluster_misspecification_index, ]), type = "l",
        col = c("red", "yellow", "green", rep("black", length(cluster_misspecification_index))))