Rfast::sourceR(path = "GPvecchia_mod//R//")
library(GPvecchia); library(fields)
km.fromscratch <- function(X, k){
  p <- ncol(X)  # number of parameters
  n <- nrow(X)  # number of observations
  Delta <- 1; iter <- 0; itermax <- 30
  while(Delta > 1e-4 && iter <= itermax){
    # initiation
    if(iter == 0){
      centroid <- X[sample(nrow(X), k),]
      centroid_mem <- centroid
    }
    
    # equivalent to E-step
    d <- sapply(1:k, function(c) sapply(1:n, 
                                        function(i) sum((centroid[c,] - X[i,])^2) ))
    cluster <- apply(d, 1, which.min)
    
    # equivalent to M-step
    centroid <- t(sapply(1:k, function(c) 
      apply(X[cluster == c,], 2, mean)))
    
    Delta <- sum((centroid - centroid_mem)^2)
    iter <- iter + 1; centroid_mem <- centroid
  }
  return(list(centroid = centroid, cluster = cluster))
}

# Uses EM algorithm with multivariate normal
# distribution to estimate cluster probability
mvnorm.cov.inv <- function(Sigma) {
  # Eigendecomposition of covariance matrix
  E <- eigen(Sigma)
  Lambda.inv <- diag(E$values^-1)   # diagonal matrix
  Q <- E$vectors
  return(Q %*% Lambda.inv %*% t(Q))
}

mvnorm.cov.inv.dup <- function(Sigma) {
  # Eigendecomposition of covariance matrix
  E <- eigen(Sigma)
  Lambda.inv.sqrt <- diag(E$values^-0.5)   # square root of diagonal matrix
  Q <- E$vectors
  return(tcrossprod(tcrossprod(Q, Lambda.inv.sqrt)))
}


#multivariate Gaussian pdf
mvn.pdf.i <- function(xi, mu, InvSigma, logval = TRUE){
  
  if(logval == FALSE){
    
    return(sqrt(det(InvSigma)/ (2*pi)^length(xi)) * 
             exp(-(1/2) * t(xi - mu) %*% InvSigma 
                 %*% (xi - mu)  )) 
  }else{
    
    return(log(sqrt(det(InvSigma))) - (1/2) * t(xi - mu) %*% InvSigma 
           %*% (xi - mu))
  }
}
mvn.pdf <- function(X, mu, Sigma, logval = TRUE){
  
  InvSigma = mvnorm.cov.inv.dup(Sigma)
  return(apply(X, 1, function(xi) mvn.pdf.i(as.numeric(xi), mu, InvSigma, logval = logval)))
}
  
gmm.fromscratch <- function(X, k, logProbs = FALSE){
  p <- ncol(X)  # number of parameters
  n <- nrow(X)  # number of observations
  Delta <- 1; iter <- 0; itermax <- 30
  while(Delta > 1e-4 && iter <= itermax){
    # initiation
    if(iter == 0){
      km.init <- km.fromscratch(X, k)
      mu <- km.init$centroid; mu_mem <- mu
      w <- sapply(1:k, function(i) length(which(km.init$cluster == i)))
      w <- w/sum(w)
      cov <- array(dim = c(p, p, k))
      for(i in 1:p) for(j in 1:p) for(c in 1:k) cov[i, j, c] <- 
        1/n * sum((X[km.init$cluster == c, i] - mu[c, i]) *
                    (X[km.init$cluster == c, j] - mu[c, j]))
    }
    
    # E-step
    if(logProbs == FALSE){
      
      mvn.c <- sapply(1:k, function(c) mvn.pdf(X, mu[c,], cov[,, c], logval = logProbs))
    }else{
      
      mvn.c <- sapply(1:k, function(c) mvn.pdf(X, mu[c,], cov[,, c], logval = logProbs))
      mvn.c <- t(apply(mvn.c, 1, FUN = function(x) exp(x - max(x))))
    }
    r_ic <- t(w*t(mvn.c)) / rowSums(t(w*t(mvn.c)))
    
    # M-step
    n_c <- colSums(r_ic)
    w <- n_c/sum(n_c)
    mu <- t(sapply(1:k, function(c) 1/n_c[c] * colSums(r_ic[, c] *
                                                         X)))
    for(i in 1:p) for(j in 1:p) for(c in 1:k) cov[i, j, c] <-
      1/n_c[c] * sum(r_ic[, c] * (X[, i] - mu[c, i]) * r_ic[, c] *
                       (X[, j] - mu[c, j]))
    Delta <- sum((mu - mu_mem)^2)
    iter <- iter + 1; mu_mem <- mu
  }
  return(list(softcluster = r_ic, cluster = apply(r_ic, 1,
                                                  which.max)))
}


getLocsVecchia = function(locs, ordering, m, seed){
  
  p = nrow(locs)
  orgordering = 1:p
  if(!isTRUE(ordering %in% c("random", "coord", "maxmin"))){
    cat("\tOnly random, coordinate, or maxmin ordering are allowed at this moment.\n")
  }
  if(ordering == "random"){
    if(missing(seed)){
      set.seed(1234)
    }
    ordered_indices = sample.int(p, size = p)
  }else{
    
    vecchia_obj = vecchia_specify_modified(locs = locs, ordering = ordering, 
                                            conditioning = "mra", m = m)
    ordered_indices = vecchia_obj$ord
    remove(vecchia_obj)
  }
  
  vecchia_obj = vecchia_specify_modified(locs = locs[ordered_indices, ], ordering = none, 
                                         conditioning = "mra", m = m)
  MatchOrigtoPerms = match(ordered_indices, orgordering)
  MatchPermstoOrig = order(ordered_indices)
  
  locs.objects = list()
  locs.objects$vecchia_obj = vecchia_obj
  locs.objects$orgordering = 1:p
  locs.objects$MatchOrigtoPerms = MatchOrigtoPerms
  locs.objects$MatchPermstoOrig = MatchPermstoOrig
  return(locs.objects)
}


getU = function(vecchia_obj, covmatrix){
  
  if(isTRUE(missing(vecchia_obj) || (!is.list(vecchia_obj)))){
    
    stop("The variable vecchia_obj is a list of setups which must be provided
         in order to compute sparse Cholesky factor")
  }
  sig.sel = getMatCov(vecchia_obj, cov_matrix)
  inds = Filter(function(i) !is.na(i), as.vector(t(vecchia_obj$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(vecchia_obj$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(sig.sel)))
  vals = createUcppM(ptrs, inds, cov.vals)
  LMatrix = Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  UMatrix = as.matrix(Matrix::triu(Matrix::t(Matrix::solve(LMatrix, sparse = TRUE))))
  return(UMatrix)
}


mvn.pdf.i.ichol <- function(xi, mu, U, logval = TRUE){
  
  if(logval == FALSE){
    val = as.numeric((1/(sqrt( (2*pi)^length(xi))) * prod(diag(U))) * 
                       exp( - 0.5 * crossprod(crossprod(U, matrix(xi-mu)))))#exp(-(1/2) * t(xi - mu) %*% mvnorm.cov.inv.dup(Sigma) %*% (xi - mu)  ) 
  }else if(logval == TRUE){
    
    val = sum(log(diag(U))) - 0.5 *as.numeric(crossprod(crossprod(U, matrix(xi-mu))))
  }
  
  return(val)
}


mvn.pdf.vecchia <-function(X, mu, sigma, vecchia_obj, logval = TRUE){
  
  Umatrix = getU(vecchia_obj = vecchia_obj, covmatrix = cov_matrix)
  values = apply(X, 1, function(xi) mvn.pdf.i.ichol(as.numeric(xi), mu, Umatrix, logval = logval))
  return(values)
}


gmm.fromscratch.vecchia <- function(X, k, vecchia_obj, logProbs = TRUE){
  p <- ncol(X)  # number of parameters
  n <- nrow(X)  # number of observations
  Delta <- 1; iter <- 0; itermax <- 30
  while(isTRUE(Delta > 1e-4 && iter <= itermax)){
    # initiation
    if(iter == 0){
      km.init <- km.fromscratch(X, k)
      mu <- km.init$centroid; mu_mem <- mu
      w <- sapply(1:k, function(i) length(which(km.init$cluster == i)))
      w <- w/sum(w)
      cov <- array(dim = c(p, p, k))
      for(i in 1:p) for(j in 1:p) for(c in 1:k) cov[i, j, c] <- 
        1/n * sum((X[km.init$cluster == c, i] - mu[c, i]) *
                    (X[km.init$cluster == c, j] - mu[c, j]))
    }
    
    # E-step
    mvn.c <- sapply(1:k, function(c) mvn.pdf.vecchia(X, mu[c,], cov[,, c], vecchia_obj, logval = logProbs))
    mvn.c <- t(apply(mvn.c, 1, FUN = function(x) exp(x - max(x))))
    wt_matrix = t(w*t(mvn.c))
    
    r_ic <- wt_matrix / rowSums(wt_matrix)
    
    # M-step
    n_c <- colSums(r_ic)
    w <- n_c/sum(n_c)
    mu <- t(sapply(1:k, function(c) 1/n_c[c] * colSums(r_ic[, c] *
                                                         X)))
    for(i in 1:p) for(j in 1:p) for(c in 1:k) cov[i, j, c] <-
      1/n_c[c] * sum(r_ic[, c] * (X[, i] - mu[c, i]) * r_ic[, c] *
                       (X[, j] - mu[c, j]))
    Delta <- sum((mu - mu_mem)^2)
    iter <- iter + 1; mu_mem <- mu
  }
  return(list(softcluster = r_ic, cluster = apply(r_ic, 1,
                                                  which.max)))
}

#expected conditional log-likelihood
elik = function(Ydata, k, r_ic, cov){
  
  N = nrow(Ydata); p = ncol(Ydata); G = k
  val = 0
  for(j in 1:G){
    
    Sigjinv = Matrix::solve(cov[ , , j])
    for(i in 1:N){
      
      logprobability = p/2 * log(det(Sigjinv) / (2 * pi)) -
        0.5 * (Ydata[i, , drop = FALSE]) %*% Sigjinv %*% t(Ydata[i, , drop = FALSE]) 
      val = val + r_ic[i, j] * logprobability #has to define
    }
  }
  return(val)
}

#derivative of likelihood wrt l
de_dl = function(Ydata, wmat, Sigmaj, l, X, std){
  
  Sigmajinv = (Matrix::solve(Sigmaj))
  alpha_ij =  Sigmajinv %*% Ydata #each column correspond to one observation
  
  #to adjust weight vector here
  qty1 = tcrossprod(alpha_ij[ ,1, drop = FALSE])
  for(i in 2: nrow(Ydata)){
    
    qty1 = qty1 + tcrossprod(alpha_ij[ ,1, drop = FALSE])
  }
  
  qty1 = qty1 - sum(W) * Sigmajinv
  
  qty2 = matrix(0, ncol = nrow(Ydata), nrow = nrow(Ydata))
  
  for(i in 1:nrow(Ydata)){
    for(j in 1:nrow(Ydata)){
      
      qty2[i, j] = std^2 * ((X[i] - X[j])^2/ (l^3)) * 
        exp(- ((X[i] - X[j])^2/ (2 * l^2)))
    }
  }
  val = 0
  for(i in 1:ncol(Ydata)){
    
    val = val + sum((qty1[i, ]) * (qty2[, i]))
  }
  return(val)
}


#derivative of likelihood wrt sigma_f
de_dstdf = function(Ydata, wmat, Sigmaj, l, X, std){
  
  dsitvalues = rdist(matrix(X))
  qty2 = 2 * std * exp(-(distvalues^2 / (2*l^2)))
  Sigmajinv = (Matrix::solve(Sigmaj))
  alpha_ij =  Sigmajinv %*% Ydata #each column correspond to one observation
  
  #to adjust weight vector here
  qty1 = tcrossprod(alpha_ij[ ,1, drop = FALSE])
  for(i in 2: nrow(Ydata)){
    
    qty1 = qty1 + tcrossprod(alpha_ij[ ,1, drop = FALSE])
  }
  
  qty1 = qty1 - sum(W) * Sigmajinv
  
  val = 0
  for(i in 1:ncol(Ydata)){
    
    val = val + sum((qty1[i, ]) * (qty2[, i]))
  }
  return(val)
}

#derivative of likelihood wrt sigma_n, the noise sqrt(variance)
de_dstdn = function(Ydata, wmat, Sigmaj, l, X, stdn){
  
  Sigmajinv = (Matrix::solve(Sigmaj))
  alpha_ij =  Sigmajinv %*% Ydata #each column correspond to one observation
  
  #to adjust weight vector here
  qty1 = tcrossprod(alpha_ij[ ,1, drop = FALSE])
  for(i in 2: nrow(Ydata)){
    
    qty1 = qty1 + tcrossprod(alpha_ij[ ,1, drop = FALSE])
  }
  
  qty1 = qty1 - sum(W) * Sigmajinv
  
  val = stdn * sum(diag(qty1))
  
  return(val)
}

gmm.fromscratch.v2 <- function(X, Y, k, logProbs = TRUE, seed = 1234, itern_em,
                               l_f_interval = seq(0.3, 3, by = 0.1),
                               sig_f_interval = seq(0.5, 3.5, by = 0.1),
                               sig_n_interval = seq(0.5, 3.5, by = 0.1)){
  
  distmat = fields::rdist(matrix(X))
  #this function separates clusters of GP using the length-scale and variance information
  p <- ncol(Y)  # number of parameters
  n <- nrow(Y)  # number of observations
  Delta <- 1; itern <- 0; itermax <- itern_em
  set.seed(seed)
  while(Delta > 1e-4 && itern <= itermax){
    # initiation
    if(itern == 0){
      #km.init <- km.fromscratch(Y, k)
      #mu <- km.init$centroid; mu_mem <- mu
      #w <- sapply(1:k, function(i) length(which(km.init$cluster == i)))
      #w <- w/sum(w)
      #cov <- array(dim = c(p, p, k))
      #for(i in 1:p) for(j in 1:p) for(c in 1:k) cov[i, j, c] <- 
      #  1/n * sum((Y[km.init$cluster == c, i] - mu[c, i]) *
      #              (Y[km.init$cluster == c, j] - mu[c, j]))
      
      w = rep(1/k, k)
      sig_f_old = runif(k, min = 0.8, max = 4) #GP cov parameter
      sig_n_old = runif(k, min = 0.1, max = 1) #noise parameter of GP kernel
      l_f_old = runif(k, min = 0.5, max = 3.5) #length-scale parameter of GP cov
      cov = array(dim = c(dim(distmat), k))
      for(i in 1:k) cov[ , , i] = (sig_f_old[i])^2 * exp( - distmat^2/(2 * l_f_old[i]^2)) + diag(sig_n_old[i]^2, p)
      #this completes the parameter initialization
      l_f_updates = cbind(l_f_old)
      sig_f_updates = cbind(sig_f_old)
      sig_n_updates = cbind(sig_n_old)
      
      
      # E-step
      if(logProbs == FALSE){
        
        mvn.c <- sapply(1:k, function(c) mvn.pdf(Y, rep(0, p), cov[,, c], logval = logProbs))
      }else{
        
        mvn.c <- sapply(1:k, function(c) mvn.pdf(Y, 0, cov[,, c], logval = logProbs))
        mvn.c <- t(apply(mvn.c, 1, FUN = function(x) exp(x - max(x))))
      }
      r_ic <- t(w*t(mvn.c)) / rowSums(t(w*t(mvn.c)))
      
      el_m = elik(Y, k, r_ic = r_ic, cov);
    }else{
      
      el_m = el_mplusone
    }
    
    
    #hyperparameter optimization
    # M-step
    n_c <- colSums(r_ic)
    w <- n_c/sum(n_c)
    #mu <- t(sapply(1:k, function(c) 1/n_c[c] * colSums(r_ic[, c] *
    #                                                     Y)))
    #for(i in 1:p) for(j in 1:p) for(c in 1:k) cov[i, j, c] <-
    #  1/n_c[c] * sum(r_ic[, c] * (Y[, i] - mu[c, i]) * r_ic[, c] *
    #                   (Y[, j] - mu[c, j]))
    
    #el_m = el_mplusone; el_old = el_m
    
    #have to repeat three M steps: l_f, sig_f, sig_n
    
    #doing optimization for l_f
    for(cl in 1:k){
      
      #l_f optimization
      elik_array = c()
      for(h in l_f_interval){
        
        cov[ , , cl] = 
          (sig_f_old[cl])^2 * exp( - distmat^2/(2 * h^2)) + diag(sig_n_old[cl]^2, p)
        elik_array = append(elik_array, elik(Y, k, r_ic, cov))
      }
      l_f_old[cl] = l_f_interval[which.max(elik_array)]
      remove(elik_array)
      
      cov[ , , cl] = 
        (sig_f_old[cl])^2 * exp( - distmat^2/(2 * l_f_old[cl]^2)) + diag(sig_n_old[cl]^2, p)
      
      #sig_f optimization
      elik_array = c()
      for(h in sig_f_interval){
        
        cov[ , , cl] = 
          (h)^2 * exp( - distmat^2/(2 * l_f_old[cl]^2)) + diag(sig_n_old[cl]^2, p)
        elik_array = append(elik_array, elik(Y, k, r_ic, cov))
      }
      sig_f_old[cl] = sig_f_interval[which.max(elik_array)]
      remove(elik_array)
      cov[ , , cl] = 
        (sig_f_old[cl])^2 * exp( - distmat^2/(2 * l_f_old[cl]^2)) + diag(sig_n_old[cl]^2, p)
      
      #sig_n optimization
      elik_array = c()
      for(h in sig_n_interval){
        
        cov[ , , cl] = 
          (sig_f_old[cl])^2 * exp( - distmat^2/(2 * l_f_old[cl]^2)) + diag(h^2, p)
        elik_array = append(elik_array, elik(Y, k, r_ic, cov))
      }
      sig_n_old[cl] = sig_n_interval[which.max(elik_array)]
      remove(elik_array)
      cov[ , , cl] = 
        (sig_f_old[cl])^2 * exp( - distmat^2/(2 * l_f_old[cl]^2)) + diag(sig_n_old[cl]^2, p)
    }
    
    for(cl in 1:k){
      
      cov[ , , cl] = (sig_f_old[cl])^2 * exp( - distmat^2/(2 * l_f_old[cl]^2)) + diag(sig_n_old[cl]^2, p)
    }
    
    # E-step
    if(logProbs == FALSE){
      
      mvn.c <- sapply(1:k, function(c) mvn.pdf(Y, rep(0, p), cov[,, c], logval = logProbs))
    }else{
      
      mvn.c <- sapply(1:k, function(c) mvn.pdf(Y, 0, cov[,, c], logval = logProbs))
      mvn.c <- t(apply(mvn.c, 1, FUN = function(x) exp(x - max(x))))
    }
    r_ic <- t(w*t(mvn.c)) / rowSums(t(w*t(mvn.c)))
    
    el_mplusone = elik(Ydata = Y, r_ic = r_ic, k = k, cov = cov)
    Delta <- (el_m - el_mplusone)^2
    itern <- itern + 1;
  }
  return(list(softcluster = r_ic, 
              cluster = apply(r_ic, 1, which.max),
              l_f = l_f_old,
              sig_f = sig_f_old,
              sig_n = sig_n_old))
}