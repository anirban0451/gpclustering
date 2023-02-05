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
      mvn.c <- t(apply(mvn.c, 1, FUN = function(x) exp(x - (0.5 * (max(x) + min(x))))))
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
    
    val = log(diag(U)) - 0.5 *as.numeric(crossprod(crossprod(U, matrix(xi-mu))))
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