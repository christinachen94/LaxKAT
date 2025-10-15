library(parallel)
library(CompQuadForm)
library(data.table)

# skat if the predicted values have already been computed
skat_hat <- function(Y, Y.hat, G, model, inference=T, X=NULL) {
  
  n <- length(Y)
  s <- norm(t(G) %*% (Y - Y.hat), type='2')^2                                                                 # compute SKAT statistic (with identity kernel)
  
  if (inference) {
    m <- max(1, ncol(X))
    
    if (model == 'gaussian') {V <- norm(Y - Y.hat, type='2')^2 / (n - m) * diag(n)}                           # compute V, see [Lin et al. pg 91]
    else if (model == 'binomial') {V <- diag(Y.hat * (1 - Y.hat))}
  
    P0 <- V - V %*% X %*% solve(t(X) %*% V %*% X) %*% t(X) %*% V                                              # compute P0, see [Lin et al. pg 91]
    eig <- eigen(P0, symmetric=T)
    indices <- which(eig$values > 0)                                                                          # omit very small negative eigenvalues due to numerical errors
    P0.sqrt <- eig$vectors[, indices] %*% diag(sqrt(eig$values[indices])) %*% t(eig$vectors[, indices])
  
    coefficients <- svd(P0.sqrt %*% G, nv=0)$d^2
    p.value <- davies(s, coefficients)$Qq                                                                     # compute p-value according to eq (6) in [Lin et al.]
    
    return(list('statistic'=s, 
                'p.value'=p.value))
  } 
  else {
    return(list('statistic'=s, 
                'p.value'=NA))}
}

# function to compute SKAT-like statistic for each ROI
quotient <- function(Y,                                                                                       # length-n vector of responses
                     Y.hat,                                                                                   # length-n vector of predicted responses from covariate-only model
                     G,                                                                                       # matrix of imaging data
                     model,  
                     inference=T,                                                                                    
                     X) {                                                                                     # covariates                                                                                  
  
  skat <- skat_hat(Y, Y.hat, G, model, inference, X)
  
  if (model == 'gaussian') {delta <- rep(1, length(Y))} 
  else if (model == 'binomial') {delta <- as.numeric(Y.hat * (1 - Y.hat))}
  
  if (ncol(G) == 1) {denom <- delta %*% G^2} 
  else {denom <- delta %*% apply(G, 1, function(x) norm(x, type='2')^2)}
  
  if (inference & (model == 'gaussian')) {
    
    return(list('statistic'=skat$statistic / denom, 
                'p.value'=skat$p.value))} 
  else {
    return(list('statistic'=skat$statistic / denom,
                'p.value'=NA))}
}

# function for running permutations in parallel
# returns length-q vector of ROI-specific SKAT-like statistics
run.permutation.quotient <- function(i,                                                                       # index to track permutation number
                                     labelvec,                                                                # vector of ROI number labels for vertices
                                     u,                                                                       # length-q vector of unique ROI labels
                                     Y,                                                                       # length-n vector of responses
                                     Y.hat,                                                                   # length-n vector of predicted responses from covariate-only model
                                     G,                                                                       # n * p matrix of imaging data
                                     model) {                                                                 # Gaussian or binomial 
  
  p <- sample(1:length(Y))                                                                                    # apply the same permutation to responses and predicted responses
  y <- Y[p]
  y.hat <- Y.hat[p]
  
  return(sapply(u, function(x) quotient(y, 
                                        y.hat, 
                                        G[, which(labelvec == x)], 
                                        model, 
                                        inference=F)$statistic))                        
}

laxkat <- function(img, 
                   atlas,                                                                                     # dataframe identifying the ROI for each vertex (one row for each vertex)
                   n.iter=1000,                                                                               # number of permutations for estimating the null distribution
                   response, 
                   model,                                                                                     # Gaussian or binomial
                   covariates=c('Age', 'Sex')) {
  
  img1 <- img %>%                                                                                             # subset to vertices present in atlas
    dplyr::select(-setdiff(names(.)[grep('^R[0-9]', names(.))],               
                           paste0('R', atlas$vertex)))
  
  labelvec <- atlas$label[match(names(img1)[grep('^R[0-9]', names(img1))], paste0('R', atlas$vertex))]        # generate vector of ROI number labels for vertices
  u <- sort(unique(labelvec))
  
  G <- as.matrix(img1[grep('^R[0-9]', names(img1))])
  Y <- img1[[response]]
  n <- length(Y)
  
  if (length(covariates) > 0) {
    null.model <- glm(as.formula(paste0(response, '~', paste(covariates, collapse='+'))),                     # fit covariate-only model
                      family=model, 
                      data=img1)
    
    Y.hat <- null.model$fitted.values
    X <- model.matrix(null.model)
  } 
  
  else {
    Y.hat <- rep(mean(Y), n)
    X <- rep(1, n)
  }
  
  lax.results <- as.data.frame(rbindlist(lapply(u, function (x) quotient(Y, 
                                                                         Y.hat, 
                                                                         G[, which(labelvec == x)], 
                                                                         model, 
                                                                         inference=T, 
                                                                         X))))
  
  lax.statistic <- max(lax.results$statistic)
  
  # obtain table of ROI-specific SKAT statistics across permutations
  # each row corresponds to one ROI
  perm.results <- do.call(cbind, mclapply(1:(n.iter - 1), 
                                          run.permutation.quotient, 
                                          labelvec, 
                                          u, 
                                          Y, 
                                          Y.hat,
                                          G,
                                          model,
                                          mc.cores=as.numeric(Sys.getenv('LSB_DJOB_NUMPROC'))))
  
  # obtain LaxKAT statistic for each permutation (maximum of each column in perm.results)
  lax.null.dist <- apply(perm.results, 2, max) 
  
  unnorm.global.p.value <- (1 + sum(lax.null.dist >= lax.statistic)) / n.iter
  
  local.p.values <-  (1 + rowSums(perm.results >= lax.results$statistic)) / n.iter
  
  return(list('unnorm.global.p.value'=unnorm.global.p.value,
              'local.p.values'=local.p.values,
              'labels'=u))
}

