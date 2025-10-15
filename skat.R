library(CompQuadForm)

skat.full <- function(img,                                                                                    # dataframe including response, covariates, and imaging data (one row for each sample)
                 atlas=F,                                                                                     # dataframe identifying the ROI for each vertex (one row for each vertex)
                 response, 
                 model='gaussian',                                                                            # Gaussian or binomial
                 covariates=c('Age', 'Sex')) {
  
  if (length(atlas) > 1) {img1 <- img}
  else {
    img1 <- img %>% 
      dplyr::select(-setdiff(names(.)[grep('^R[0-9]', names(.))],                                             
                             paste0('R', atlas$vertex)))
  }
  
  G <- as.matrix(img1[grep('^R[0-9]', names(img1))])                                                          # extract imaging data matrix
  Y <- img1[[response]]
  
  if (length(covariates) > 0) {
    null.model <- glm(as.formula(paste0(response, '~', paste(covariates, collapse='+'))),                     # fit covariate-only model
                      family=model, 
                      data=img1)
    
    Y.hat <- null.model$fitted.values
    X <- model.matrix(null.model)
    m <- ncol(X)
  } else {
    Y.hat <- rep(mean(Y), length(Y))
    X <- rep(1, length(Y))
    m <- 1
  }
  
  skat <- norm(t(G) %*% (Y - Y.hat), type='2')^2                                                              # compute SKAT statistic (with identity kernel)
  
  if (model == 'gaussian') {
    V <- norm(Y - Y.hat, type='2')^2 / (length(Y) - m) * diag(length(Y))                                      # compute V, see [Lin et al. pg 91]
  } else if (model == 'binomial') {
    V <- diag(Y.hat * (1 - Y.hat))
  }
  
  P0 <- V - V %*% X %*% solve(t(X) %*% V %*% X) %*% t(X) %*% V                                                # compute P0, see [Lin et al. pg 91]
  eig <- eigen(P0, symmetric=T)
  indices <- which(eig$values > 0)                                                                            # omit very small negative eigenvalues due to numerical errors
  P0.sqrt <- eig$vectors[, indices] %*% diag(sqrt(eig$values[indices])) %*% t(eig$vectors[, indices])
  
  p.value <- davies(skat, svd(P0.sqrt %*% G, nv=0)$d^2)$Qq                                                    # compute p-value according to eq (6) in [Lin et al.]
  
  return(p.value)
}