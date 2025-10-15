pst <- function(img,                                                                                          # dataframe including response, covariates, and imaging data (one row for each sample)
                atlas,                                                                                        # dataframe identifying the ROI for each vertex (one row for each vertex)
                response, 
                model='gaussian',                                                                             # Gaussian or binomial
                covariates=c('Age', 'Sex')) {
  
  img1 <- img %>% 
    dplyr::select(-setdiff(names(.)[grep('^R[0-9]', names(.))],                                               # subset to vertices present in atlas
                           paste0('R', atlas$vertex))) 
  
  labelvec <- atlas$label[match(names(img1)[grep('^R[0-9]', names(img1))],                                    # generate vector of ROI number labels for vertices
                                paste0('R', atlas$vertex))]
  u <- sort(unique(labelvec))
  r <- length(u)                                                                                              # number of ROIs (dimension of subspace L)
  
  G <- as.matrix(img1[grep('^R[0-9]', names(img1))])
  Y <- img1[[response]]
  
  if (length(covariates) > 0) {
    null.model <- glm(as.formula(paste0(response, '~', paste(covariates, collapse='+'))),                       # fit covariate-only model                     
                      family=model,                                                                             
                      data=img1)
    Y.hat <- null.model$fitted.values
    X <- model.matrix(null.model)
    m <- ncol(X)
    n <- nrow(X)
  } else {
    
    Y.hat <- rep(mean(Y), length(Y))
    X <- rep(1, length(Y))
    m <- 1
    n <- length(Y)
  }
  
  scores <- t(G) %*% (Y - Y.hat) / n                                                                          # compute scores, see [Vandekar et al. pg 823]
  Gamma <- diag((as.vector(Y) - Y.hat)^2)                                                                     # define Gamma, see [Vandekar et al. pg 823]
  Gt.Gamma.X <- t(G) %*% Gamma %*% X
  Omega.hat <- (t(G) %*% Gamma %*% G - Gt.Gamma.X %*% solve(t(X) %*% Gamma %*% X) %*% t(Gt.Gamma.X)) / n      # estimate covariance via eq (19) in [Vandekar et al.]
  
  Q <- apply(sapply(u, function(x) labelvec == x), 2, function(x) x / sqrt(sum(x)))                           # the i-th column of the p * r matrix Q equals the indicator vector for the i-th ROI, see [Vandekar et al. pg 822]
  rotated.scores <- t(Q) %*% scores
  pst <- n * t(rotated.scores) %*% solve(t(Q) %*% Omega.hat %*% Q) %*% rotated.scores                         # compute PST statistic via eq (10) in [Vandekar et al.]
  
  if (model == 'gaussian') {
    p.value <- pf((r * (n - m) / pst - r) / (n - m - r), n - m - r, r)                                        # compute p-value via eq (18) in [Vandekar et al.]
  } else if (model == 'binomial') {
    p.value <- pchisq(pst, r, lower.tail=F)                                                                   # compute p-value via chi2 approximation in [Vandekar et al. pg 822]
  }
  
  return(p.value)                                             
}