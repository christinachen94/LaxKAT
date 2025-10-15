library(tidyverse)

# assume that all vertices in image are present in atlas
simulate.beta.normal <- function(img, 
                                 atlas,                     
                                 signal.regions,                                                              # vector of ROI numbers indicating the signal ROIs
                                 tau2) {                                                                      # betas for signal vertices are i.i.d. N(0, tau2)
  
  labelvec <- atlas$label[match(names(img)[grep('^R[0-9]', names(img))],                                      # generate vector of ROI number labels for vertices  
                                paste0('R', atlas$vertex))]  
  
  beta <- rep(0, length(labelvec))
  signal.vertices <- labelvec %in% signal.regions
  beta[signal.vertices] <- rnorm(sum(signal.vertices), 0, sqrt(tau2))
  
  return(beta)
}

# assume that all vertices in image are present in atlas
simulate.beta.flat <- function(img, 
                               atlas, 
                               signal.regions,                                                                # vector of r ROI numbers indicating the signal ROIs
                               magnitude,                                                                     # length-r vector indicating the magnitude of the beta across each signal ROI
                               signs) {                                                                       # length-r vector of +/- 1s indicating the sign of the beta across each signal ROI
  
  labelvec <- atlas$label[match(names(img)[grep('^R[0-9]', names(img))],                                      # generate vector of ROI number labels for vertices  
                                paste0('R', atlas$vertex))]  
  
  beta <- rep(0, length(labelvec))
  
  for (l in 1:length(signal.regions)) {
    signal.vertices <- which(labelvec == signal.regions[l])
    beta[signal.vertices] <- rep(signs[l] * magnitude[l], length(signal.vertices))
  }
  
  return(beta)
}