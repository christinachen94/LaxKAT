library(tidyverse)
library(parallel)

tau2 <- as.numeric(commandArgs(TRUE)[[1]])
beta.shape <- commandArgs(TRUE)[[2]]

signal.regions <- c(35, 36, 37, 38, 39)
shape <- ''

if (beta.shape == 'normal') {shape <- ''} else if (beta.shape == 'flat') {shape <- '_flat'}

dir <- file.path('~/LaxKAT/revision', paste0('dk_binary_', paste(signal.regions, collapse='_'), shape))

atlas <- read.table('~/LaxKAT/dk_right_atlas.txt', header=T)

img <- readRDS('/home/c13/Data/n839_thickness_demog.rds') %>%                                                 # load FreeSurfer data
  filter(Diagnosis == 'MCI' & !is.na(ADNI_MEM)) %>%                                                           # filter to patients w/ MCI diagnosis & memory score available
  dplyr::select(!matches('^L[0-9]')) %>%                                                                      # remove left hemisphere vertices
  dplyr::select(!(matches('^R[0-9]') & where(~ any(. == 0, na.rm=T)))) %>%                                    # remove right hemisphere vertices with 0/NA
  dplyr::select(-setdiff(names(.)[grep('^R[0-9]', names(.))],
                         paste0('R', atlas$vertex)))

labelvec <- atlas$label[match(names(img)[grep('^R[0-9]', names(img))],
                              paste0('R', atlas$vertex))]

signal.indices <- which(labelvec %in% signal.regions)

compute_naive_p_values <- function(j, img, betas, responses) {
  
  beta <- betas[, j]
  response <- responses[, j]
  
  G <- img[grep('^R[0-9]', names(img))]
  G[, signal.indices] <- G[, signal.indices] + outer(response, beta[signal.indices]) 
  
  cortical_thickness <- as.data.frame(t(as.matrix(G)))
  cortical_thickness$label <- labelvec
  
  roi_averages <- as.data.frame(t(as.matrix(cortical_thickness %>%
                                              group_by(label) %>%
                                              summarize(across(everything(), mean, na.rm=T)) %>%
                                              dplyr::select(-label))))
  
  roi_averages$Sex <- img$Sex
  roi_averages$Age <- img$Age
  roi_averages$ADNI_MEM <- response
  
  sapply(1:34, function(x) summary(glm(as.formula(paste0('ADNI_MEM ~ Age + Sex + V', toString(x))),
                                       family=binomial(link='logit'), 
                                       data=roi_averages))$coefficients[paste0('V', toString(x)), 'Pr(>|z|)'])
}

responses <- read.table(file.path(dir, paste0('response_', format(tau2, scientific=F), '.txt')))
  
betas <- read.table(file.path(dir, paste0('beta_', format(tau2, scientific=F), '.txt')))
  
naive_p <- do.call(cbind, mclapply(1:1000, compute_naive_p_values, img, betas, responses, 
                                   mc.cores=as.numeric(Sys.getenv('LSB_DJOB_NUMPROC'))))
  
write.table(naive_p, file.path(dir, paste0('naive_p_values_unadjusted_', format(tau2, scientific=F), '.txt')),
            row.names=F,
            col.names=F)
  

