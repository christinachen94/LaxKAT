setwd('/home/c13/Scripts/LaxKat')

library(parallel)

source('~/Scripts/LaxKat/skat.R')
source('~/Scripts/LaxKat/laxkat_power_setup.R')

tau2.list <- c(0, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01)
# tau2.list <- c(0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01)
signal.regions <- c(35, 36, 37, 38, 39)
shape <- ''
model <- 'binomial'

tau2 <- as.numeric(commandArgs(TRUE)[[1]])
beta.shape <- commandArgs(TRUE)[[2]]

if (beta.shape == 'normal') {shape <- ''} else if (beta.shape == 'flat') {shape <- '_flat'}

dir <- file.path('~/LaxKAT/revision', paste0('dk_binary_', paste(signal.regions, collapse='_'), shape))

atlas <- read.table('~/LaxKAT/dk_right_atlas.txt', header=T)

img <- readRDS('/home/c13/Data/n839_thickness_demog.rds') %>%                                                 # load FreeSurfer data
  filter(Diagnosis == 'MCI' & !is.na(ADNI_MEM)) %>%                                                           # filter to patients w/ MCI diagnosis & memory score available
  dplyr::select(!matches('^L[0-9]')) %>%                                                                      # remove left hemisphere vertices
  dplyr::select(!(matches('^R[0-9]') & where(~ any(. == 0, na.rm=TRUE)))) %>%                                 # remove right hemisphere vertices with 0/NA
  dplyr::select(-setdiff(names(.)[grep('^R[0-9]', names(.))], 
                         paste0('R', atlas$vertex)))

labelvec <- atlas$label[match(names(img)[grep('^R[0-9]', names(img))],
                              paste0('R', atlas$vertex))]

signal.indices <- which(labelvec %in% signal.regions)

skat.parallel <- function(j, betas, responses) {
  
  beta <- betas[, j]
  response <- responses[, j]
  
  img1 <- img[grep('^R[0-9]', names(img))]
  img1[, signal.indices] <- img1[, signal.indices] + outer(response, beta[signal.indices])
  
  img1$Sex <- img$Sex
  img1$Age <- img$Age
  img1$ADNI_MEM <- response
  
  skat.full(img1, atlas, 'ADNI_MEM', model=model, covariates=c('Age', 'Sex'))
}

responses <- read.table(file.path(dir, paste0('response_', format(tau2, scientific=F), '.txt')))
  
betas <- read.table(file.path(dir, paste0('beta_', format(tau2, scientific=F), '.txt')))
  
write.table(unlist(mclapply(1:1000, skat.parallel, betas, responses, mc.cores=as.numeric(Sys.getenv('LSB_DJOB_NUMPROC')))), 
            file=file.path(dir, paste0('skat_', format(tau2, scientific=F), '.txt')), 
            row.names=F, 
            col.names=F)
  
