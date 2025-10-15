library(tidyverse)

source('~/Scripts/LaxKat/laxkat_variations.R')
source('~/Scripts/LaxKat/laxkat_power_setup.R')

signal.regions <- c(35, 36)
j <- as.numeric(Sys.getenv("LSB_JOBINDEX"))
tau2 <- as.numeric(commandArgs(TRUE)[[1]])
beta.shape <- commandArgs(TRUE)[[2]]

sig2 <- 1

atlas <- read.table('~/LaxKAT/dk_right_atlas.txt', header=T)

img <- readRDS('~/Data/n839_thickness_demog.rds') %>%                                                         # load FreeSurfer data
  filter(Diagnosis == 'MCI' & !is.na(ADNI_MEM)) %>%                                                           # filter to patients w/ MCI diagnosis & memory score available
  dplyr::select(!matches('^L[0-9]')) %>%                                                                      # remove left hemisphere vertices
  dplyr::select(!(matches('^R[0-9]') & where(~ any(. == 0, na.rm=TRUE)))) %>%                                 # remove right hemisphere vertices with 0/NA
  dplyr::select(-setdiff(names(.)[grep('^R[0-9]', names(.))], 
                         paste0('R', atlas$vertex)))

labelvec <- atlas$label[match(names(img)[grep('^R[0-9]', names(img))], 
                              paste0('R', atlas$vertex))]  

u <- sort(unique(labelvec))                                                            

img$ADNI_MEM <- rnorm(nrow(img))

if (beta.shape == 'normal') {
  
  beta <- simulate.beta.normal(img, atlas, signal.regions, tau2)
  shape <- ''
  
} else if (beta.shape == 'flat') {
  
  beta <- simulate.beta.flat(img, atlas, signal.regions, rep(tau2, length(signal.regions)), rep(1, length(signal.regions)))
  shape <- '_flat'
}

img[grep('^R[0-9]', names(img))] <- outer(img$ADNI_MEM, beta) + img[grep('^R[0-9]', names(img))]

result <- laxkat(img, atlas, 1000, response='ADNI_MEM', model='gaussian')

write.table(x=result$unnorm.global.p.value,
            file=paste0('~/LaxKAT/paper/desikan_killiany_', 
                        paste(signal.regions, collapse='_'), 
                        shape,
                        '_1000_runs_unassociated/unnormalized_global_p_value_',
                        format(tau2, scientific=F), 
                        '_', 
                        toString(j), 
                        '.txt'), 
            row.names=F, 
            col.names=F)

write.table(x=result$local.p.values,
            file=paste0('~/LaxKAT/paper/desikan_killiany_', 
                        paste(signal.regions, collapse='_'), 
                        shape,
                        '_1000_runs_unassociated/local_p_values_',
                        format(tau2, scientific=F), 
                        '_', 
                        toString(j), 
                        '.txt'), 
            row.names=F, 
            col.names=F)

write.table(img$ADNI_MEM,
            file=paste0('~/LaxKAT/paper/desikan_killiany_', 
                        paste(signal.regions, collapse='_'), 
                        shape,
                        '_1000_runs_unassociated/response_', 
                        format(tau2, scientific=F), 
                        '_', 
                        toString(j),
                        '.txt'),
            row.names=F,
            col.names=F)

write.table(beta,
            file=paste0('~/LaxKAT/paper/desikan_killiany_', 
                        paste(signal.regions, collapse='_'), 
                        shape,
                        '_1000_runs_unassociated/beta_', 
                        format(tau2, scientific=F), 
                        '_', 
                        toString(j), 
                        '.txt'),
            row.names=F,
            col.names=F)
