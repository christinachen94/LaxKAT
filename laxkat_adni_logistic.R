library(tidyverse)

source('~/Scripts/LaxKat/laxkat_variations.R')

atlas <- read.table('~/LaxKAT/dk_right_atlas.txt', header=T)

img <- readRDS('~/Data/n839_thickness_demog.rds') %>%                                                         # load FreeSurfer data
  filter(Diagnosis == 'CN' & (!is.na(ADNI_MEM))) %>% 
  mutate(Sex = as.numeric(Sex == 'M')) %>%
  dplyr::select(!matches('^L[0-9]')) %>%                                                                      # remove left hemisphere vertices
  dplyr::select(!(matches('^R[0-9]') & where(~ any(. == 0, na.rm=TRUE)))) %>%                                 # remove right hemisphere vertices with 0/NA
  dplyr::select(-setdiff(names(.)[grep('^R[0-9]', names(.))], 
                         paste0('R', atlas$vertex)))

img_vertices <- img %>% dplyr::select(matches('^R[0-9]'))

labelvec <- atlas$label[match(names(img_vertices), paste0('R', atlas$vertex))]  
labels <- sort(unique(labelvec))

laxkat.result <- laxkat(img,
                        atlas,
                        10000,
                        response='Sex',
                        model='binomial',
                        covariates=c('Age', 'ADNI_MEM'))

write.table(laxkat.result$local.p.values,
            file='~/LaxKAT/paper/adni_sex_laxkat_p_values_unadjusted_age_score_covariate.txt',
            row.names=F,
            col.names=F)

write.table(laxkat.result$unnorm.global.p.value, 
            file='~/LaxKAT/paper/adni_sex_laxkat_unnormalized_global_p_value_age_score_covariate.txt', 
            row.names=F, 
            col.names=F)

# roi-level marginal tests
roi_averages <- sapply(labels,
                       function(x) rowMeans(img_vertices %>% dplyr::select(which(labelvec == x))))

marginal.p.values <- sapply(1:length(labels), function(x) tail(summary(glm(img$Sex ~ img$Age + img$ADNI_MEM + roi_averages[, x],
                                                               family='binomial'))$coefficients,
                                                   n=1)[, 'Pr(>|z|)'])

write.table(marginal.p.values,
            file='~/LaxKAT/paper/adni_sex_naive_p_values_unadjusted_age_score_covariate.txt',
            row.names=F,
            col.names=F)

which(p.adjust(scan('~/LaxKAT/paper/adni_sex_laxkat_p_values_unadjusted_age_score_covariate.txt', what=numeric(), quiet=T), method='holm') <= 0.05)
which(p.adjust(scan('~/LaxKAT/paper/adni_sex_naive_p_values_unadjusted_age_score_covariate.txt', what=numeric(), quiet=T), method='holm') <= 0.05)
