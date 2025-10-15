library(tidyverse)

# function to create basis based on structural atlas
atlas.table <- function(rois, x) {
  
  read.table(rois[x], 
             header=FALSE, 
             skip=2, 
             col.names=c('vertex', 'x', 'y', 'z', 'label')) %>%
    mutate(label=!!x) %>%                                                                                     # label with ROI number
    mutate(name=gsub('.label', '', basename(rois[!!x])))                                                      # label with ROI name
}

# create basis for Desikan-Killiany atlas
# construct dataframe with columns: vertex, coordinates, ROI number, and ROI label
aparcrois <- list.files('~/Data/aparc', full.names=T)
aparcrois <- aparcrois[-grep('corpuscallosum|unknown', aparcrois)]
aparc <- do.call(rbind, lapply(1:length(aparcrois), function(x) atlas.table(aparcrois, x)))

write.table(aparc %>% filter(grepl('^rh\\.', name)), '~/LaxKAT/dk_right_atlas.txt', row.names=F, quote=F)

write.table(data.frame(label=1:68, 
                       name=(aparc %>% group_by(label) %>% dplyr::slice(1) %>% ungroup())$name), 
            '~/LaxKAT/dk_dictionary.txt', 
            row.names=F,
            quote=F)