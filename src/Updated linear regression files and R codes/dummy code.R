library(dplyr)

# test set 
demotestset <- read.csv("GWAS_75testset_all_afterBH_probe5_v2_converted.csv") # change probe no. accordingly

ncol(demotestset)

demotestset_dc <- demotestset %>% 
  mutate(
    mother_chinesevmalay = ifelse(mother_ethnicity_cat %in% c(1, 3), 0, as.numeric(mother_ethnicity_cat == 2)),
    mother_chinesevindian = ifelse(mother_ethnicity_cat %in% c(1, 2), 0, as.numeric(mother_ethnicity_cat == 3))
  )

ncol(demotestset_dc)

write.csv(demotestset_dc,"GWAS_75testset_all_afterBH_probe5_v2_converted_dc.csv", row.names = FALSE)
