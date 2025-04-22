library(tidyverse)
library(swfscMisc)

load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data/id.key.rda")

samps_w_SI <- read_csv("/Users/Shared/KKMDocuments/Documents/Karen/Structure/Pseudorca/Emma Carroll collaboration/all.samps.w.stable.isotopes.csv")
Pcra_DNA_concentration <- read_csv("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data-raw/Pcra_DNA_concentration.csv")
archive_data <- read_csv("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data-raw/archive-dump_07Aug2024.csv") |> 
  mutate(LABID = paste0('z0', zero.pad(LABID)))
social_cluster_assignments <- read_csv("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data-raw/social.cluster.assignments.csv")

Animals_w_SI <- left_join(samps_w_SI, 
                        filter(Pcra_id_key, id.type == "CRC.ID") |> 
                          rename(catalog_id = alt.id) |> 
                          select(Animal.ID, catalog_id)
                        ) |> 
  rename(CRC.ID = catalog_id) |> 
  left_join(social_cluster_assignments |> 
              rename(CRC.ID = ID, cluster = Cluster_Louvain)) |> 
#  select(Animal.ID) |> 
  distinct()

###    NEED TO ADD SOCIAL CLUSTER INFO TO THIS   #########
DQuant_samples_w_SI <- Pcra_DNA_concentration |> mutate(LabID = zero.pad(tis_LabID),
                                 LabID = paste0('z0', LabID)) |>
  left_join(
    Pcra_id_key |> 
      filter(id.type == 'LABID') |> 
      rename(LabID = alt.id) |> 
      select(Animal.ID, LabID)) |> 
#  filter(Animal.ID %in% Animals_w_SI$Animal.ID) |> 
  right_join(Animals_w_SI) |> 
  left_join(
    archive_data |> rename(LabID = LABID) |> 
      select(LabID, FIELDID)
  )
write.csv(DQuant_samples_w_SI, file = "DQuant_samples_w_SI.csv")

load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data/Pcra.sample.data.rda")
samp.dat <- samp.dat |> left_join(filter(Pcra_id_key, id.type == 'CRC.ID')) |> 
  select(-id.type) |> 
  rename(CRC.ID = alt.id)

lcWGS_samps <- left_join(
  right_join(samp.dat,
  read_csv("results/Pc_aligned_mean_depths.csv") |> 
  filter(Species == 'Pseudorca crassidens') |> 
    select(c('LABID', 'Mean coverage'))),
  read_csv("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data-raw/Pcra_Tissue_Vial_data.csv") |> 
    rename(LABID = dbo_vw_Tissue_Vial_LabID) |> 
    mutate(LABID = paste0('z0', zero.pad(LABID))),
  by = 'LABID'
) |> 
  left_join(CR.haps)
# |> 
#  filter(Tissue_Type == 'Skin & Blubber') #|> 
  #filter(TQuant >= 2)
write.csv(lcWGS_samps, file = 'for_Emma_lcWGS_samps.csv')
