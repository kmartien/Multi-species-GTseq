library(tidyverse)

samps <- read.csv(file = "qry_results/Steno.lcWGS.sample.picking.csv")

main.samples <- filter(samps, Rank == 1)

substitutes <- do.call(rbind, lapply(unique(main.samples$category), function(id){
  main.subs <- filter(samps, category == id) %>% 
    filter(Rank %in% c(2,3)) %>% arrange(Rank) %>%
    select(LABID)
  backups <- filter(samps, category == id) %>% 
    filter(is.na(Rank)) %>% select(LABID)
  return(data.frame(category = id, subs = paste(main.subs$LABID, collapse = ", "), alts = paste(backups$LABID, collapse = ", ")))
}))

main.samples <- left_join(main.samples, substitutes)

write.csv(main.samples, file = "data-raw/Steno.lcWGS.final.sample.list.csv")
