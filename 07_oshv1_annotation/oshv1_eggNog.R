library(readxl)
library(tidyverse)


protein <- read_excel("Downloads/out.emapper.annotations_oshv1_proteins.xlsx", skip = 2) %>% 
  select(query, Description) %>% 
  filter(Description != "-") %>% 
  separate(query, into = c("query", "rest"), sep = "_") %>%
  select(query, Description)



dna <- read_excel("Downloads/out.emapper.annotations_oshv1_dna.xlsx", skip = 2) %>% 
  select(query, Description) %>% 
  filter(Description != "-") %>% 
  separate(query, into = c("query", "rest"), sep = "::") %>%
  select(query, Description)


oshv1_eggNOG <- left_join(dna, protein, by = "query", suffix = c("_dna_eggNOG", "_protein_eggNOG"))

write_tsv(oshv1_eggNOG, "Downloads/oshv1_eggNOG_combined")