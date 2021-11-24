# Set working directory
setwd('~/Dropbox/Work stuff/Projects/Savina_size_scaling_project/1. Data/NLA2012')
setwd('C:/Users/gjoni/Dropbox/Savina_size_scaling_project/1. Data/NLA2012')
setwd('C:/Users/Thomasm/Dropbox/Work stuff/Projects/Savina_size_scaling_project/1. Data/NLA2012')

# Libraries 
library(dplyr)
# library(purrr)
# library(broom)
library(janitor)
se <- function(x) sd(x)/sqrt(length(x))

# Read in data files
chemdat <-  read.csv('nla2012_waterchem_wide.csv')
physdat <- read.csv('nla2012_wide_profile_08232016.csv') %>%
  clean_names()
phytodat <- read.csv('nla2012_wide_phytoplankton_count_2020_04_25.csv') %>%
  clean_names()
zoodat_meta <- read.csv('nla2012_zooptaxa_wide_10272015.csv') %>%
  select(TARGET_TAXON, FFG) %>%
  rename(feeding_group = FFG)
zoodat_orig <- read.csv('nla2012_wide_zooplankton_count_2020_04_30.csv') %>%
  left_join(., zoodat_meta, by = 'TARGET_TAXON') %>%
  clean_names()

### I. Chemical data processing ####
# Selecting chemical data to merge using UID as identified 

chemdat <- select(chemdat, c(UID, PTL_RESULT, NITRATE_N_RESULT, SILICA_RESULT)) %>%
  # filter(visit_no == 1) %>%
  group_by(UID) %>%
  summarise(phos_total_mean = mean(PTL_RESULT, na.rm = TRUE),
            nitrate_mean = mean(NITRATE_N_RESULT, na.rm = TRUE),
            silicate_mean = mean(SILICA_RESULT, na.rm = TRUE)) %>%
  rename(uid = UID)

### II. Physical data processing
physdat <- physdat %>% 
  filter(visit_no == 1) %>%
  select(., c(uid, temperature)) %>%
  group_by(uid) %>%
  summarise(temp_mean = mean(temperature, na.rm = TRUE))

### III. Phytoplankton data processing  (using UID as identifier as the chem data doesnt have site_ID variable)
# filter data only to Visit no 1 and biovol > 0

# Phytoplankton mean biovolume and total density per lake
# replace NAs in unit_cells by 1
# use filter() to remove unit_cells values below 1 (for future: decide whether to keep them based on Richard's response)
# use mutate() to define a new variable natural_unit_density = density / unit_cells
# modify the calculation of phyto_mean_biovol to use natural_unit_density instead of density
phytodat <- phytodat %>% 
  mutate(unit_cells = tidyr::replace_na(unit_cells, 1)) %>%
  filter(., visit_no == 1) %>%
  filter(., density != 'NA') %>%
  filter(., biovolume != 'NA') %>%
  filter(., biovolume > 0) %>% 
  # filter(., unit_cells > 0.99) %>% 
  mutate(natural_unit_density = density / unit_cells,
         biovol_per_natural_unit = biovolume / natural_unit_density)

phytodatsummary <- phytodat %>%
  group_by(uid) %>%
  summarise(phyto_mean_biovol = sum(biovolume) / sum(natural_unit_density),
            # phyto_se_biovol = se(biovolume),
            phyto_se_biovol = sd(rep(biovolume, natural_unit_density)) / sqrt(sum(natural_unit_density)),
            phyto_total_density = sum(natural_unit_density))

### IV. Zooplankton data processing

# Biomass and density per species
# Grouping by target_taxon because some species are in multiple size classes
zoodatsummary <- zoodat_orig %>%
  filter(., visit_no == 1) %>%
  filter(., sample_type == 'ZOFN' | sample_type == 'ZOCN') %>%
  filter(., feeding_group != 'PARA') %>%
  filter(., feeding_group != 'PRED') %>%
  # filter(., feeding_group == 'HERB') %>%
  mutate(biomass = as.numeric(as.character(na_if(biomass, '#DIV/0!')))) %>%
  mutate(density = as.numeric(as.character(na_if(density, '#DIV/0!')))) %>%
  filter(., density != 'NA') %>%
  filter(., biomass != 'NA') %>%
  filter(., biomass > 0) %>%
  group_by(target_taxon, uid) %>%
  summarise(biomass = sum(biomass),
            density = sum(density),
            site_id = site_id[1]) %>%
  ungroup() %>%
  group_by(uid) %>%
  summarise(zoo_mean_biomass = sum(biomass) / sum(density),
            zoo_se_biomass = sd(rep(biomass, density)) / sqrt(sum(density)),
            zoo_total_density = sum(density))


### V. Merging CCSR datasets

ccsr_mergedat <- left_join(phytodatsummary, zoodatsummary, by = "uid") %>%
  left_join(., physdat, by ="uid") %>%
  left_join(., chemdat, by = "uid")

# Add Site_ID to the table, so both Site_ID and UID are present together
siteid <- unique(phytodat[,c('uid','site_id')])
ccsr_mergedat <- left_join(ccsr_mergedat, siteid, by = "uid")

### Save it as csv 
write.csv(ccsr_mergedat, paste0('C:/Users/Thomasm/Dropbox/Work stuff/Projects/Savina_size_scaling_project/3. analysis/ccsr_colonybased_summaries_',                    
                                gsub('-', '_', Sys.Date()), '.csv'), row.names = FALSE)

write.csv(ccsr_mergedat, paste0('~/Dropbox/Work stuff/Projects/Savina_size_scaling_project/3. analysis/ccsr_colonybased_summaries_',                    
                                gsub('-', '_', Sys.Date()), '.csv'), row.names = FALSE)

