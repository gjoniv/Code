# # Appendix A1 # # 
# Cell-based CCSR 

# Set working directory
setwd('~/Dropbox/Savina_size_scaling_project/3. analysis/')

library(modelsummary)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(dplyr)
library(ggeffects)
library(interactions)
library(ggthemes)

# Dataset
# Replace Nitrate 0 Values with half of apparent detection limit
# Replace temperature values above 50 with NAs
dat <- read.csv("ccsr_cellbased_summaries_2020_10_28.csv") %>%
  select(-c(site_id, uid)) %>% 
  filter(phyto_mean_biovol > 0) %>%
  mutate(temp_mean = replace(temp_mean, temp_mean > 50, NA),
         nitrate_mean = log10(replace(nitrate_mean, nitrate_mean == 0, 0.001/2)),
         silicate_mean = log10(silicate_mean),
         zoo_total_biomass = log10(zoo_total_density * zoo_mean_biomass),
         zoo_mean_biomass = log10(zoo_mean_biomass),
         zoo_total_density = log10(zoo_total_density),
         phyto_mean_biovol = log10(phyto_mean_biovol),
         phyto_total_density = log10(phyto_total_density),
         phyto_se_biovol = log10(phyto_se_biovol))
dat2 <- data.frame(scale(dat[,c(1,3:10)], scale = FALSE), phyto_se_biovol = dat$phyto_se_biovol)
dat3 <- dat %>%
  mutate(nitrate_mean = replace(nitrate_mean, nitrate_mean == min(nitrate_mean, na.rm = TRUE), NA))
dat3 <- data.frame(scale(dat3[,c(1,3:10)], scale = FALSE), phyto_se_biovol = dat3$phyto_se_biovol)

# FILTERING OUTLIER LAKE
dat <- dat %>%
  filter(phyto_mean_biovol > 0) 

# Histogram of the Nitrate_mean
# hist(log10(dat$nitrate_mean))
# min(dat$nitrate_mean, na.rm = TRUE)
# min(dat$nitrate_mean, na.rm = TRUE, dat$nitrate_mean > 0)

# FITTING PREDICTED INTERACTIONS INCLUDING PREDATION

mod_full <- lm(phyto_total_density ~ phyto_mean_biovol + 
                 temp_mean + 
                 nitrate_mean + 
                 zoo_total_biomass + 
                 phyto_mean_biovol:temp_mean + 
                 phyto_mean_biovol:zoo_total_biomass + 
                 phyto_mean_biovol:nitrate_mean + 
                 temp_mean:zoo_total_biomass + 
                 temp_mean:nitrate_mean + 
                 zoo_total_biomass:nitrate_mean +
                zoo_total_biomass:nitrate_mean:temp_mean +
                 phyto_mean_biovol:nitrate_mean:temp_mean +
                 zoo_total_biomass:phyto_mean_biovol:temp_mean +
                 zoo_total_biomass:nitrate_mean:phyto_mean_biovol +
                 phyto_mean_biovol:zoo_total_biomass:nitrate_mean:temp_mean, 
               data = dat, 
               weights = (1/phyto_se_biovol))

summary(mod_full)


# # FIGURE Appendix A1 # # 
# Plotting Three-Way Interaction Effects of Regression Models

pr <- ggpredict(mod_full, type = "fixed", 
                terms = c("phyto_mean_biovol",
                          "zoo_total_biomass [1, 2, 3]",
                          "nitrate_mean [-3.2, -2.2, -1.2]", "temp_mean [10, 20, 30]"))

plot(pr, limits = c(0, 8), show.legend = FALSE, colors = "viridis", add.data = TRUE) 



