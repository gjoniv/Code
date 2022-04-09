# Log-transform SE as well as mean? 
# Nitrate (main effect) DOES NOT matter if data with truncated nitrate values are omitted. Vice versa if those data are included.

# # # APENDIX OF THE PAPER # # # 
# Colony-based CCSR 

# Set working directory
setwd('~/Dropbox/Work stuff/Projects/Savina_size_scaling_project/3. analysis/')
setwd('C:/Users/Thomasm/Dropbox/Work stuff/Projects/Savina_size_scaling_project/3. analysis')
setwd('~/Dropbox/Savina_size_scaling_project/3. analysis/')

library(modelsummary)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(dplyr)
library(ggeffects)


# Dataset
# Replace Nitrate 0 Values with half of apparent detection limit
# Replace temperature values above 50 with NAs
dat <- read.csv("ccsr_colonybased_summaries_2020_11_09.csv") %>%
  select(-c(site_id, uid, zoo_se_biomass)) %>%
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

# Histogram of the Nitrate_mean
# hist(log10(dat$nitrate_mean))
# min(dat$nitrate_mean, na.rm = TRUE)
# min(dat$nitrate_mean, na.rm = TRUE, dat$nitrate_mean > 0)

# FILTERING OUTLIER LAKE
dat <- dat %>%
  filter(phyto_mean_biovol > 0) 

# FITTING PREDICTED INTERACTIONS
mod_final <- lm(phyto_total_density ~ phyto_mean_biovol + 
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
                  phyto_mean_biovol:nitrate_mean:temp_mean, 
                  #  zoo_total_biomass:phyto_mean_biovol:temp_mean + 
                  #  zoo_total_biomass:nitrate_mean:phyto_mean_biovol +
                  # phyto_mean_biovol:zoo_total_biomass:nitrate_mean:temp_mean, 
                data = dat)

summary(mod_final)


bbmle::AICtab(mod_2way, mod_final, mod_all)
# mod_final is still best

sjPlot::tab_model(mod_final)

# library(texPreview)
# library(equatiomatic)
equatiomatic::extract_eq(mod_final)


#  #APPENDIX A2 
# Plotting Two-Way Interaction Effects of Regression Models
interact_plot(mod_final, pred = phyto_mean_biovol, 
              modx = temp_mean, 
              mod2= nitrate_mean, 
              plot.points = TRUE, 
              point.size = 2,
              interval = FALSE,
              modx.values = c(10, 20, 30), 
              mod2.values = c(-3, 0), 
              colors = c("darkblue", "seagreen", "yellow"), 
              x.label = "Log Cell Volume", 
              y.label = "Log Density", 
              legend.main = "T (°C)") + theme_few() + 
  scale_color_viridis(option = "C")

#  #APPENDIX A3
# Plotting Two-Way Interaction Effects of Regression Models
interact_plot(mod_final, pred = phyto_mean_biovol, 
              mod2 = temp_mean, 
              modx= nitrate_mean, 
              mod2.values = c(10,  30), 
              modx.values = c(-3, -1.5, 0),
              colors = c("orange", "purple", "darkblue"), 
              plot.points = TRUE, 
              point.size = 2, 
              interval = FALSE,
              x.label = "Log Cell Volume", 
              y.label = "Log Density", 
              legend.main = "N (mg/L)") + theme_few() + 
  scale_color_viridis(option = "D", direction = -1)

