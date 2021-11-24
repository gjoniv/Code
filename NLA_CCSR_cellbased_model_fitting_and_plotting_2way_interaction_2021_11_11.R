# # # MAIN RESULTS OF THE PAPER # # # 
# Cell-based CCSR 

# Set working directory
# setwd('~/Dropbox/Savina_size_scaling_project/3. analysis/')
setwd("C:/Users/Thomasm/Dropbox/Work stuff/Projects/Lab group/Savina_size_scaling_project/3. analysis")

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


## FILTERING OUTLIER LAKE
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
                  #zoo_total_biomass:nitrate_mean + 
                  # zoo_total_biomass:nitrate_mean:temp_mean + 
                  phyto_mean_biovol:nitrate_mean:temp_mean, 
                # zoo_total_biomass:phyto_mean_biovol:temp_mean + 
                # zoo_total_biomass:nitrate_mean:phyto_mean_biovol +
                # phyto_mean_biovol:zoo_total_biomass:nitrate_mean:temp_mean, 
                data = dat)

summary(mod_final)

# Just for reference: two models with just all 2-way interactions and with all possible interactions
mod_all <- lm(phyto_total_density ~ phyto_mean_biovol + 
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
              data = dat)


mod_2way <- lm(phyto_total_density ~ phyto_mean_biovol + 
                 temp_mean + 
                 nitrate_mean + 
                 zoo_total_biomass + 
                 phyto_mean_biovol:temp_mean + 
                 phyto_mean_biovol:zoo_total_biomass + 
                 phyto_mean_biovol:nitrate_mean + 
                 temp_mean:zoo_total_biomass + 
                 temp_mean:nitrate_mean, 
               #zoo_total_biomass:nitrate_mean + 
               # zoo_total_biomass:nitrate_mean:temp_mean + 
               #phyto_mean_biovol:nitrate_mean:temp_mean, 
               # zoo_total_biomass:phyto_mean_biovol:temp_mean + 
               # zoo_total_biomass:nitrate_mean:phyto_mean_biovol +
               # phyto_mean_biovol:zoo_total_biomass:nitrate_mean:temp_mean, 
               data = dat)

summary(mod_2way)


bbmle::AICtab(mod_2way, mod_final, mod_all)
# mod_final is still best

sjPlot::tab_model(mod_final)

# library(texPreview)
# library(equatiomatic)
equatiomatic::extract_eq(mod_final)


# mod_final <- lm(phyto_total_density ~ phyto_mean_biovol + temp_mean + nitrate_mean + 
#              zoo_total_biomass, data = dat, 
#              weights = (1/phyto_se_biovol))


pr <- ggpredict(mod_final, type = "fixed", 
                terms = c("phyto_mean_biovol",
                          "temp_mean [8, 20, 32]",
                          # "zoo_total_biomass [1, 2, 3]",
                          "nitrate_mean [-2.4, -1, 0]"))

plot(pr, limits = c(0, 10), show.legend = FALSE) 

# # # MAIN FIGURES OF THE PAPER

# # FIGURE 1
# Plotting CCSR Regression Model at cell level 
ggplot(data = dat, 
       aes(x = phyto_mean_biovol, 
           y = phyto_total_density)) + geom_point(color='black', 
                                                  alpha = 0.5) + geom_smooth(method = "lm", 
                                                                             color='black', 
                                                                             se = FALSE) + labs(x = "Log Cell Volume", 
                                                                                y = "Log Abundance") + theme_few()
# # FIGURE 2
# Plotting Two-Way Interaction Effects of Regression Models
interact_plot(mod_final, pred = phyto_mean_biovol, 
              modx = temp_mean, 
              mod2= nitrate_mean, 
              plot.points = TRUE, 
              point.size = 2,
              interval = TRUE,
              modx.values = c(10, 20, 30), 
              mod2.values = c(-3.2, -2.2, -1.2), 
              colors = c("blue", "green", "red"), 
              x.label = "Log Cell Volume", 
              y.label = "Log Abundance", 
              legend.main = "Temperature") + theme_few()

# # FIGURE 3
# Plotting Two-Way Interaction Effects of Regression Models 
interact_plot(mod_final, pred = phyto_mean_biovol, 
              mod2 = temp_mean, 
              modx= nitrate_mean, 
              mod2.values = c(10, 20, 30), 
              modx.values = c(-3.2, -2.2, -1.2), 
              plot.points = TRUE, 
              point.size = 2, 
              interval = TRUE,
              colors = c("darkgreen", "green", "darkorange"), 
              x.label = "Log Cell Volume", 
              y.label = "Log Abundance", 
              legend.main = "Nuntrients") + theme_few()
