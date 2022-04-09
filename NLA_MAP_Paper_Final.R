library(tidyverse)
library(ggthemes)
library(sp)
library(rgdal)
library(viridis)
library(dplyr)
library(purrr)
library(broom)
library(janitor)

setwd('~/Dropbox/Savina_size_scaling_project/1. data/NLA2012')
dat <- read.csv("nla2012_wide_siteinfo_08232016.csv")
chemdat <-  read.csv('nla2012_waterchem_wide.csv') %>%
  clean_names()
physdat <- read.csv('nla2012_wide_profile_08232016.csv') %>%
  clean_names()

chemdat <- select(chemdat, c(uid, nitrate_n_result)) %>%
  # filter(visit_no == 1) %>%
  group_by(uid) %>%
  summarise(nitrate_mean = mean(nitrate_n_result, na.rm = TRUE))

physdat <- physdat %>% 
  filter(visit_no == 1) %>%
  select(., c(uid, temperature)) %>%
  group_by(uid) %>%
  summarise(temp_mean = mean(temperature, na.rm = TRUE))

dat <- dat %>%
  rename(uid = UID)


dat2 <- dat %>%
  left_join(chemdat, by = "uid")

dat3 <- dat2 %>%
  left_join( physdat, by = "uid")

dat3 <- dat3 %>%
  mutate(temp_mean = replace(temp_mean, temp_mean > 50, NA)) %>%
  mutate(nitrate_mean = log10(replace(nitrate_mean, nitrate_mean == 0, 0.001/2))) %>%
  filter(., nitrate_mean != 'NA') %>%
  filter(., temp_mean != 'NA') %>%
  filter(nitrate_mean < 1) 

states <- map_data("state")
coordinates(states) <- c("long", "lat")
proj4string(states) <- CRS("+proj=longlat +ellps=clrk66")
summary(states)
states2 <- spTransform(states, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96"))
summary(states2)
states2 <- as.data.frame(states2)

ggplot() +
  geom_polygon(data=states2, aes(x=long, y=lat, group=group),
               color="black", show.legend = F, fill="#FFFFFF") +
  geom_point(data = dat3, mapping = aes(x=XCOORD, y=YCOORD, col=temp_mean), 
             inherit.aes = FALSE, alpha = 1, size = 1.5, alpha=0.5) +
  scale_color_viridis(option = "C") +
  theme_void()

ggplot() +
  geom_polygon(data=states2, aes(x=long, y=lat, group=group),
               color="black", show.legend = F, fill="#FFFFFF") +
  geom_point(data = dat3, mapping = aes(x=XCOORD, y=YCOORD, col=nitrate_mean), 
             inherit.aes = FALSE, alpha = 1, size = 1.5, alpha=0.5) +
  scale_color_viridis(option = "D", direction = -1) +
  theme_void()
