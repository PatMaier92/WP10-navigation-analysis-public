# ############################################################################ #
# ############################################################################ #
#                                                                              #
# ------------------------- WP10 Starmaze data ------------------------------- #
# Script_02_GMDA_to_Rdata                                                      #
# Author: Patrizia Maier                                                       #
#                                                                              #
# ############################################################################ #
# ############################################################################ #                                   


# ------------------------------------------------------------------------------
# ::: LOAD PACKAGES ::: #
# ------------------------------------------------------------------------------

library(tidyverse)

# ############################################################################ #  

# ------------------------------------------------------------------------------
# ::: PROCESS RESULTS FILES AFTER SCORING ::: #
# ------------------------------------------------------------------------------

#  read-in data
protocol <- readline("Enter protocol folder: ")
my_path <- paste("GMDA/Data/", protocol, "/", sep="")
pattern <- list.files(path=my_path, pattern="_Summary\\.csv$")
file_list <- paste(my_path, pattern, sep="")


# process gmda data 
gmda_data <- file_list %>%
  purrr::map_df(read_csv,
                col_names=c("Measure Type", "Filename", "Measure", "Score", "Score_2"), 
                col_type=cols(`Measure Type`=col_character(),
                              Filename=col_character(),
                              Measure=col_character(),
                              Score=col_double(),
                              Score_2=col_double()),
                skip=9, 
                n_max=8) %>% 
  # correct for delimiter error in raw data 
  unite(Score, Score_2, col="Score", sep=".", na.rm=T) %>% 
  # add id and type info from Filename
  separate(Filename, sep="_", into=c("id", "Rotation"))
  

# process bdr data 
bdr_data <- file_list %>%
  purrr::map_df(read_csv, 
                col_names=c("Measure Type", "Filename", "Measure", "Score"), 
                col_type=cols(`Measure Type`=col_character(),
                              Filename=col_character(),
                              Measure=col_character(),
                              Score=col_double()),
                skip=21, 
                n_max=10) %>% 
  # add id and type info from Filename
  separate(Filename, sep="_", into=c("id", "Rotation"))


# # individual landmark gmda data
# gmda_ind_data <- file_list %>%
#   purrr::map_df(read_csv,
#                 col_names=c("Measure Type", "Filename", "Landmark", "Measure", "Score"),
#                 cols(`Measure Type`=col_character(),
#                      Filename=col_character(),
#                      Landmark=col_character(),
#                      Measure=col_character(),
#                      Score=col_double()),
#                 skip=35,
#                 n_max=35) %>% 
#   # add id and type info from Filename
#   separate(Filename, sep="_", into=c("id", "Rotation"))


# check for "bad" theta values (> 10)
threshold <- 10
# threshold <- as.numeric(readline("Enter theta threshold: "))
bad_theta <- bdr_data %>% filter(Measure=="theta") %>% filter(abs(Score) > threshold) 
bad_theta


# select and combine data 
temp_bdr <- bdr_data %>% 
  select(!c(`Measure Type`, Rotation)) %>%
  filter(Measure %in% c("r", "theta"))

data_gmda <- gmda_data %>% 
  select(!c(`Measure Type`, Rotation)) %>% 
  filter(!Measure %in% c("SQRT(Canonical Organization)", "Canonical Organization", "Rotational Bias", "Scaling Bias", "Num Landmarks Missing")) %>% 
  mutate(Score=as.numeric(Score),
         Measure=case_when(Measure=="Canonical Accuracy" ~ "CanAcc",
                           Measure=="Distance Accuracy" ~ "DistAcc",
                           Measure=="Angle Accuracy" ~ "AngleAcc",
                           TRUE ~ "NA")) %>% 
  full_join(temp_bdr, by=c("id", "Measure", "Score")) %>% 
  mutate(id=as.numeric(id),
         group=factor(case_when(id<12000 | (id>20000 & id<22000) ~ "YoungKids",
                                id<15000 | (id>20000 & id<25000) ~ "OldKids",
                                TRUE ~ "YoungAdults"), levels=c("YoungKids", "OldKids", "YoungAdults"))) %>%
  arrange(id) %>% 
  relocate("group", .after="id") %>% 
  rename(gmda_measure=Measure, score=Score)


# save as Rdata 
date <- readline("Enter date and scoring protocol info: ")

out_file_R <-  paste("GMDA/Data/wp10_GMDA_data_", date, ".Rdata", sep="")
save(data_gmda, file=out_file_R)


# ############################################################################ #  

# clear workspace
rm(list = ls())