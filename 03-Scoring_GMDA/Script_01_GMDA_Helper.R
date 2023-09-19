# ############################################################################ #
# ############################################################################ #
#                                                                              #
# ------------------------- WP10 Starmaze data ------------------------------- #
# Script_01_GMDA_Helper                                                        #
# Author: Patrizia Maier                                                       #
#                                                                              #
# ############################################################################ #
# ############################################################################ #


# ------------------------------------------------------------------------------
# ::: LOAD PACKAGES ::: #
# ------------------------------------------------------------------------------

library(tidyverse)

# ############################################################################ #


# ::: Preprocessing for GMDA scoring: get goal objects for id ::: # 

# read-in data
my_path <- "../WP10_data/WP10_results/"
load(paste(my_path, "wp10_post_navigation_data.Rdata", sep=""))


# define function
giveMeGoals <- function(data){
  ID = readline(prompt = "Enter id: ")
  
  temp <- data %>% 
    filter(id==as.numeric(ID) & trial==4) %>% 
    select(id, obj_MA, obj_MC, obj_MI)
  
  print(temp)
}


# call function 
giveMeGoals(pt_data)


# ############################################################################ #

# clear workspace
rm(list = ls())