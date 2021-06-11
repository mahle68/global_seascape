# Scripts for step selection function analysis
# This is script 2 of x for reproducing the results of Nourani et al 2021, ProcB.
# Elham Nourani, PhD. Jun.10. 2021
#-----------------------------------------------------------------

#to do: add the wind support functions to functions.R

library(tidyverse)
library(lubridate)


setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/")

# ---------- STEP 1: load data #####
load("2021/ssf_input_ann_cmpl_60_15.RData") #ann_cmpl; This dataframe includes used and alternative steps and can be reproduced using step_generation.R