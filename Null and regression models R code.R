# Null and regression models in ecology

rm(list=ls())
getwd()
setwd("C:/Users/charl/Documents/University/Year 5 Masters/Semester 1/Biodiversity Under Pressure/8. Week Eight/Week Eight Computer Practical 1")
getwd()

library(vegan)
library(gam)
library(scam)
library(gdm)
library(car)
library(betapart)
library(tidyverse)
library(zetadiv)


# 1. NULL MODELS: PERMUTATION ALGORITHMS
load("islands_Soc_Haw_null_models_practical.RData")