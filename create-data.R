#Create dataframes to use in the model, from qfish data 
# CJ Brown 09-04-2020

#Redspot king prawn extract 
# https://qfish.fisheries.qld.gov.au/query/b1e142ce-a3c8-4bdb-9cb6-39bd0c9a423e/table?customise=True#
#Red-spot king prawns: 28711048, 28711910, 28711908


rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

cdat <- read.csv("data-raw/line-fishing-species-kg-GBR-by-fin-years.csv", 
                 na.strings = c("", "N/A"), skip = 1) 

edat <- read.csv("data-raw/line-fishing-days-GBR-by-fin-years.csv", 
                 na.strings = c("", "N/A"), skip = 1)

rspdat <- read.csv("data-raw/red-spot-3codes-GBR_tonnes.csv") %>%
  rename(days = Days, catch = Tonnes) %>%
  mutate(COMMON_NAME = "RSK prawn")
gsdat <- read.csv("data-raw/goldband-snapper.csv") %>%
  rename(days = Days, catch = Tonnes) %>%
  mutate(COMMON_NAME = "Goldband snapper")


caab <- read.csv("../caab_dump_latest.csv") %>%
  mutate(SPCODEX = paste0("X", SPCODE)) %>%
  select(SPCODEX, COMMON_NAME, SCIENTIFIC_NAME, FAMILY)

#
# Data wrangling 
#


cdat <- cdat %>% 
  pivot_longer(cols = -1, names_to = "SPCODEX", values_to = "catch") %>%
  left_join(caab)


# Identify most caught species 

spcounts <- cdat %>% 
  #merge goldband and Goldband snapper 
  mutate(SPCODEX = ifelse(SPCODEX == "X37346901", 
                              "X37346002", SPCODEX)) %>% 
  group_by(SPCODEX) %>%
  summarize(tcatch = sum(catch, na.rm = TRUE)) %>%
  arrange(desc(tcatch)) %>%
  left_join(caab)

spcounts[1:25,] %>% data.frame()

#Some of the top most caught with at least genus level ID
spkeep <- data.frame(COMMON_NAME = c(
  "coral trout", #plectropomus/variola, rank 1 by total catch 
  "Redthroat Emperor", #Lethrinus miniatus, rank 2
  "Saddletail Snapper", #Lutjanus malabaricus, rank 6
  "Red Emperor", #Lutjanus sabae, rank 7
  # "Goldband snapper", #Pristipomoides multidens 
  # + Pristipomoides multidens & Pristipomoides typus ('goldband snapper')
  # rank 8
  "Spangled Emperor", #Lethrinus nebulosus, rank 13
  # "goldband snapper", 
  "Stripey snapper" #Lutjanus carponotatus, rank 24, but high dependence spp
))

#Join to cdat
cdat2 <- cdat %>% inner_join(spkeep) %>%
  filter(FinancialYear != "Grand Total") %>%
  mutate(year = as.numeric(substring(FinancialYear, 1,4)))

edat <- edat %>% 
  pivot_longer(cols = -1, names_to = "SPCODEX", values_to = "days")

dat <- inner_join(edat, cdat2) %>%
  mutate(catch = catch/1000) %>%
  select(-SPCODEX, -SCIENTIFIC_NAME, -FAMILY) %>%
  bind_rows(rspdat) %>%
  bind_rows(gsdat) %>%
  filter(year < 2019) %>%
  filter(year > 1989)


ggplot(dat) + 
  aes(x = year, y = catch, color = COMMON_NAME) +
  geom_line() + 
  facet_wrap(~COMMON_NAME, scales = "free")
  
  

save(dat, file = "data-raw/qfish-dat.rda")
