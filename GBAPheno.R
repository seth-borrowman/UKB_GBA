library(tidyverse)
library(lubridate)

setwd("Z:/UKB Research/GBA1")

cardio <- read.csv("GBACardio.csv")
covar <- read.csv("GBACovar.csv")
hemat <- read.csv("GBAHemat.csv")
hepat <- read.csv("GBAHepat.csv")
musc <- read.csv("GBAMusc.csv")
neuro <- read.csv("GBANeuro.csv")
ocular <- read.csv("GBAOcular.csv")
ICD9 <- read.csv("TBI_ICD9.csv")
ICD10 <- read.csv("TBI_ICD10.csv")

# Looks like no ICD9 instances of TBI
rm(ICD9)

### Cardio

names(cardio) <- c("IID", "J84", "I27", "I70", "I60", "I61", "I62", "I63",
                   "I64")
cardio <- cardio %>%
    mutate(AnyCardio = case_when(
        (J84 != "") |
            (I27 != "") |
            (I70 != "") |
            (I60 != "") |
            (I61 != "") |
            (I62 != "") |
            (I63 != "") ~ 1,
        .default = 0
    ))

## Hemat

names(hemat) <- c("IID", "D64",
                  "Hemoglobin_i0", "Hemoglobin_i1", "Hemoglobin_i2",
                  "D69",
                  "Platelets_i0", "Platelets_i1", "Platelets_i2")
hemat$Hemoglobin <- rowMeans(hemat[,c(3:5)], na.rm = T)
hemat$Hemoglobin[which(is.nan(hemat$Hemoglobin))] <- NA
hemat$Platelets <- rowMeans(hemat[,c(7:9)], na.rm = T)
hemat$Platelets[which(is.nan(hemat$Platelets))] <- NA

sex <- covar %>% select(eid, p31)
names(sex) <- c("IID", "p31")
hemat <- left_join(hemat, sex, by = "IID")

hemat <- hemat %>%
    mutate(AnyHemat = case_when(
        (D64 != "") |
            (D69 != "") |
            (Hemoglobin < 11.9 & p31 == "Female") |
            (Hemoglobin < 13.6 & p31 == "Male") ~ 1,
        .default = 0
    ))

### Hepat

names(hepat) <- c("IID", "K74",
                  "LiverVol_i2", "LiverVol_i3",
                  "SpleenVol_i2", "SpleenVol_i3",
                  "E55",
                  "VitD_i0", "VitD_i1")
hepat$LiverVol <- rowMeans(hepat[,c(3:4)], na.rm = T)
hepat$SpleenVol <- rowMeans(hepat[,c(5:6)], na.rm = T)
hepat$vitD <- rowMeans(hepat[,8:9], na.rm = T) # nmol/L

height <- covar %>% select(eid, p50_i0, p50_i1, p50_i2, p50_i3)
height$height <- rowMeans(height[,2:5], na.rm = T)
height <- height %>% select(eid, height)
names(height) <- c("IID", "height")
hepat <- left_join(hepat, height, by = "IID")

hepat <- hepat %>%
    mutate(LiverVol = LiverVol / height,
           SpleenVol = SpleenVol / height)
hepat <- hepat %>%
    mutate(LiverZ = (LiverVol - mean(LiverVol, na.rm = T)) / sd(LiverVol, na.rm = T),
           SpleenZ = (SpleenVol - mean(SpleenVol, na.rm = T)) / sd(SpleenVol, na.rm = T))
# Really not enough data to justify using this

hepat <- hepat %>%
    mutate(AnyHepat = case_when(
        (K74 != "") |
            (E55 != "") |
            (vitD < 30) ~ 1,
        .default = 0
    ))

### Musc

names(musc) <- c("IID", "M85", "M87", "M80", "M81", "P94", "M40", "M41",
                 "10yoHeight_i0", "10yoHeight_i1", "10yoHeight_i2")
musc <- musc %>%
    mutate(AnyMusc = case_when(
        (M85 != "") |
            (M87 != "") |
            (M80 != "") |
            (M81 != "") |
            (P94 != "") |
            (M40 != "") |
            (M41 != "") ~ 1,
        .default = 0
    ))

### Neuro

names(neuro) <- c("IID", "G11", "G20", "G21", "G22", "G23", "G24", "G25", "G40",
                  "F00", "F01", "F03", "G91", "G56", "G57", "G58", "G59", "G60",
                  "G61", "G62", "G63", "H90", "H91")
neuro <- neuro %>%
    mutate(AnyNeuro = case_when(
        (G11 != "") |
            (G20 != "") |
            (G21 != "") |
            (G22 != "") |
            (G23 != "") |
            (G24 != "") |
            (G25 != "") |
            (G40 != "") |
            (F00 != "") |
            (F01 != "") |
            (F03 != "") |
            (G91 != "") |
            (G56 != "") |
            (G57 != "") |
            (G58 != "") |
            (G59 != "") |
            (G60 != "") |
            (G61 != "") |
            (G62 != "") |
            (G63 != "") |
            (H90 != "") |
            (H91 != "") ~ 1,
        .default = 0
    ))

neuro <- neuro %>%
    mutate(Parkinson = case_when(
        (G20 != "" & G21 == "" & G22 == "" & G23 == "") ~ 1,
        .default = 0
    ))

### Ocular

names(ocular) <- c("IID", "H25", "H26", "H27", "H28", "H34", "H35", "H49",
                   "H50", "H55")
ocular <- ocular %>%
    mutate(AnyOcular = case_when(
        (H25 != "") |
            (H26 != "") |
            (H27 != "") |
            (H28 != "") |
            (H34 != "") |
            (H35 != "") |
            (H49 != "") |
            (H50 != "") |
            (H55 != "") ~ 1,
        .default = 0
    ))


### Covariates

names(covar) <- c("IID", "YOB", "DODeath_i0", "DODeath_i1", "Townsend",
                  "Sex", "BMI_i0", "BMI_i1", "BMI_i2", "BMI_i3",
                  "Alcohol_i0", "Alcohol_i1", "Alcohol_i2", "Alcohol_i3",
                  "EvrSmoke_i0", "EvrSmoke_i1", "EvrSmoke_i2", "EvrSmoke_i3",
                  "PackYr_i0", "PackYr_i1", "PackYr_i2", "PackYr_i3",
                  "Height_i0", "Height_i1", "Height_i2", "Height_i3")
covar$DODeath <- covar$DODeath_i0 # None in i1 not in i0
covar$BMI <- rowMeans(covar[,7:10], na.rm = T)
covar <- covar %>%
    mutate(Alcohol = case_when(
        is.na(Alcohol_i0) == F & Alcohol_i0 != "Prefer not to answer" ~
            Alcohol_i0,
        is.na(Alcohol_i1) == F & Alcohol_i1 != "Prefer not to answer" ~
            Alcohol_i1,
        is.na(Alcohol_i2) == F & Alcohol_i2 != "Prefer not to answer" ~
            Alcohol_i2,
        is.na(Alcohol_i3) == F & Alcohol_i3 != "Prefer not to answer" ~
            Alcohol_i3,
        Alcohol_i0 == "Prefer not to answer" | Alcohol_i1 ==
            "Prefer not to answer" | Alcohol_i2 == "Prefer not to answer" |
            Alcohol_i3 == "Prefer not to answer" ~ "Prefer not to answer",
        .default = NA
    ))
covar$Alcohol[which(covar$Alcohol == "")] <- "Prefer not to answer"
covar <- covar %>%
    mutate(EverSmoke = case_when(
        is.na(EvrSmoke_i0) == F ~ EvrSmoke_i0,
        is.na(EvrSmoke_i1) == F ~ EvrSmoke_i1,
        is.na(EvrSmoke_i2) == F ~ EvrSmoke_i2,
        is.na(EvrSmoke_i3) == F ~ EvrSmoke_i3,
        .default = NA
    ))
covar$PackYears <- rowMeans(covar[,19:22], na.rm = T)
covar$PackYears[which(covar$EverSmoke == "Yes" &
                          is.nan(covar$PackYears))] <- median(covar$PackYears,
                                                              na.rm = T)

### TBI

covar <- covar %>%
    mutate(TBI = case_when(
        IID %in% ICD10$eid ~ 1,
        .default = 0
    ))


### Make final dataframe

covar <- covar %>%
    select(IID, YOB, Townsend, Sex, DODeath, BMI, Alcohol, EverSmoke, PackYears,
           TBI)

### Output

cardio <- cardio %>% select(IID, AnyCardio)
hemat <- hemat %>% select(IID, AnyHemat)
hepat <- hepat %>% select(IID, AnyHepat)
musc <- musc %>% select(IID, AnyMusc)
neuro <- neuro %>% select(IID, AnyNeuro, Parkinson)
ocular <- ocular %>% select(IID, AnyOcular)

dflist <- list(covar, cardio, hemat, hepat, musc, neuro, ocular)
out <- dflist %>% reduce(full_join, by = "IID")
out <- out[order(out$IID), ]

write_csv(out, "Phenotypes.csv")
