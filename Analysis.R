library(tidyverse)

setwd("Z:/UKB Research/GBA1")
load("CleanVariantsWrkspc.RData")
pheno <- read_csv("Phenotypes.csv")

pheno$Sex <- factor(pheno$Sex, labels = c("Female", "Male"))
pheno$Sex <- relevel(pheno$Sex, ref = "Female")
pheno$BMI <- cut(pheno$BMI,
              c(0, 18.5, 25, 30, Inf),
              c("Underweight", "Normal", "Overweight", "Obese"),
              include.lowest = T)
pheno$BMI <- relevel(pheno$BMI, ref = "Normal")
pheno$Alcohol <- factor(pheno$Alcohol,
                        labels = c("Never",
                                   "Prefer not to answer",
                                   "Special occasions only",
                                   "One to three times a month",
                                   "Once or twice a week",
                                   "Three to four times a week",
                                   "Daily or almost daily"))
pheno$EverSmoke <- factor(pheno$EverSmoke, labels = c("No", "Yes"))
pheno$PackYears[which(is.na(pheno$PackYears))] <- 0
pheno <- pheno %>% select(-c(EverSmoke, AgeDeath))
pheno$Parkinson <- factor(pheno$Parkinson, levels = c(0, 1))
pheno <- na.omit(pheno)
pheno <- droplevels(pheno)

pheno <- pheno %>% filter(IID %in% plinkPath$IID)
plinkPath <- na.omit(plinkPath)
plinkPath <- plinkPath %>% filter(IID %in% pheno$IID)
pheno <- pheno[order(pheno$IID),]
plinkPath <- plinkPath[order(plinkPath$IID),]
plinkPath <- droplevels(plinkPath)
removecols <- c()
for (i in 2:ncol(plinkPath)){
    if (nlevels(plinkPath[,i]) < 2){
        removecols <- append(removecols, as.integer(i))
    }
}
plinkPath <- plinkPath[,-removecols]

park <- merge(pheno, plinkPath, by = "IID") %>%
    select(-c(AnyCardio, AnyHemat, AnyHepat, AnyMusc, AnyNeuro, AnyOcular, IID))

y <- as.matrix(as.numeric(park$Parkinson)) - 1
x <- park %>% select(-Parkinson)
x <- model.matrix(~ . -1, data = x)

mod <- glm(y ~ x, family = "binomial")

