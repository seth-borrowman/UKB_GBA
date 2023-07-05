library(tidyverse)

setwd("Z:/UKB Research/GBA1")
load("CleanVariantsWrkspc.RData")
pheno <- read_csv("Phenotypes.csv")

# Factor columns
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
pheno <- pheno %>% select(-c(EverSmoke, DODeath))
pheno$Parkinson <- factor(pheno$Parkinson, levels = c(0, 1))
pheno$TBI <- factor(pheno$TBI, levels = c(0, 1))

# Remove NA's and unneeded levels
pheno <- na.omit(pheno)
pheno <- droplevels(pheno)

# Filter to matching IID
pheno <- pheno %>% filter(IID %in% plinkPath$IID)
plinkPath <- na.omit(plinkPath)
plinkPath <- plinkPath %>% filter(IID %in% pheno$IID)
pheno <- pheno[order(pheno$IID),]
plinkPath <- plinkPath[order(plinkPath$IID),]
plinkPath <- droplevels(plinkPath)

# Remove columns with <2 levels
removecols <- c()
for (i in 2:ncol(plinkPath)){
    if (nlevels(plinkPath[,i]) < 2){
        removecols <- append(removecols, as.integer(i))
    }
}
plinkPath <- plinkPath[,-removecols]

### Parkinson's

park <- merge(pheno, plinkPath, by = "IID") %>%
    select(-c(AnyCardio, AnyHemat, AnyHepat, AnyMusc, AnyNeuro, AnyOcular, IID))

y_park <- as.matrix(as.numeric(park$Parkinson)) - 1
x_park <- park %>% select(-Parkinson)
x_park <- model.matrix(~ . -1, data = x_park)
x_park <- x_park[,-3] # Keeps female in there for some reason - need to remove

mod_park <- glm(y_park ~ x_park, family = "binomial")

### Any Gaucher related outcome

y_any <- pheno %>% 
    select(c(AnyCardio, AnyHemat, AnyHepat, AnyMusc, AnyNeuro, AnyOcular))
y_any <- as.matrix(y_any)
x_any <- merge(pheno, plinkPath, by = "IID") %>%
    select(-c(AnyCardio, AnyHemat, AnyHepat, AnyMusc, AnyNeuro, AnyOcular, IID,
              Parkinson))
x_any <- model.matrix(~ . -1, data = x_any)
x_any <- x_any[,-3] # Keeps female in there for some reason - need to remove

mod_any <- joinet::joinet(Y = y_any, X = x_any, family = "binomial", trace.it = T)
