library(tidyverse)
library(survival)
library(ggsurvfit)

setwd("Z:/UKB Research/GBA1")
load("CleanVariantsWrkspc.RData")
pheno <- read_csv("Phenotypes.csv")

# Factor columns
path_Vars <- pathVars[which(pathVars$ClinSig == "Pathogenic" |
                                   pathVars$ClinSig == "Pathogenic/Likely pathogenic" |
                                   pathVars$ClinSig == "Likely pathogenic" |
                                   pathVars$ClinSig == "Pathogenic/Likely pathognic; risk factor" |
                                    pathVars$REVEL >= 0.75),] # Filter further to REVEL > 0.75
plinkPath <- plinkPath[ ,which(colnames(plinkPath) %in% c("IID", path_Vars$Variant))]
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
pheno1 <- pheno
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
coef(summary(mod_park))[which(coef(summary(mod_park))[,4] < (0.05/56)),]

### Any Gaucher related outcome

y_any <- pheno %>% 
     select(c(AnyCardio, AnyHemat, AnyHepat, AnyMusc, AnyNeuro, AnyOcular))
# y_any <- as.matrix(y_any)
# x_any <- merge(pheno, plinkPath, by = "IID") %>%
#     select(-c(AnyCardio, AnyHemat, AnyHepat, AnyMusc, AnyNeuro, AnyOcular, IID,
#               Parkinson))
# x_any <- model.matrix(~ . -1, data = x_any)
# x_any <- x_any[,-3] # Keeps female in there for some reason - need to remove

# Are outcomes correlated?
cormatrix <- cor(y_any, method = "pearson")
round(cormatrix, 2)
caret::confusionMatrix(y_any)

# mod_any <- joinet::joinet(Y = y_any, X = x_any, family = "binomial", trace.it = T)
# coef(mod_any)
# weights(mod_any)


### only Clinvar
path_clinvar <- pathVars[which(pathVars$ClinSig == "Pathogenic" |
                        pathVars$ClinSig == "Pathogenic/Likely pathogenic" |
                        pathVars$ClinSig == "Likely pathogenic" |
                        pathVars$ClinSig == "Pathogenic/Likely pathognic; risk factor")]
park_clinvar <- park[,c(1:8, which(names(park) %in% path_clinvar$Variant))]

mod_park_clinvar <- glm(Parkinson ~ ., data = park_clinvar, family = "binomial")
coef(summary(mod_park_clinvar))[which(coef(summary(mod_park_clinvar))[,4] < (0.05/56)),]

### Survival
survivors <- pheno1 %>%
    select(IID, YOB, DODeath) %>%
    merge(., plinkPath, by = "IID") %>%
    select(IID, YOB, DODeath, rs1671872221, `1:155235769:G>A`, rs80356771)
survivors <- survivors %>%
    mutate(time = DODeath - ymd(YOB, truncated = 2L)) %>%
    mutate(status = case_when(
        is.na(DODeath) ~ 0,
        .default = 1
    ))
survivalAnalysis <- survdiff(Surv(time, status) ~ rs1671872221, data = survivors)
names(survivors)[5] <- "A155235769"
survivalAnalysis1 <- survdiff(Surv(time, status) ~ A155235769, data = survivors)
survivalAnalysis2 <- survdiff(Surv(time, status) ~ rs80356771, data = survivors)

survfit2(Surv(time, status) ~ rs1671872221, data = survivors) %>% 
    ggsurvfit() +
    labs(
        x = "Survival time (days)",
        y = "Overall survival probability"
    )

survfit2(Surv(time, status) ~ rs80356771, data = survivors) %>% 
    ggsurvfit() +
    labs(
        x = "Survival time (days)",
        y = "Overall survival probability"
    )
