library(tidyverse)

setwd("Z:/Chandra Lab/002 Chandra Lab Projects/UKB Research/GBA1/Biochem")

biomarkers <- read_csv("HearingAndBiochem.csv") %>%
    dplyr::select(-c(p2247_i0, p2247_i1, p2247_i2, p2247_i3)) %>%
    dplyr::select(contains(c("eid", "i0"))) %>%
    rename(IID = eid)

covariates <- read_csv("covariates.csv")
covariates$Sex <- factor(covariates$Sex, levels = c("Female", "Male"))
covariates$Alcohol <- factor(covariates$Alcohol,
                             labels = c("Never",
                                        "Prefer not to answer",
                                        "Special occasions only",
                                        "One to three times a month",
                                        "Once or twice a week",
                                        "Three to four times a week",
                                        "Daily or almost daily"))
covariates$TBI <- factor(covariates$TBI,
                         levels = c(0, 1),
                         labels = c("No", "Yes"))

rsid_files <- list.files(path = getwd(),
                         pattern = glob2rx('plink*csv'))

summary <- data.frame(Variant = c("delete"),
                      Estimate = c(0),
                      `Std. Error` = c(0),
                      `t value` = c(0),
                      `Pr(>|t|)` = c(0),
                      Marker = c("delete"))
for (i in 2:ncol(biomarkers)) {
    outcome <- biomarkers[,c(1,i)]
    marker <- colnames(outcome[2])
    colnames(outcome)[2] <- "Biomarker"
    for (file in rsid_files) {
        genotype <- read_csv(file)
        genotype[,2] <- lapply(genotype[,2], as.factor)
        independent <- merge(covariates, genotype)
        formod <- merge(outcome, independent) %>%
            dplyr::select(-IID)
        try(mod <- glm(Biomarker ~ ., data = formod))
        modsummary <- c(colnames(genotype)[2], coef(summary(mod))[24,], marker)
        summary <- rbind(summary, modsummary)
    }
}

summary <- summary[which(summary[,1] != "delete"),]
write.csv(summary, "summary.csv", row.names = F)
