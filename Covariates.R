library(tidyverse)

setwd("Z:/UKB Research/GBA1")
load("CleanVariantsWrkspc.RData")
pheno <- read_csv("Phenotypes.csv")
pca <- read_csv("GBA_PCA.csv")
ethnicity <- read_csv("Ethnicity.csv")

### Factor columns ----

IIDs <- which(colnames(plinkPath) %in% c("IID", pathVars$Variant))
plinkPath <- plinkPath[ ,IIDs]
pheno$Sex <- factor(pheno$Sex, labels = c("Female", "Male"))
pheno$Sex <- relevel(pheno$Sex, ref = "Female")
pheno$BMI_cat <- cut(pheno$BMI,
                 c(0, 18.5, 25, 30, Inf),
                 c("Underweight", "Normal", "Overweight", "Obese"),
                 include.lowest = T)
pheno$BMI_cat <- relevel(pheno$BMI_cat, ref = "Normal")
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
pheno$Parkinson <- factor(pheno$Parkinson, levels = c(0, 1))
pheno$TBI <- factor(pheno$TBI, levels = c(0, 1))
# Create ethnicity variables
ethnicity <- ethnicity %>%
    mutate(Ethnicity = case_when(
        all(is.na(p21000_i0), is.na(p21000_i1), is.na(p21000_i2),
            is.na(p21000_i3)) ~ NA,
        all(is.na(p21000_i0), is.na(p21000_i1), is.na(p21000_i2)) ~ p21000_i3,
        all(is.na(p21000_i0), is.na(p21000_i1)) ~ p21000_i2,
        is.na(p21000_i0) ~ p21000_i1,
        .default = p21000_i0
    )) %>%
    mutate(ethnic_group = case_when(
        Ethnicity %in% c("White", "British", "Irish",
                         "Any other white background") ~ "White",
        Ethnicity %in% c("Mixed", "White and Black Caribbean",
                         "White and Asian", "White and Black African",
                         "Any other mixed background") ~ "Mixed",
        Ethnicity %in% c("Asian or Asian British", "Indian", "Pakistani",
                         "Bangladeshi", "Chinese",
                         "Any other Asian background") ~ "Asian",
        Ethnicity %in% c("Black or Black British", "Caribbean", "African",
                         "Any other Black background") ~ "Black",
        Ethnicity == "Other ethnic group" ~ "Other",
        Ethnicity == "Do not know" ~ "Do not know",
        Ethnicity == "Prefer not to answer" ~ "Prefer not to answer",
        is.na(Ethnicity) ~ NA
    )) %>%
    select(eid, ethnic_group) %>%
    rename(IID = eid)
# Add ethnicity to phenotype data
pheno <- merge(pheno, ethnicity)
# Add 10 PCs
pca <- rename(pca, IID = eid)
pheno <- merge(pheno, pca)
# Select columns
pheno1 <- pheno
pheno <- pheno %>% dplyr::select(-c(EverSmoke, DODeath, AnyCardio, AnyHemat,
                                    AnyHepat, AnyMusc, AnyNeuro, AnyOcular,
                                    ethnic_group, BMI_cat, Parkinson))

### Clean up ----
# Remove NA's and unneeded levels
pheno <- na.omit(pheno)
pheno <- droplevels(pheno)
plinkPath <- na.omit(plinkPath)

# Filter to matching IID
pheno <- pheno %>% filter(IID %in% plinkPath$IID)
plinkPath <- plinkPath %>% filter(IID %in% pheno$IID)
pheno <- pheno[order(pheno$IID),]
plinkPath <- plinkPath[order(plinkPath$IID),]
plinkPath <- droplevels(plinkPath)

# Remove columns with <2 levels
removecols <- c()
for (i in 2:ncol(plinkPath)) {
    if (nlevels(plinkPath[,i]) < 2) {
        removecols <- append(removecols, as.integer(i))
    }
}
plinkPath <- plinkPath[,-removecols]
pathVars <- pathVars[which(pathVars$Variant %in% colnames(plinkPath)),]

# Export for PheWAS ----
for (i in 2:ncol(plinkPath)) {
    new <- plinkPath %>% select(IID, names(plinkPath)[i])
    name <- names(new)[2] %>% gsub(":", "_", .) %>% gsub(">", "_", .)
    write_csv(new, sprintf("PheWAS/plink_%s.csv", name), append = F)
}
write_csv(pheno, "PheWAS/covariates.csv", append = F)

save.image("AfterAnalysis.RData")

### Create table 1 ----
pheno1 <- pheno1 %>%
    mutate(AnyVar = case_when(
        IID %in% AnyVars$IID ~ 1,
        .default = 0
    ))
phecode <- read_csv("PheWAS/phecode.csv")
excl_codes <- phecode_exclude$exclusion_criteria[which(phecode_exclude$code == 332)]
excl <- which(colnames(phecode) %in% excl_codes)
                   
excl_phecode <- phecode[,c(1, excl)] %>% filter_all(any_vars(. %in% c(T)))
phecode <- phecode %>%
    mutate(Parkinson = case_when(
        `332` == T  ~ 1,
        .default = 0
    )) %>%
    mutate(Parkinson = case_when(
        id %in% excl_phecode$id ~ 0,
        .default = Parkinson
    ))
pheno1 <- pheno1 %>%
    select(-Parkinson) %>%
    rename(id = IID) %>%
    merge(., select(phecode, c(id, Parkinson)))
pheno1$Parkinson <- factor(pheno1$Parkinson)

pheno1 <- pheno1[which(pheno1$id %in% pheno$IID),]
table1::table1(~ YOB + Townsend + Sex + BMI + Alcohol + EverSmoke + PackYears +
                   TBI + Parkinson + ethnic_group | AnyVar, data = pheno1)

### Plot PCA ----
# Plot those having any variant on PCA
plot <- ggplot(pheno[which(pheno$AnyVar != 1),], aes(x = p22009_a1,
                                                      y = p22009_a2)) +
    geom_point(alpha = 0.1, shape = 16, size = 1) +
    geom_point(data = pheno[which(pheno$AnyVar == 1),],
               aes(x = p22009_a1, y = p22009_a2),
               color = "red", alpha = 0.3, shape = 16, size = 1) +
    theme_minimal() +
    xlab("PC1") +
    ylab("PC2") +
    ggtitle("Principal components of those having any GBA variant")
plot


df <- data.frame(noVar = c(438621, 3522),
                 Var = c(10249 , 141))
test <- chisq.test(df)
test
t.test(pheno1$PackYears[which(pheno1$AnyVar == 0)], pheno1$PackYears[which(pheno1$AnyVar == 1)])
wilcox.test(pheno1$PackYears[which(pheno1$AnyVar == 0)], pheno1$PackYears[which(pheno1$AnyVar == 1)])
