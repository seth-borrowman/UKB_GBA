library(tidyverse)
library(latex2exp)
library(ggrepel)

setwd("Z:/UKB Research/GBA1")
load("CleanVariantsWrkspc.RData")
pheno <- read_csv("Phenotypes.csv")
pca <- read_csv("GBA_PCA.csv")

# Factor columns
path_Vars <- pathVars[
    which(pathVars$ClinSig == "Pathogenic" |
          pathVars$ClinSig == "Pathogenic/Likely pathogenic" |
          pathVars$ClinSig == "Likely pathogenic" |
          pathVars$ClinSig == "Pathogenic/Likely pathognic; risk factor" |
          pathVars$REVEL >= 0.6),] # More strict than CleanVariants.R
plinkPath <- plinkPath[ ,which(colnames(plinkPath) %in% c("IID",
                                                          path_Vars$Variant))]
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
pheno <- pheno %>% dplyr::select(-c(EverSmoke, DODeath))
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
try(plinkPath <- plinkPath[,-removecols], silent = T)


### Parkinson's
colnames(pca)[1] <- "IID"
park <- merge(pheno, plinkPath, by = "IID") %>%
    merge(., pca) %>%
    dplyr::select(-c(AnyCardio, AnyHemat, AnyHepat, AnyMusc,
                     AnyNeuro, AnyOcular, IID)) %>%
    na.omit()

y_park <- as.matrix(as.numeric(park$Parkinson)) - 1
x_park <- park %>% dplyr::select(-Parkinson)
x_park <- model.matrix(~ . -1, data = x_park)
# Keeps female in there for some reason - need to remove
x_park <- x_park[,-which(colnames(x_park) == "SexFemale")]

# Association analysis with logistic regression for each SNV
summary <- matrix(rep(0, 5*(ncol(x_park)-24)), # -14 to skip covariates
                  ncol = 5, nrow = ncol(x_park)-24) %>%
    as.data.frame()
colnames(summary) <- c("Variant", "Estimate", "Std. Error", "z value",
                       "p")
for (i in 15:(ncol(x_park)-10)) {
    variants <- seq(15, ncol(x_park)-10)
    variants <- variants[which(variants != i)]
    x_new <- x_park[,-variants]
    mod <- glm(y_park ~ x_new, family = "binomial")
    summary[(i - 14), 1] <- colnames(x_new)[15]
    summary[(i - 14), 2:5] <- unname(coef(summary(mod))[16,])
}
summary[,2:5] <- sapply(summary[,2:5], as.numeric)
summary <- summary %>%
    # Bonferroni correction based on alpha 0.05
    mutate(Bonferroni = if_else(p <= (0.05/(nrow(summary)+1)), # +1 for any var
                                TRUE, FALSE)) %>%
    # Benjamini-Hochberg FDR with alpha 0.05
    mutate(FDR = if_else(p <= ((match(p,
                        summary$p[order(summary$p, decreasing = F)]) / 
                            (nrow(summary)+1))*0.05),
                        TRUE, FALSE)) %>%
    # Create label for plotting
    mutate(Label = gsub("`", "", substr(Variant, 1, nchar(Variant) - 1)))

# Find and set FDR based on lowest passing p-value
max_FDR_p <- max(summary[which(summary$FDR == T), 5], na.rm = T)
summary <- summary %>%
    mutate(FDR = if_else(summary$p <= max_FDR_p, T, F))
summary <- merge(summary, pathVars[,c(1, 3)], by.x = "Label", by.y = "Variant")


### Having any GBA variant
any_x <- pheno %>% dplyr::select(-c(AnyCardio, AnyHemat, AnyHepat, AnyMusc,
                                    AnyNeuro, AnyOcular)) %>%
    mutate(AnyVar = case_when(
        IID %in% AnyVars$IID ~ 1,
        .default = 0
    )) %>%
    merge(., pca)
any_x <- dplyr::select(any_x, -IID)
mod_any <- glm(Parkinson ~ ., data = any_x, family = "binomial")

# Interaction b/t TBI and having variant?
tbi_any <- glm(Parkinson ~ YOB + Townsend + Sex + BMI + Alcohol +
                   PackYears + p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 +
                   p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 +
                   p22009_a10 +
                   TBI*AnyVar,
                   data = any_x, family = "binomial") # no

# Add for plotting
mod_any_sum <- c(rep("Any Variant", 2), unname(coef(summary(mod_any))[16,]),
                      T, T, NA)
summary[nrow(summary)+1,] <- mod_any_sum
summary[,c(3:6, 9)] <- summary[,c(3:6, 9)] %>% apply(., 2, as.numeric)

### Plot results
plot <- ggplot(na.omit(summary), aes(x = Pos, y = -log10(p))) + 
    geom_point() +
    geom_hline(aes(yintercept = -log10(0.05/nrow(summary)),
               color = "Bonferroni", linetype = "Bonferroni")) +
    geom_hline(aes(yintercept = -log10(max_FDR_p), color = "FDR",
                   linetype = "FDR")) +
    scale_color_manual(values = c("Bonferroni" = "red", "FDR" = "blue"),
                       name = "Legend") +
    scale_linetype_manual(values = c("Bonferroni" = "solid", "FDR" = "dashed"),
                          name = "Legend") +
    geom_label_repel(
        data=. %>% mutate(label = ifelse((Bonferroni == T) | (FDR == T),
        as.character(Label), '')), aes(label = label), size = 3.5,
        label.size = NA, box.padding = 0.6, min.segment.length = 0.25) +
    ggtitle("GBA variants associated with Parkinson's Disease") +
    ylab(TeX("$-log_{10}p$")) +
    xlab("Chr. 1") +
    labs(shape = "FDR", color = "") +
    theme(panel.grid.major.x = element_blank(),
          legend.key = element_rect(fill = "transparent"),
          legend.title = element_blank(),
          panel.background = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_line(color = "lightgrey",
                                            linetype = "dotted"),
          panel.grid.minor.y = element_line(color = "lightgrey",
                                            linetype = "dotted"),
          axis.line.y.left = element_line(color = "black"),
          axis.line.x.bottom = element_line(color = "black"))
    
plot

summary <- summary %>%
    mutate(sigOR = case_when(
        Bonferroni == T | FDR == T ~ 1,
        .default = 0
    ))
summary[nrow(summary), 1] <- " Any Variant" # Make it show up last on plot
plot2 <- ggplot(summary[which(summary$sigOR == 1),],
                aes(x = exp(Estimate), y = Label)) +
    geom_point() +
    geom_errorbarh(aes(xmax = exp(Estimate) + 1.96 * exp(`Std. Error`),
                       xmin = exp(Estimate) - 1.96 * exp(`Std. Error`),
                       height = 0.3)) +
    xlim(c(-2, 160)) +
    ylab("Variant") +
    xlab("Odds Ratio Estimate") +
    ggtitle("Parkinson's odds ratios of selected variants with 95% CI") +
    theme(panel.grid.major.x = element_line(color = "lightgrey",
                                            linetype = "dotted"),
          legend.key = element_rect(fill = "transparent"),
          legend.title = element_blank(),
          panel.background = element_blank(),
          axis.ticks.y = element_blank(),
          #axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.line.y.left = element_line(color = "black"),
          axis.line.x.bottom = element_line(color = "black"))
plot2

### Export for PheWAS and other analysis
to_export <- plinkPath %>%
    dplyr::select(IID, summary$Label[which(summary$FDR == T &
                                        summary$Label != " Any Variant")])

anyVar <- pheno %>% dplyr::select(-c(AnyCardio, AnyHemat, AnyHepat, AnyMusc,
                                     AnyNeuro, AnyOcular)) %>%
    mutate(AnyVar = case_when(
        IID %in% AnyVars$IID ~ 1,
        .default = 0
    )) %>%
    dplyr::select(IID, AnyVar)
to_export <- merge(to_export, anyVar)

for (i in 2:ncol(to_export)) {
    new <- to_export %>% select(IID, names(to_export)[i])
    name <- names(new)[2] %>% gsub(":", "_", .) %>% gsub(">", "_", .)
    write_csv(new, sprintf("PheWAS\\plink_%s.csv", name), append = F)
}

path_Vars1 <- path_Vars[which(path_Vars$Variant %in%
                                   summary$Label[which(summary$sigOR == 1)]),]
write_csv(path_Vars1, "SelectedVariants.csv")
write_csv(summary, "AnalysisSummary.csv")

save.image("AfterAnalysis.RData")

### Plot PCA
# Plot those having any variant on PCA
plot3 <- ggplot(any_x[which(any_x$AnyVar != 1),], aes(x = p22009_a1,
                                                      y = p22009_a2)) +
    geom_point(alpha = 0.1, shape = 16, size = 1) +
    geom_point(data = any_x[which(any_x$AnyVar == 1),],
               aes(x = p22009_a1, y = p22009_a2),
               color = "red", alpha = 0.3, shape = 16, size = 1)+
    theme_minimal()+
    xlab("PC1")+
    ylab("PC2")+
    ggtitle("Principal components of those having any GBA variant")
plot3

# Create new variable highlighting selected variants
park <- park %>%
    mutate(WhichVar = case_when(
        rs80356771 == 1 ~ "rs80356771",
        `1:155235769:G>A` == 1 ~ "1:155235769:G>A",
        rs76763715 == 1 | rs76763715 == 0 ~ "rs76763715",
        rs75548401 == 1 | rs75548401 == 0 ~ "rs75548401",
        `1:155236366:C>T` == 1 ~ "1:155236366:C>T",
        rs381427 == 1 ~ "rs381427",
        rs381418 == 1 ~ "rs381418",
        rs61748906 == 1 ~ "rs61748906",
        rs1671872221 == 1 ~ "rs1671872221",
        .default = "None"
    ))
park$WhichVar <- factor(park$WhichVar)

# Plot selected variants on PCA
plot4 <- ggplot(park, aes(x = p22009_a1, y = p22009_a2)) +
    geom_point(alpha = 0.1, shape = 16, size = 1) +
    geom_point(data = park[which(park$WhichVar != "None"),],
               aes(x = p22009_a1, y = p22009_a2,
                   color = WhichVar), alpha = 0.5, shape = 16, size = 1) +
    theme_minimal()+
    xlab("PC1") +
    ylab("PC2")+
    ggtitle("Principal components of selected variants")+
    scale_color_discrete(name = "Variant")
plot4
