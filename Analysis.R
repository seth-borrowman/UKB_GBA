library(tidyverse)
library(latex2exp)
library(ggrepel)

setwd("Z:/UKB Research/GBA1")
load("CleanVariantsWrkspc.RData")
pheno <- read_csv("Phenotypes.csv")

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
try(plinkPath <- plinkPath[,-removecols], silent = T)


### Parkinson's

park <- merge(pheno, plinkPath, by = "IID") %>%
    select(-c(AnyCardio, AnyHemat, AnyHepat, AnyMusc, AnyNeuro, AnyOcular, IID))
 
y_park <- as.matrix(as.numeric(park$Parkinson)) - 1
x_park <- park %>% select(-Parkinson)
x_park <- model.matrix(~ . -1, data = x_park)
# Keeps female in there for some reason - need to remove
x_park <- x_park[,-which(colnames(x_park) == "SexFemale")]

# Association analysis with logistic regression for each SNV
summary <- matrix(rep(0, 5*(ncol(x_park)-14)), # -14 to skip covariates
                  ncol = 5, nrow = ncol(x_park)-14) %>%
    as.data.frame()
colnames(summary) <- c("Variant", "Estimate", "Std. Error", "z value",
                       "p")
for (i in 15:ncol(x_park)) {
    variants <- seq(15, ncol(x_park))
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
any_x <- pheno %>% select(-c(AnyCardio, AnyHemat, AnyHepat, AnyMusc, AnyNeuro,
                             AnyOcular)) %>%
    mutate(AnyVar = case_when(
        IID %in% AnyVars$IID ~ 1,
        .default = 0
    ))
mod_any <- glm(Parkinson ~ ., data = any_x, family = "binomial")
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
    select(IID, summary$Label[which(summary$FDR == T &
                                        summary$Label != " Any Variant")])

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
