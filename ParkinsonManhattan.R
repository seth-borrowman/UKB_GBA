library(tidyverse)
library(ggrepel)
library(latex2exp)
library(logistf)

setwd("~/LSS/lss_chndr/UKB Research/GBA1")

### Set up workspace ----
covariate <- read_csv("PheWAS/covariates.csv")
plinkPath <- read_csv("plinkPath.csv")
phenotypes <- read_csv("PheWAS/phecode.csv")
pathVars <- read_csv("pathVars.csv")
phenotypes <- phenotypes[,-1]
excl_codes <- PheWAS::phecode_exclude$exclusion_criteria[which(PheWAS::phecode_exclude$code == 332)]
excl <- which(colnames(phenotypes) %in% excl_codes)

### Set exclusions to match PheWAS package----
excl_phecode <- phenotypes[,c(1, excl)] %>% filter_all(any_vars(. %in% c(T)))
phenotypes <- phenotypes %>%
    mutate(Parkinson = case_when(
        `332` == T  ~ 1,
        .default = 0
    )) %>%
    mutate(Parkinson = case_when(
        id %in% excl_phecode$id ~ 0,
        .default = Parkinson
    )) %>%
    select(c(id, Parkinson)) %>%
    rename(IID = id)

summary <- data.frame()

### Create functions for tryCatch ----
tryfun <- function(data, var_name) {
    mod <- logistf(Parkinson ~ ., data = data,
                   control = logistf.control(maxstep = -1))
    beta <- unname(mod$coefficients[24])
    beta_ci_low <- unname(mod$ci.lower[24])
    beta_ci_upper <- unname(mod$ci.upper[24])
    OR <- exp(beta)
    OR_ci_low <- exp(beta_ci_low)
    OR_ci_upper <- exp(beta_ci_upper)
    pval <- unname(mod$prob[24])
    sum_line <- c(var_name, beta, beta_ci_low, beta_ci_upper,
                  OR, OR_ci_low, OR_ci_upper, pval)
    return(sum_line)
}
errorfun <- function(var_name) {
    sum_line <- c(var_name, rep(1, 7))
    return(sum_line)
}

### Test association of each variant ----
for (i in 2:ncol(plinkPath)) {
    genotype <- plinkPath[,c(1, i)]
    varname <- colnames(genotype)[2]
    message(paste(varname,
                  as.character(as.integer(i) - 1),
                  "/", ncol(plinkPath) - 1))
    colnames(genotype)[2] <- "Variant"
    genotype <- genotype %>%
        mutate(Variant = case_when(
            Variant == 2 ~ 0,
            Variant == 1 | Variant == 0 ~ 1
        )) %>%
        mutate(Variant = factor(Variant))
    df <- merge(covariate, phenotypes) %>%
        merge(., genotype) %>%
        select(-IID) %>%
        mutate(Sex = factor(Sex)) %>%
        mutate(Alcohol = factor(Alcohol)) %>%
        mutate(TBI = factor(TBI)) %>%
        mutate(Parkinson = factor(Parkinson))
    message("Fitting model")
    # Some have had trouble converging - use try/catch
    tryCatch(sum_line <- tryfun(df, varname),
             error = (sum_line <- errorfun(varname)))
    print(sum_line)
    summary <- rbind(summary, sum_line)
}
# Rename summary columns something useful
colnames(summary) <- c("variant", "beta", "beta_ci_low", "beta_ci_upper",
                       "OR", "OR_ci_low", "OR_ci_upper", "pval")

### Set significance cutoffs using Bonferroni and Bejamini-Hochberg FDR ----
summary <- summary %>%
    mutate(Bonferroni = if_else(pval <= (0.05/nrow(summary)), T, F)) %>%
    mutate(FDR = if_else(pval <= ((match(pval,
                                         summary$pval[order(summary$pval, decreasing = F)]) /
                                       nrow(summary)) * 0.05), T, F))
max_FDR_p <- max(summary[which(summary$FDR == T), 8])
summary <- summary %>%
    mutate(FDR = if_else(summary$pval <= max_FDR_p, T, F))

### Write to file ----
write_csv(summary, "ManhattanSummary.csv")

### Plot ----
# Get variant position
summary <- merge(summary, pathVars[,c(1,3)], by.x = "variant", by.y = "Variant")
# Manhattan plot
ggplot(summary, aes(x = Pos, y = -log10(pval))) +
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
        data = . %>% mutate(label = ifelse((Bonferroni == T) | (FDR == T),
                                           as.character(variant), "")), aes(label = label), size = 3.5,
        label.size = NA, box.padding = 0.6, min.segment.length = 0.25) +
    ggtitle("GBA1 variants associated with Parkinson's disease") +
    ylab(TeX("$-log_{10}p$")) +
    xlab("Chromosome 1") +
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
ggsave("ParkinsonManhattan.png")
