library(tidyverse)
library(ggrepel)
library(latex2exp)



### Set up workspace ----
setwd("Z:/UKB Research/GBA1/")
pathVars <- read_csv("pathVars.csv")
filelist <- list.files(path = "PheWAS/csv/")

### Pull p-val from files ----
summary <- data.frame()
for (file in filelist) {
    currentfile <- read_csv(paste("PheWAS/csv/", file, sep = ""))
    park <- currentfile[which(currentfile$phenotype == "332"),]
    if (substr(file, 1, 1) == "1") {
        snp <- strsplit(file, split = "_")
        snp <- paste(snp[[1]][1], ":", snp[[1]][2], ":",
                     snp[[1]][3], ">", snp[[1]][4], sep = "")
        park$snp <- snp
    }
    newline <- merge(park, pathVars,
                     by.x = "snp", by.y = "Variant")
    summary <- rbind(summary, newline)
}
summary <- summary %>%
    rename(pval = p) %>%
    rename(Variant = snp)

bon <- (0.05 / 214) / nrow(currentfile)

# Fix SNP name
summary$Variant[which(summary$Variant == "1:155236366:C>T")] <- "rs1064648"

### Write to file ----
write_csv(summary, "ManhattanSummary.csv")

### Plot ----
# Manhattan plot
ggplot(summary, aes(x = Pos, y = -log10(pval))) +
    geom_point() +
    geom_hline(aes(yintercept = -log10(bon),
                   color = "Bonferroni adjusted alpha",
                   linetype = "Bonferroni adjusted alpha")) +
    scale_color_manual(values = c("Bonferroni adjusted alpha" = "red"),
                       name = "Legend") +
    scale_linetype_manual(values = c("Bonferroni adjusted alpha" = "solid"),
                          name = "Legend") +
    geom_label_repel(
        data = . %>% mutate(label = ifelse(-log10(pval) > 5,
                                           as.character(Variant), "")),
        aes(label = label), size = 3.5,
        label.size = NA, box.padding = 0.6, min.segment.length = 0.25) +
    ylab(TeX("$-log_{10}p$")) +
    xlab("Chromosome 1 Position") +
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
# Save
ggsave("ParkinsonManhattan.png", dpi = 600)
