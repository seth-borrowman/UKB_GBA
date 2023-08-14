library(ggrepel)
library(PheWAS)

setwd("Z:/UKB Research/GBA1/PheWAS/")

pheno <- read.csv('pheno_icd10_long.csv')
covariate <- read.csv('covariates.csv')
phenotypes <- createPhenotypes(
    pheno, min.code.count = 1, add.phecode.exclusions = T, translate = T,
    vocabulary.map = PheWAS::phecode_map_icd10
)
rsid_files <- list.files(path = 'Z:/UKB Research/GBA1/PheWAS/',
                        pattern = glob2rx('plink*csv'))

# Create table to aggregate significant results from all variants
sig_phewas = data.frame()

# Create directories to write PheWAS results and plots for each variant
ifelse(!dir.exists('png'), dir.create('png'), FALSE)
ifelse(!dir.exists('csv'), dir.create('csv'), FALSE)

# Run PheWAS for each variant separately
for (file in rsid_files){
    rsid <- strsplit(file, '.csv')[[1]]
    rsid <- strsplit(rsid, 'plink_')[[1]][2]
    print(rsid)
    # Load geno data
    geno_file <- paste('Z:/UKB Research/GBA1/PheWAS/', file, sep = "")
    geno_data <- read.csv(geno_file)
    names(geno_data)[1] <- "id"
    results <- phewas(
        phenotypes, geno_data, cores = detectCores(),
        significance.threshold = c('p-value', 'bonferroni', 'fdr'),
        covariates = covariate, additive.genotypes = T,
        MASS.confint.level = 0.95)
    
    # Add PheWAS descriptions
    results_d <- addPhecodeInfo(results)
    
    # Get significant results
    res <- results_d[results_d$bonferroni & !is.na(results_d$p),]
    print("Results")
    print(res)
    sig_phewas <- rbind(sig_phewas, res)
    
    # Re-create the same threshold for plots as is used in bonferroni column
    sig_p <- 0.05/(nrow(results_d[!is.na(results_d$p),]))
    print("sig_p")
    print(sig_p)
    print(nrow(results_d[!is.na(results_d$p),]))
    
    # Exctract significant result and save it as csv
    results_d <- results_d[!is.na(results_d$p),]
    results_d <- results_d[order(results_d$group),]
    results_d$order_num <- rep(1:nrow(results_d))
    results_d$order_num <- factor(results_d$order_num, levels = results_d$order_num)
    write.csv(results_d, sprintf('csv/%s_phewas.csv', rsid), row.names = FALSE)
    
    # Create Manhattan plot annotating significant phenotypes
    # Significant phenotypes are defined by passing Bonferroni correction
    png(filename=sprintf('png/%s_phewas.png', rsid), , width = 1400, height = 800)
    options(ggrepel.max.overlaps = Inf) 
    man_plot <- ggplot(
        results_d, 
        aes(x=order_num, y=-log(p))) + geom_point(aes(col=group, size=OR)) + theme_classic() + theme(
            axis.text.x = element_blank(), panel.grid.minor=element_line(colour = 'grey', linetype = 'dashed'), 
            axis.ticks=element_blank()
        ) + labs(color = 'Category', size = 'Effect size', x = rsid, y = '-log(p-value)') + geom_text_repel(
            data=. %>% mutate(label = ifelse((p < sig_p) & (bonferroni == TRUE), as.character(description), '')), aes(label = label), size = 3,
            box.padding = unit(0.7, 'lines')
        ) + geom_hline(yintercept=-log(sig_p), color='red', linewidth=1, alpha=0.5) 
    print(man_plot)
    dev.off()
}

write.csv(sig_phewas, 'significant_phewas.csv', row.names = FALSE)

sig_phewas_tb <- tibble::tibble(sig_phewas)
sig_phewas_agg <- dplyr::count(sig_phewas_tb, description, sort = TRUE)

write.csv(sig_phewas_agg, 'significant_phewas_agg.csv', row.names = FALSE)