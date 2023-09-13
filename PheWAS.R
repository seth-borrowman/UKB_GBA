library(ggrepel)
library(PheWAS)

### Change PheWAS package's logistic regression to Firth adjusted----
phe_as <-
    function(phe.gen, additive.genotypes=T,min.records=20,return.models=F,confint.level=NA, my.data, my.covariates) {
        if(!missing(my.data)) data=my.data
        if(!missing(my.covariates)) covariates=my.covariates
        #Retrieve the targets for this loop
        phe_o=phe.gen[[1]]
        phe=phe_o
        gen=phe.gen[[2]]
        gens=gen
        adjustment=phe.gen[[3]]
        #Subset the data
        d=data[,na.omit(unlist(c(gen,phe,covariates,adjustment)))]
        #Turn adjustment into a string, if not NA
        if(!is.na(adjustment[1])) {adjustment=paste(adjustment,collapse=",")}
        else {adjustment=NA_character_} #Make sure it is a character NA for aggregation
        #Alter phe_o if necessary for the regression formula
        if(suppressWarnings(!is.na(as.numeric(phe_o)))) {
            phe=paste0("pheno_",phe_o)
            names(d)[2]=phe
        }
        #Exclude the exclusions for the target phenotype
        d=d[!is.na(d[[phe]]),]
        n_no_snp=sum(is.na(d[[gen]]))
        #Exclude rows with missing data
        d=na.omit(d)
        n_total=nrow(d)
        n_cases=NA_integer_
        n_controls=NA_integer_
        allele_freq=NA_real_
        HWE_pval=NA_real_
        or=NA_real_
        se=NA_real_
        p=NA_real_
        beta=NA_real_
        type=NA_character_
        note=""
        model=NA
        model_name=sprintf("No model: %s ~ %s",phe,paste0(names(d)[-1],collapse = " + "))
        if(n_total<min.records) {
            note=paste(note,"[Error: <",min.records," complete records]")
        } else if(length(unique(na.omit(d[[phe]])))<=1 | length(unique(na.omit(d[[gen]]))) <=1) {
            note=paste(note,"[Error: non-varying phenotype or genotype]")
        } else {
            if(additive.genotypes) {
                if(class(d[[gen]]) %in% c("numeric","integer")){
                    allele_freq=sum(d[[gen]])/(2*n_total)
                }
                if(class(d[[gen]]) %in% c("numeric","integer") & sum(!(na.omit(d[[gen]]) %in% 0:2))==0) {
                    P=allele_freq
                    Q=1-allele_freq
                    AA=sum(d[[gen]]==2)
                    xAA=P^2*n_total
                    Aa=sum(d[[gen]]==1)
                    xAa=2*P*Q*n_total
                    aa=sum(d[[gen]]==0)
                    xaa=Q^2*n_total
                    HWE_pval=pchisq((AA-xAA)^2/(xAA)+(Aa-xAa)^2/(xAa)+(aa-xaa)^2/(xaa),1)
                } else {note=paste(note,"[Warning: Genotype is not coded 0,1,2, but additive.genotypes was TRUE.]")}
            } 
            #Check if genotype was available
            #Check if phenotype is logical (boolean)
            if(class(d[[phe]]) %in% c("logical")) {
                type = "logistic"
                #Create the logistic model
                n_cases=sum(d[[phe]])
                n_controls=n_total-n_cases
                if(n_cases<min.records|n_controls<min.records) {note=paste(note,"[Error: <",min.records," cases or controls]")}
                else {
                    model = logistf::logistf(as.formula(paste(phe," ~ .")), data=d, dataout=F, alpha=confint.level)
                    modsum= summary(model)
                    model_name=paste0(as.character(terms(model))[c(2,1,3)],collapse=" ")
                    #If the models did not converge, report NA values instead.
                    if(model$converged) {
                        #Find the observed genotype columns
                        gen_expansion=attr(model.matrix(my.formula, data=d),"assign")
                        gen_list=which(gen_expansion %in% 1:length(gen))
                        gen_expansion=gen_expansion[gen_list]
                        
                        #Find the rows with results that gets merged across all loops
                        gens=names(model$coef)[gen_list]
                        or=exp(model$coef[gen_list])
                        beta=model$coef[gen_list]
                        se=(diag(model$var)^0.5)[gen_list]
                        p=model$prob[gen_list]
                    } else {
                        note=paste(note,"[Error: The model did not converge]")
                    }
                }
            } else {
                type = "linear"
                if(n_total<min.records) {
                    note=paste(note,"[Error: <",min.records," records with phenotype and genotype]")
                } else {
                    model = glm(as.formula(paste(phe," ~ .", sep="", collapse="")), data=d)
                    modsum= summary(model)
                    model_name=paste0(as.character(terms(model))[c(2,1,3)],collapse=" ")
                    
                    #If the models did not converge, report NA values instead.
                    if(model$converged) {
                        #Find the rows with results that gets merged across all loops
                        gen_list=grep(gen,row.names(modsum$coef))
                        gens=row.names(modsum$coef)[gen_list]
                        beta=modsum$coef[gen_list,1]
                        se=modsum$coef[gen_list,2]
                        p=modsum$coef[gen_list,4]
                    } else {
                        note=paste(note,"[Error: The model did not converge]")
                    }
                }
            }
        }
        
        output=data.frame(phenotype=phe_o,snp=gens,
                          adjustment=adjustment,
                          beta=beta, SE=se,
                          OR=or,
                          p=p, type=type,
                          n_total=n_total, n_cases=n_cases, n_controls=n_controls,
                          HWE_p=HWE_pval,allele_freq=allele_freq,n_no_snp=n_no_snp, 
                          note=note, stringsAsFactors=F)
        
        #Add confidence intervals if requested.
        if(!is.na(confint.level)) {
            if(!is.na(model)[1]){
                suppressMessages(conf<-confint(model,c("(Intercept)",gens),level=confint.level))
                lower=conf[-1,1]
                upper=conf[-1,2]
                if(type=="logistic") {
                    lower=exp(lower)
                    upper=exp(upper)
                }
            } else {
                lower=NA_real_
                upper=NA_real_
            }
            output$lower=lower
            output$upper=upper
            
            output=output[,c("phenotype","snp","adjustment","beta","SE",
                             "lower","upper","OR","p","type",
                             "n_total","n_cases","n_controls",
                             "HWE_p","allele_freq","n_no_snp","note")]
        }
        
        #If the complete models were requested, add them as well.
        if(return.models) {
            attributes(output)$model=model
            attributes(output)$model_name=model_name
        }
        attributes(output)$successful.phenotype=ifelse(is.na(p),NA,phe_o)
        attributes(output)$successful.genotype=ifelse(is.na(p),NA,gen)
        #Return this to the loop to be merged.
        output
    }

### Import data ----
setwd("Z:/UKB Research/GBA1/PheWAS/")
pheno <- read.csv('pheno_icd10_long.csv')
covariate <- read.csv('covariates.csv')
# Make things easier for creating phenotypes
names(covariate)[1] <- "id"
sex <- covariate %>%
    select(id, Sex) %>%
    mutate(Sex = case_when(
        Sex == "Male" ~ "M",
        Sex == "Female" ~ "F"
    ))
# Remove sex from further analysis - causes issues
covariate <- covariate %>% select(-Sex)
### Create phenotypes ----
phenotypes <- createPhenotypes(
    pheno, min.code.count = 1, add.phecode.exclusions = T, translate = T,
    vocabulary.map = PheWAS::phecode_map_icd10,
    full.population.ids = covariate$id,
    id.sex = sex)
dplyr::write_csv(phenotypes, "phecode.csv")

# Set up folders ----
rsid_files <- list.files(path = getwd(),
                        pattern = glob2rx('plink*csv'))

# Create table to aggregate significant results from all variants
sig_phewas = data.frame()

# Create directories to write PheWAS results and plots for each variant
ifelse(!dir.exists('png'), dir.create('png'), FALSE)
ifelse(!dir.exists('csv'), dir.create('csv'), FALSE)

# Run PheWAS for each variant ----
for (file in rsid_files) {
    print(Sys.time())
    rsid <- strsplit(file, '.csv')[[1]]
    rsid <- strsplit(rsid, 'plink_')[[1]][2]
    print(rsid)
    # Load geno data
    geno_file <- paste('Z:/UKB Research/GBA1/PheWAS/', file, sep = "")
    geno_data <- read.csv(geno_file)
    names(geno_data)[1] <- "id"
    geno_data[,2] <- -(geno_data[,2] - 2)
    results <- phewas(
        phenotypes, geno_data, cores = detectCores(),
        significance.threshold = c('p-value', 'bonferroni', 'fdr'),
        covariates = covariate, additive.genotypes = T,
        alpha = (0.05/length(rsid_files)))
    
    # Add PheWAS descriptions
    results_d <- addPhecodeInfo(results)
    
    # Get significant results
    res <- results_d[results_d$bonferroni & !is.na(results_d$p),]
    print("Results")
    print(res)
    sig_phewas <- rbind(sig_phewas, res)
    
    # Re-create the same threshold for plots as is used in bonferroni column
    sig_p <- (0.05/length(rsid_files))/(nrow(results_d[!is.na(results_d$p),]))
    print("sig_p")
    print(sig_p)
    print(nrow(results_d[!is.na(results_d$p),]))
    
    # Exctract significant result and save it as csv
    results_d <- results_d[!is.na(results_d$p),]
    results_d <- results_d[order(results_d$group),]
    # Add try() because this step sometimes fails and kills the loop
    try(results_d$order_num <- rep(1:nrow(results_d)), silent = T)
    try(results_d$order_num <- factor(results_d$order_num,
                                  levels = results_d$order_num), silent = T)
    write.csv(results_d, sprintf('csv/%s_phewas.csv', rsid), row.names = FALSE)
    
    # Create Manhattan plot annotating significant phenotypes
    # Significant phenotypes are defined by passing Bonferroni correction
    png(filename = sprintf('png/%s_phewas.png', rsid),
        width = 1400, height = 800)
    options(ggrepel.max.overlaps = Inf) 
    man_plot <- ggplot(
        results_d, 
        aes(x = order_num, y = -log(p))) +
        geom_point(aes(col = group, size = OR)) +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            panel.grid.minor = element_line(colour = 'grey',
                                            linetype = 'dashed'), 
            axis.ticks = element_blank()) +
        labs(color = 'Category', size = 'Effect size',
             x = rsid, y = '-log(p-value)') +
        geom_text_repel(
            data = . %>%
                mutate(label = ifelse((p < sig_p) & (bonferroni == TRUE),
                                      as.character(description), '')),
            aes(label = label), size = 3, box.padding = unit(0.7, 'lines')) +
        geom_hline(yintercept = -log(sig_p), color = 'red',
                   linewidth = 1, alpha = 0.5) 
    print(man_plot)
    dev.off()
}

# Output aggregate results ----
write.csv(sig_phewas, 'significant_phewas.csv', row.names = FALSE)

sig_phewas_tb <- tibble::tibble(sig_phewas)
sig_phewas_agg <- dplyr::count(sig_phewas_tb, description, sort = TRUE)

write.csv(sig_phewas_agg, 'significant_phewas_agg.csv', row.names = FALSE)
