library(tidyverse)

setwd("Z:/UKB Research/GBA1")

# Import data
var <- read_delim("variant_summary.txt", delim = "\t")
plink <- read_delim("GBA1exome.raw", delim = "\t")

# Remove plink columns that are all NA or there is no variation
plink <- Filter(function(x)!all(is.na(x)), plink)
removecols <- c()
for (i in 1:ncol(plink)) {
    if (min(plink[,i], na.rm = T) == max(plink[,i], na.rm = T)) {
        removecols <- append(removecols, as.integer(i))
    }
}
plink <- plink[,-removecols]
rm(removecols)

# Limit ClinVar reference to Chr1 and correct assembly
varch1 <- var[which(var$Chromosome == 1 & var$Assembly == "GRCh38"),]
rm(var)
gc()

# Create new df with only relevant columns
reference <- varch1 %>%
    select(c(Chromosome, Start, ReferenceAlleleVCF, AlternateAlleleVCF,
             `RS# (dbSNP)`, ClinicalSignificance))

# Create names from ClinVar to match Plink output
reference$PlinkName <- paste(reference$Chromosome, ":", reference$Start, ":",
                             reference$ReferenceAlleleVCF, ":",
                             reference$AlternateAlleleVCF, "_",
                             reference$ReferenceAlleleVCF, sep = "")
# Match up Plink names w/ ClinVar
plinknames <- names(plink)
plinknames <- plinknames[-c(1:3)] %>%
    as.data.frame()
colnames(plinknames) <- "PlinkName"

names <- left_join(plinknames, reference, by = "PlinkName")

# Redo Plink names in a way I like
for (i in 1:nrow(names)) {
    a <- strsplit(names$PlinkName[i], split = ":")
    b <- strsplit(a[[1]][4], split = "_")
    a[[1]][4] <- b[[1]][1]
    a <- paste(a[[1]][1], ":", a[[1]][2], ":",
               a[[1]][3], ">", a[[1]][4], sep = "")
    names$PlinkName[i] <- a
}

# Create new name for variants choosing ClinVar rs# over Plink format
for (i in 1:nrow(names)) {
    if (is.na(names$`RS# (dbSNP)`[i]) | names$`RS# (dbSNP)`[i] == -1) {
        names$NewName[i] <- names$PlinkName[i]
    } else {
        names$NewName[i] <- paste("rs", names$`RS# (dbSNP)`[i], sep = "")
    }
}
names$NewName <- make.unique(names$NewName)

# Make plink have new names
a <- names$NewName
a <- append(a, c("FID", "IID", "Sex"), after = 0)
colnames(plink) <- make.unique(a) # adds .# to duplicate column names
plink <- select(.data = plink, -c(FID, Sex))

# Clean up
rm(b)
rm(plinknames)
rm(reference)
rm(varch1)
rm(a)
rm(i)
gc()

# Import REVEL predictions
revel <- read_csv("revel_grch38_chrom_01_152188850_156640555.csv")

# Get positions for all variants
for (i in 1:nrow(names)) {
    if (is.na(names[i, 2])) {
        names[i, 2] <- 1
        names[i, 3] <- substr(names[i, 1], 3, 11)
        a <- strsplit(names[i, 1], split = ":")
        b <- strsplit(a[[1]][3], split = ">")
        names[i, 4] <- b[[1]][1]
        names[i, 5] <- b[[1]][2]
    }
}

# Join ClinVar and REVEL predictions
names$Start <- as.double(names$Start)
pathogenicity <- left_join(names,
                           revel,
                           by = c("Start" = "grch38_pos",
                                  "ReferenceAlleleVCF" = "ref",
                                  "AlternateAlleleVCF" = "alt"))

# Clean up pathogenicity prediction table
pathogenicity <- pathogenicity %>%
    select(NewName, Chromosome, Start, ReferenceAlleleVCF, AlternateAlleleVCF,
           aaref, aaalt, ClinicalSignificance, REVEL)
colnames(pathogenicity) <- c("Variant", "Chr", "Pos", "Ref", "Alt",
                             "Ref_AA", "Alt_AA", "ClinSig", "REVEL")
pathogenicity$ClinSig <- factor(pathogenicity$ClinSig)

# Clean up
rm(a)
rm(b)
rm(names)
rm(revel)
rm(i)
gc()

# Filter to only pathogenic variants
table(pathogenicity$ClinSig)
pathVars <- pathogenicity[which(pathogenicity$ClinSig == "Likely pathogenic" |
                        pathogenicity$ClinSig == "Pathogenic" |
                        pathogenicity$ClinSig == "Pathogenic/Likely pathogenic" |
                        pathogenicity$ClinSig == "Pathogenic/Likely pathogenic; risk factor" |
                        pathogenicity$REVEL > 0.6),]

# Filter plink and factor so homozyg. major allele is reference
plinkPath <- cbind(plink$IID, plink[,which(names(plink) %in% pathVars$Variant)])
names(plinkPath)[1] <- "IID"
plinkPath[,2:ncol(plinkPath)] <- plinkPath[, 2:ncol(plinkPath)] %>%
    lapply(factor, levels = c(2, 1, 0))
plinkPath <- na.omit(plinkPath)


# Find allele freq of minor variants in cohort
pathVars$CohortFreq <- rep(0, nrow(pathVars))
for (i in 1:nrow(pathVars)) {
    pathVars$CohortFreq[i] <- (0.5 * unname(table(plinkPath[,(i+1)])[2]) +
        unname(table(plinkPath[,(i+1)])[2])) / nrow(plinkPath)
}

# Count how many people have each variant
pathVars$NumHet <- rep(-1, nrow(pathVars))
pathVars$NumHom <- rep(-1, nrow(pathVars))
for (i in 2:ncol(plinkPath)) {
    pathVars$NumHet[i-1] <- length(which(plinkPath[,i] == "1"))
    pathVars$NumHom[i-1] <- length(which(plinkPath[,i] == "0"))
}
pathVars <- pathVars %>% mutate(Total = NumHet + NumHom)

# Count how many people in the cohort have any variants
AnyVars <- plinkPath %>% filter(if_any(where(is.factor), \(val) val %in% c("0", "1")))
    # 22301
for(i in 2:ncol(plinkPath)) {
    if(names(plinkPath)[i] != pathVars$Variant[i-1]){
        print(names(plinkPath)[i])
    }
}
rm(plink)
rm(pathogenicity)
save.image("Z:/UKB Research/GBA1/CleanVariantsWrkspc.RData")
