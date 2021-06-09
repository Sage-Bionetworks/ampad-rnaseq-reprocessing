###### MSBB Meta Data Cleaning and Counts combining
# Metadata and counts for sageseqr pipeline are assembled and 
# pushed to synapse
library(readxl)
library(synapser)
library(data.table)
library(dplyr)
library(tidyr)
library(digest)
library(stringr)

counts <- 'syn21544664';
used_synids <- counts
counts <- synapser::synGet(counts)$path %>%
  read.table(header=T, sep='\t', check.names = F, row.names = 1)

counts <- counts[ row.names(counts)[ 
  !(row.names(counts) %in% 
      c('N_unmapped','N_multimapping','N_noFeature','N_ambiguous')
  )
],]

#Biospecimin and assay don't align with hB_RNA_10892
# rename to hB_RNA_10892_resequenced
colnames(counts)[colnames(counts) == 'hB_RNA_10892'] <-
  'hB_RNA_10892_resequenced'

# Convert rownames of counts from tracking id to ensemble gene id
# tmp = data.frame(Gene.ID = rownames(counts)) %>%
#   dplyr::mutate(ID = Gene.ID) %>%
#   tidyr::separate(ID, c('ensembl_gene_id', 'position'), sep = '\\.')
#  rownames(tmp) = tmp$Gene.ID
# rownames(counts) = tmp[rownames(counts), 'ensembl_gene_id']

# Get sample ids
# SampleID = data.frame(SampleID = colnames(counts),
#                       ID = colnames(counts)) %>% 
#   tidyr::separate(ID, c('A','B','ID'), sep = '_') %>%
#   dplyr::select(SampleID, ID) %>%
#   dplyr::mutate(ID = as.numeric(ID))


# Get technical metrics
metrics <- 'syn21544666'
used_synids <- c(used_synids,metrics)
metrics <- synapser::synGet(metrics)$path %>%
  read.table(sep='\t',header=T, row.names=1, stringsAsFactors = F)
colnames(metrics) <- gsub('__', '_', colnames(metrics))
colnames(metrics)[colnames(metrics) == 'sample'] <- 'specimenID'
metrics$specimenID <- row.names(metrics)

# Get 90+ age
metadata_age = 'syn10156693'
used_synids <- c(used_synids,metadata_age)
metadata_age <- synGet(metadata_age)$path %>%
  read.table(sep='\t', header=T, stringsAsFactors = F)

# Get clinical metadata
clinical <- 'syn6101474'
used_synids <- c(used_synids,clinical)
clinical <- synGet(clinical)$path %>%
  read.csv(header = T, stringsAsFactors = F)

# RECODE THE CERAD TO MATCH ROSMAP
# MSBB ####################
## 1 = Normal
## 2 = Definite AD
## 3 = Probable AD
## 4 = Possible AD
###########################
# RosMap ##################
## 1 = Definite AD
## 2 = Probable AD
## 3 = Possible AD
## 4 = No AD/Normal
## Cases (CERAD =c(1,2))
## Controls (CERAD =c(3,4))
###########################
clinical$ceradsc <- 0
clinical[ clinical$CERAD == 1, ]$ceradsc <- 4
clinical[ clinical$CERAD == 2, ]$ceradsc <- 1
clinical[ clinical$CERAD == 3, ]$ceradsc <- 2
clinical[ clinical$CERAD == 4, ]$ceradsc <- 3
clinical <- clinical[,colnames(clinical)[!(colnames(clinical) %in% 'CERAD')]]

#########################################
#Recoad ethnicity to match  ROSMAP:
#1=White
#2=Black, Negro, African-American
#3=Native American, Indian
#4=Eskimo 
#5=Aleut
#6=Asian or Pacific Island
#8=REFUSAL
#9=DON'T KNOW
##########################################
#MSBB:
# A: Asian
# B: Black
# H: Hispanic
# U: Unknown
# W: White
##########################################
clinical[ clinical$race == 'W', ]$race <- 1
clinical[ clinical$race == 'B', ]$race <- 2
clinical[ clinical$race == 'A', ]$race <- 6
clinical[ clinical$race == 'U', ]$race <- 9
clinical$spanish <- 2
clinical[ clinical$race == 'H', ]$spanish <- 1
clinical[ clinical$race == 'H', ]$race <- 9

##### Code Diagnosis #####
# Diagnosis
clinical$diagnosis <- 'OTHER'
clinical[ clinical$braaksc >= 4 &
            clinical$ceradsc <= 2 &
            clinical$CDR >= 1  , ]$diagnosis <- 'AD'
clinical[ clinical$braaksc <= 3 &
            clinical$ceradsc >= 3 &
            clinical$CDR <= 0.5  , ]$diagnosis <- 'CT'


#Braaksc/Braak
colnames(clinical)[colnames(clinical) == 'Braak'] <- 'braaksc'

# apoe Allele:
clinical$apoe4_allele <- NA
clinical[ clinical$apoeGenotype %in% c(22,23,33),]$apoe4_allele <- 0
clinical[ clinical$apoeGenotype %in% c(24,34),]$apoe4_allele <- 1
clinical[ clinical$apoeGenotype %in% c(44),]$apoe4_allele <- 2


# conver pmi time interval from minutes to hours
clinical$pmi <- clinical$pmi/60

# Get biospecimin metadata
biospecimin <- 'syn21893059'
used_synids <- c(used_synids,biospecimin)
biospecimin <- synGet(biospecimin)$path %>%
  read.csv(header = T, stringsAsFactors = F)
biospecimin <- biospecimin[ biospecimin$assay=='rnaSeq',]
biospecimin<-biospecimin[biospecimin$specimenID %in% colnames(counts),]
biospecimin$tissue[biospecimin$tissue=='superior temporal gyrus'] <- 'STG'
biospecimin$tissue[biospecimin$tissue=='frontal pole'] <- 'FP'
biospecimin$tissue[biospecimin$tissue=='inferior frontal gyrus'] <- 'IFG'
biospecimin$tissue[biospecimin$tissue=='parahippocampal gyrus'] <- 'PHG'

#Toss sample swaps
biospecimin <- biospecimin[ biospecimin$exclude == 'false',]
counts <- counts[,biospecimin$specimenID]
clinical <- clinical[ clinical$individualID %in% biospecimin$individualID, ]

# Get assay metadata
assay <- 'syn22447899'
used_synids <- c(used_synids,assay)
assay <- synGet(assay)$path %>%
  read.csv(header = T, stringsAsFactors = F)
assay <- assay[, !(colnames(assay) %in% c('totalReads', 'mapped', 'rRNA.rate')) ]
assay$RIN2 <- assay$RIN^2
assay<-assay[assay$specimenID %in% colnames(counts),]

#Merge the metadata frames
metrics <- metrics[metrics$specimenID %in% biospecimin$specimenID,]
comb <- dplyr::right_join(biospecimin,
                          assay[assay$specimenID %in% biospecimin$specimenID,], 
                          by = 'specimenID') %>%
  dplyr::right_join(clinical, by ='individualID') %>%
  dplyr::right_join(metrics, by = 'specimenID')

comb$Tissue.APOE4 <- paste0(  comb$tissue, '.', comb$apoe4_allele)

# Find Tossable IDs:
# Remove the factors that don't matter:
remove_na <- NULL
length_issue <- NULL
all_zeros <- NULL
for(name in colnames(comb)){
  tab <- table(is.na(comb[,name]))
  if('TRUE' %in% names(tab)){
    if( as.numeric(table(is.na(comb[,name]))['TRUE']) == dim(comb)[1]) {
      remove_na <- c(remove_na,name)
    }
  }
  tab <- table(comb[,name])
  if('0' %in% names(tab)){
    if(as.numeric(tab['0']) == dim(comb)[1]) {
      all_zeros <- c(all_zeros,name)
    }
  }
  if(length(tab)==1){
    length_issue <- c(length_issue, name)
  }
}

# Variables that are all NAs
message('Variables that are all NAs:')
message(paste0('Variables that are all NAs: ', 
               paste(remove_na,collapse = ', ')
)
)

# Variables that are all Zeros
message(paste0('Variables that are all Zeros: ',
               paste(all_zeros,collapse = ', ')
)
)

message(paste0('Variables that are all the same value: ',
               paste(length_issue,collapse = ', ')
)
)

toss <- c(remove_na,length_issue,all_zeros)
toss <- toss[!duplicated(toss)]
keeps <- colnames(comb)[ !(colnames(comb) %in% toss)]
keeps_AsmRsm <- keeps[ 
  grepl('AlignmentSummaryMetrics_',keeps) | 
  grepl('RnaSeqMetrics_',keeps) ]
keeps <- keeps[!(keeps %in% keeps_AsmRsm)]

colnames(comb)[ colnames(comb) == 'ageDeath'] <- 'age_death'

total_metadata <- comb[ , c("specimenID", "individualID", "tissue", 
                            "BrodmannArea", "sex", "race","spanish", 
                            "ethnicity", "age_death", "braaksc", "CDR", 
                            "plaqueMean", "ceradsc", "diagnosis", "RIN", "RIN2", 
                            "sequencingBatch", "barcode",  "pmi", "apoeGenotype", 
                            "apoe4_allele", "Tissue.APOE4", keeps_AsmRsm,toss)]
metadata <- comb[ , c("specimenID", "individualID", "tissue", 
                      "BrodmannArea", "sex", "race","spanish", 
                      "ethnicity", "age_death", "braaksc", "CDR", 
                      "plaqueMean", "ceradsc", "RIN", "RIN2", "diagnosis", 
                      "sequencingBatch", "barcode",  "pmi", "apoeGenotype", 
                      "apoe4_allele", "Tissue.APOE4", keeps_AsmRsm)] 

sageseqr_censored <- comb[ , c("specimenID", "individualID", "tissue", 
                               "BrodmannArea", "sex", "race","spanish", 
                               "ethnicity", "age_death", "braaksc", "CDR", 
                               "plaqueMean", "ceradsc", "RIN", "RIN2", "diagnosis",
                               "sequencingBatch", "barcode",  "pmi", "apoeGenotype", 
                               "apoe4_allele", "Tissue.APOE4", 
                               'AlignmentSummaryMetrics_PCT_PF_READS_ALIGNED',
                               'RnaSeqMetrics_PCT_INTRONIC_BASES',
                               'RnaSeqMetrics_PCT_INTERGENIC_BASES',
                               'RnaSeqMetrics_PCT_CODING_BASES')] 

row.names(metadata_age) <- metadata_age$individualIdentifier
sageseqr_uncensored <- sageseqr_censored

sageseqr_uncensored$age_death <- metadata_age[sageseqr_uncensored$individualID,]$AOD
metadata$age_death <- metadata_age[metadata$individualID,]$AOD
total_metadata$age_death <- metadata_age[total_metadata$individualID,]$AOD

counts <- counts[,total_metadata$specimenID]

################################################################################
## --   Push to Synapse   --  ##
#Upload full file and SageSeqr input version to synapse - internal Sage Location:
internal_parentid <- 'syn25872211'
folder_loc <- 'syn25872212'

#Set Activity
activity <- synapser::synGetEntity(folder_loc)

##Set Annotations:
all.annotations = list(
  dataType = c('clinical','geneExpression'),
  resourceType = 'metadata',
  metadataType = 'analytical covariates',
  isModelSystem = 'FALSE',
  isMultiSpecimen = 'TRUE',
  fileFormat = 'csv',
  grant = 'U01AG046152',
  species = 'Human',
  organ = 'brain',
  tissue = c('frontal pole ',
               'inferior frontal gyrus',
               'parahippocampal gyrus',
               'superior temporal gyrus'
  ),
  study = c('MSBB','rnaSeqReprocessing'), 
  consortium = 'AMP-AD',
  assay = 'rnaSeq'
)

#Github Code Pull
thisRepo <- githubr::getRepo(
  repository = "Sage-Bionetworks/ampad-rnaseq-reprocessing",
  ref="branch",
  refName='main'
)
thisFile <- githubr::getPermlink(
  repository = thisRepo,
  repositoryPath=paste0('code/metadata_preprocessing/',
                        'msbb_preprocessing.R'
  )
)

activityName = 'Full MSBB Metadata'
activityDescription = 'Codified and Recoded Metadata'

write.csv(total_metadata,
          file = 'Full_MSBB_RNASeq_Covariates.csv',
          row.names = F,
          quote = F
)
ENRICH_OBJ <- synapser::synStore( synapser::File( 
  path='Full_MSBB_RNASeq_Covariates.csv',
  name = 'MSBB Full Covariates',
  parentId=activity$properties$id ),
  used = synids_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
file.remove("Full_MSBB_RNASeq_Covariates.csv")

## Upload cleaned metadata file and SageSeqr input version to synapse - internal Sage Location:
activityDescription = 'Cleaned Codified and Recoded Metadata'

write.csv(metadata,
          file = 'Cleaned_MSBB_RNASeq_Covariates.csv',
          row.names = F,
          quote = F
)
ENRICH_OBJ <- synapser::synStore( synapser::File( 
  path='Cleaned_MSBB_RNASeq_Covariates.csv',
  name = 'MSBB Cleaned Covariates',
  parentId=activity$properties$id ),
  used = synids_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
file.remove("Cleaned_MSBB_RNASeq_Covariates.csv")

## Upload Sageseqr file and SageSeqr input version to synapse - internal Sage Location:
activityDescription = 'Sageseqr Input Metadata'

write.csv(sageseqr_uncensored,
          file = 'Sageseqr_MSBB_RNASeq_Uncensored_Covariates.csv',
          row.names = F,
          quote = F
)

ENRICH_OBJ <- synapser::synStore( synapser::File( 
  path='Sageseqr_MSBB_RNASeq_Uncensored_Covariates.csv',
  name = 'MSBB Sageseqr Covariates',
  parentId=activity$properties$id ),
  used = synids_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
file.remove("Sageseqr_MSBB_RNASeq_Uncensored_Covariates.csv")


## Upload Sageseqr file to synapse
sageseqrmeta_parentid <- 'syn25808143'
activity <- synapser::synGetEntity(sageseqrmeta_parentid)
activityName = 'Sageseqr Input Metadata'
activityDescription = 'Cleaned Codified and Recoded Metadata'

write.csv(sageseqr_censored,
          file = 'Sageseqr_MSBB_RNASeq_Covariates_Censored.csv',
          row.names = F,
          quote = F
)

ENRICH_OBJ <- synapser::synStore( synapser::File( 
  path='Sageseqr_MSBB_RNASeq_Covariates_Censored.csv',
  name = 'MSBB Ages Censored Sageseqr Input Covariates',
  parentId=activity$properties$id ),
  used = synids_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
file.remove("Sageseqr_MSBB_RNASeq_Covariates_Censored.csv")

## Upload Sageseqr counts
counts_used <- c('syn21544664')
counts_parentid <- 'syn25808173'
activity <- synapser::synGetEntity(counts_parentid)
activityName = 'Sageseqr Input Counts'
activityDescription = 'MSBB RNASeq Counts for input to SageSeqr'

all.annotations.expression = list(
  dataType = 'geneExpression',
  dataSubtype = 'processed',
  resourceType = 'experimentalData',
  assay = 'rnaSeq',
  nucleicAcidSource = "bulk cell",
  isModelSystem = 'FALSE',
  isConsortiumAnalysis = 'TRUE',
  isMultiSpecimen = 'TRUE',
  fileFormat = 'txt',
  grant = 'U01AG046152',
  species = 'Human',
  organ = 'brain',
  tissue = c('frontal pole ',
             'inferior frontal gyrus',
             'parahippocampal gyrus',
             'superior temporal gyrus'
  ),
  study = c('MSBB','rnaSeqReprocessing'), 
  consortium = 'AMP-AD'
)

counts_write <- count
counts_write$feature <- row.names(counts_write)
counts_write <- counts_write[,c('feature',
                                colnames(counts_write)[!(colnames(counts_write) %in% 'feature')]
)
]
write.table(counts_write,
            file = 'MSBB_counts.txt',
            row.names = F,
            col.names = T,
            quote = F,
            sep = '\t'
)

ENRICH_OBJ <- synapser::synStore( synapser::File( 
  path='MSBB_counts.txt',
  name = 'MSBB Sageseqr Input Counts',
  parentId=activity$properties$id ),
  used = counts_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations.expression)
file.remove("MSBB_counts.txt")




