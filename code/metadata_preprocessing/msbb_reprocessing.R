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
colnames(clinical)[ colnames(clinical) == 'Braak'] <- 'braaksc'

clinical$diagnosis <- 'OTHER'
clinical[ clinical$braaksc >= 4 &
            clinical$ceradsc <= 2 &
            clinical$CDR >= 1  , ]$diagnosis <- 'AD'
clinical[ clinical$braaksc <= 3 &
            clinical$ceradsc >= 3 &
            clinical$CDR <= 0.5  , ]$diagnosis <- 'CT'

# apoe Allele:
used_synids <- c(used_synids,biospecimin,'syn26452711')

apoe_fix <- read.csv(synapser::synGet('syn26452711')$path)
row.names(apoe_fix) <- apoe_fix$individualID
apoe_fix <- apoe_fix[apoe_fix$Exclude == FALSE,]
clinical$apoe4_allele <- NA
apoe_fix$final_apoe <- NA
apoe_fix$final_apoe <- apoe_fix$WGS_e4_dosage
apoe_fix[is.na(apoe_fix$WGS_e4_dosage),]$final_apoe <- apoe_fix[is.na(apoe_fix$WGS_e4_dosage),]$WES_e4_dosage

for( i in 1:dim(clinical)[1]){
  if(clinical[i,]$individualID %in% row.names(apoe_fix) ){
    clinical[i,]$apoe4_allele <- apoe_fix[clinical[i,]$individualID,]$final_apoe
  }
}

# conver pmi time interval from minutes to hours
clinical$pmi <- clinical$pmi/60

# Get biospecimin metadata
region_aide <- read.csv(synapser::synGet('syn6100548')$path, header = T)
region_aide <- region_aide[ region_aide$fileType == 'bam',]
row.names(region_aide) <- region_aide$sampleIdentifier

biospecimin <- 'syn21893059'
used_synids <- c(used_synids,biospecimin,'syn6100548')
biospecimin <- synGet(biospecimin)$path %>%
  read.csv(header = T, stringsAsFactors = F)
biospecimin <- biospecimin[ biospecimin$assay=='rnaSeq',]
biospecimin<-biospecimin[biospecimin$specimenID %in% colnames(counts),]

biospecimin[ biospecimin$specimenID %in% row.names(region_aide),]$BrodmannArea <-
  as.numeric(
    gsub('BM','',region_aide[ biospecimin[ biospecimin$specimenID %in% row.names(region_aide),]$specimenID,]$BrodmannArea)
  )

biospecimin$tissue[biospecimin$BrodmannArea==22] <- 'STG'
biospecimin$tissue[biospecimin$BrodmannArea==10] <- 'FP'
biospecimin$tissue[biospecimin$BrodmannArea==44] <- 'IFG'
biospecimin$tissue[biospecimin$BrodmannArea==36] <- 'PHG'

#Toss sample swaps
biospecimin <- biospecimin[ biospecimin$exclude == 'FALSE',]
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


sageseqr_uncensored <- sageseqr_uncensored[,
                                           colnames(sageseqr_uncensored)[
                                             !(colnames(sageseqr_uncensored) %in%
                                                 c('braaksc', 'CDR', 'plaqueMean', 'ceradsc', 'apoeGenotype',
                                                   'Tissue.APOE4','individualID', 'BrodmannArea','barcode'))
                                           ]]
sageseqr_uncensored <- sageseqr_uncensored[ !is.na(sageseqr_uncensored$tissue), ]

sex_swapped <- c( 'hB_RNA_12901', 'hB_RNA_12934', 'BM_36_296', 'hB_RNA_10622',
                  'hB_RNA_10622_L43C014', 'hB_RNA_10702', 'hB_RNA_10702_E007C014', 'hB_RNA_7765', 'hB_RNA_7765_resequenced',
                  'hB_RNA_7995', 'hB_RNA_8015', 'hB_RNA_8025', 'hB_RNA_8025_E007C014', 'hB_RNA_8385', 'hB_RNA_8305')
expression_outliers <- c('hB_RNA_10577', 'hB_RNA_12964', 'hB_RNA_13144', 'hB_RNA_13397', 'hB_RNA_10567', 'hB_RNA_10577', 'hB_RNA_7855', 'hB_RNA_9005')
total_toss <- c(sex_swapped,expression_outliers)

################################################################################
## --   Push to Synapse   --  ##
#Upload full file and SageSeqr input version to synapse - internal Sage Location:
internal_parentid <- 'syn25872211'
folder_loc <- 'syn25872212'

#Set Activity
activity <- synapser::synGet(folder_loc)

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
                        'msbb_reprocessing.R'
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
  used = used_synids,
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
  used = used_synids,
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
  name = 'MSBB  Ages Uncensored Sageseqr Input Covariates',
  parentId=activity$properties$id ),
  used = used_synids,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
file.remove("Sageseqr_MSBB_RNASeq_Uncensored_Covariates.csv")


## Upload Sageseqr file to synapse
sageseqrmeta_parentid <- 'syn25808143'
activity <- synapser::synGet(sageseqrmeta_parentid)
activityName = 'Sageseqr Input Metadata'
activityDescription = 'Cleaned Codified and Recoded Metadata'


children <- synapser::synGetChildren(activity$properties$id)$asList()
cesoredmeta_used <- NULL
for(i in 1:length(children)) {
  if(children[[i]]$name == "MSBB  Ages Uncensored Sageseqr Input Covariates") {
    cesoredmeta_used <- children[[i]]$id
  }
}
sageseqr_censored <- sageseqr_censored[!is.na(sageseqr_censored$individualID), ]
write.csv(sageseqr_censored,
          file = 'Sageseqr_MSBB_RNASeq_Covariates_Censored.csv',
          row.names = F,
          quote = F
)

ENRICH_OBJ <- synapser::synStore( synapser::File(
  path='Sageseqr_MSBB_RNASeq_Covariates_Censored.csv',
  name = 'MSBB Ages Censored Sageseqr Input Covariates',
  parentId=activity$properties$id ),
  used = cesoredmeta_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
file.remove("Sageseqr_MSBB_RNASeq_Covariates_Censored.csv")

#Store Tissue Specific Files:
for( tis in names(table(sageseqr_censored$tissue))){

  meta_write <- sageseqr_uncensored[as.character(sageseqr_uncensored$tissue) == tis,]
  meta_write <- meta_write[ !(as.character(meta_write$specimenID) %in% total_toss), ]
  meta_write <- meta_write[ !is.na(meta_write$tissue),]
  meta_write <- meta_write[complete.cases(meta_write),]
  write.csv(meta_write[ , colnames(meta_write)[!(colnames(meta_write) %in% 'tissue')]],
            file = paste0('Sageseqr_MSBB_', tis,'_RNASeq_Covariates.csv'),
            row.names = F,
            quote = F
  )

  ENRICH_OBJ <- synapser::synStore( synapser::File(
    path=paste0('Sageseqr_MSBB_', tis,'_RNASeq_Covariates.csv'),
    name = paste0('MSBB ', tis,' Sageseqr Input Covariates'),
    parentId=activity$properties$id ),
    used = used_synids,
    activityName = activityName,
    executed = thisFile,
    activityDescription = activityDescription
  )
  synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
  file.remove(paste0('Sageseqr_MSBB_', tis,'_RNASeq_Covariates.csv'))

}

## Upload Sageseqr counts
counts_used <- c('syn21544664')
counts_parentid <- 'syn25808173'
activity <- synapser::synGet(counts_parentid)
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

counts_write <- counts
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

for( tis in names(table(sageseqr_censored$tissue))){
  meta_write <- sageseqr_uncensored[as.character(sageseqr_uncensored$tissue) == tis,]
  meta_write <- meta_write[complete.cases(meta_write),]
  meta_write <- meta_write[ !(as.character(meta_write$specimenID) %in% total_toss), ]
  counts_wr <- counts_write[,c('feature',
                                  colnames(counts_write)[(colnames(counts_write) %in% meta_write$specimenID)]
  )
  ]
  write.table(counts_wr,
              file = paste0('MSBB_', tis, '_counts.txt'),
              row.names = F,
              col.names = T,
              quote = F,
              sep = '\t'
  )

  ENRICH_OBJ <- synapser::synStore( synapser::File(
    path=paste0('MSBB_', tis, '_counts.txt'),
    name = paste0('MSBB ', tis, ' Sageseqr Input Counts'),
    parentId=activity$properties$id ),
    used = counts_used,
    activityName = activityName,
    executed = thisFile,
    activityDescription = activityDescription
  )
  synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations.expression)
  file.remove(paste0('MSBB_', tis, '_counts.txt'))
}

