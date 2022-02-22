library(dplyr)
synapser::synLogin()
code_file_name <- 'agora_rnaseq_formating.R'
synids_used_exp <- NULL
synids_used_de <- NULL

download <- function( id ){
  
  foo <- data.table::fread(synapser::synGet(as.character(id))$path, sep = '\t', header = T) %>%
    data.frame()
  # Fix Direction
  foo$Direction <- toupper(foo$Direction)
  
  # Make Tissue
  foo$Tissue <- NA
  foo$Study <- NA
  if (all(grepl('_CBE', foo$Comparison) | grepl('_TCX',foo$Comparison))) {
    foo$Study <- 'MAYO'
    foo[grepl('_CBE', foo$Comparison),]$Tissue <- 'CBE'
    foo[grepl('_TCX', foo$Comparison),]$Tissue <- 'TCX'
  }else {
    if (all(grepl('_DLPFC', foo$Comparison) | grepl('_ACC',foo$Comparison) | grepl('_PCC',foo$Comparison))) {
      foo$Study <- 'ROSMAP'
      foo[grepl('_DLPFC', foo$Comparison),]$Tissue <- 'DLPFC'
      foo[grepl('_PCC', foo$Comparison),]$Tissue <- 'PCC'
      foo[grepl('_ACC', foo$Comparison),]$Tissue <- 'ACC'
    } else{
      if (all(grepl('_FP', foo$Comparison) | grepl('_STG',foo$Comparison) | grepl('_IFG',foo$Comparison) | grepl('_PHG',foo$Comparison))) {
        foo$Study <- 'MSSM'
        foo[grepl('_FP', foo$Comparison),]$Tissue <- 'FP'
        foo[grepl('_STG', foo$Comparison),]$Tissue <- 'STG'
        foo[grepl('_IFG', foo$Comparison),]$Tissue <- 'IFG'
        foo[grepl('_PHG', foo$Comparison),]$Tissue <- 'PHG'
      } else {
        warning( message('Data does not fit the constraints of AMP-AD tissue nomenclature'))
      }
    }
  }
  
  if (all(grepl('age_death', foo$Comparison))) {
    foo$Model <- 'Diagnosis.AOD'
    foo$Comparison <- 'AD-CONTROL'
    foo$Sex <- 'ALL'
  } else {
    if (all(grepl('_female_', foo$Comparison) | grepl('_male_', foo$Comparison))) {
      foo$Sex <- NA
      foo[ grepl('_female_', foo$Comparison),]$Sex <- 'FEMALE'
      foo[ grepl('_male_',foo$Comparison),]$Sex <- 'MALE'
      foo$Model <- 'Diagnosis.Sex'
      foo$Comparison <- 'AD-CONTROL'
        
      #comparison <- ''
    } else {
      contrasts <- paste0( 'AD_', names(table(foo$Tissue)), ' - CT_', names(table(foo$Tissue)) )
      if (all(foo$Comparison %in% contrasts)) {
        foo$Model <- 'Diagnosis'
        foo$Comparison <- 'AD-CONTROL'
        foo$Sex <- 'ALL'
      } else {
        warning( message('Data does not fit the constraints of AMP-AD tissue contrast schema'))
      }
    }
  }
  colnames(foo)[ colnames(foo) == 'percentage_gene_gc_content'] <- 'percentage_gc_content'
  colnames(foo)[ colnames(foo) == 'gene_length'] <- 'gene.length'
  
  
  foo <- foo[, c('Model', 'Tissue', 'Comparison', 'ensembl_gene_id', 'logFC', 'CI.L', 'CI.R', 'AveExpr', 't', 'P.Value', 
          'adj.P.Val', 'gene_biotype', 'chromosome_name', 'Direction', 'hgnc_symbol', 'percentage_gc_content', 'gene.length', 'Sex', 'Study')
     ]
 
  return(foo)
}

###### INGEST DE Data:
# MAYO DE
  #daignosis <- 'syn27024969'
  #aod <- 'syn27024972'
  #sex <- 'syn27024970'
mayo_de <- as.data.frame(do.call(rbind,lapply(c('syn27024969','syn27024972', 'syn27024970' ),download)))


# MSSM DE
  #daignosis <- 'syn27068762'
  #aod <- 'syn27068764'
  #sex <- 'syn27068763'
mssm_de <-  as.data.frame(do.call(rbind,lapply(c('syn27068762','syn27068764','syn27068763'),download)))

# ROSMAP DE
  #daignosis <- 'syn26967457'
  #aod <- 'syn26967459'
  #sex <- 'syn26967458'
rosmap_de <- as.data.frame(do.call(rbind,lapply(c('syn26967457','syn26967459','syn26967458'),download)))

de_stats <- as.data.frame(rbind( mayo_de,mssm_de,rosmap_de ))
synids_used_de <- c('syn27024969', 'syn27024972',  'syn27024970',
                    'syn27068762','syn27068764','syn27068763',
                    'syn26967457','syn26967459','syn26967458')
########################
#### Expression Statistics
# Pull Metadata
mayo_md <- data.table::fread(synapser::synGet('syn27024950')$path, header = T, sep = '\t') %>%
  tibble::column_to_rownames(var = 'specimenID') %>%
  as.data.frame()
rosmap_md <- data.table::fread(synapser::synGet('syn26967450')$path, header = T, sep = '\t') %>%
  tibble::column_to_rownames(var = 'specimenID') %>%
  as.data.frame()
mssm_md <- data.table::fread(synapser::synGet('syn27068753')$path, header = T, sep = '\t') %>%
  tibble::column_to_rownames(var = 'specimenID') %>%
  as.data.frame()

## Pull Expression Stats
mayo_exp <- data.table::fread(synapser::synGet('syn27024965')$path, header = T, sep = '\t') %>%
  tibble::column_to_rownames(var = 'feature') %>%
  as.data.frame()
rosmap_exp <- data.table::fread(synapser::synGet('syn26967453')$path, header = T, sep = '\t') %>%
  tibble::column_to_rownames(var = 'feature') %>%
  as.data.frame()
mssm_exp <- data.table::fread(synapser::synGet('syn27068756')$path, header = T, sep = '\t') %>%
  tibble::column_to_rownames(var = 'feature') %>%
  as.data.frame()

synids_used_exp <- c('syn27024950', 'syn27024965', 
                    'syn26967450','syn26967453',
                    'syn27068753','syn27068756')

median_finder <- function( tiss, expression, md ){
  temp <- expression[, row.names( md[ md$tissue ==tiss, ] )] 
  apply( temp, 1, summary ) %>%
    t() %>%
    as.data.frame() %>%
    mutate(tissue = tiss) %>%
    rename(minimumLogCPM = Min.) %>%
    rename(maximumLogCPM = Max.) %>%
    rename(medianLogCPM = Median) %>%
    rename(meanLogCPM = Mean) %>%
    rename(quartile1LogCPM = `1st Qu.`) %>%
    rename(quartile3LogCPM = `3rd Qu.`) %>%
    tibble::rownames_to_column(var = 'ensembl_gene_id') %>%
    return()
}

mayo_stats <- as.data.frame(do.call( 
  rbind, 
  lapply( 
    names(table(mayo_md$tissue)), 
    median_finder,
    expression = mayo_exp, 
    md = mayo_md
  )
))
rosmap_stats <- as.data.frame(do.call( 
  rbind, 
  lapply( 
    names(table(rosmap_md$tissue)), 
    median_finder,
    expression = rosmap_exp, 
    md = rosmap_md
  )
))
mssm_stats <- as.data.frame(do.call( 
  rbind, 
  lapply( 
    names(table(mssm_md$tissue)), 
    median_finder,
    expression = mssm_exp, 
    md = mssm_md
  )
))

exp_stats <- as.data.frame(rbind(mayo_stats, rosmap_stats, mssm_stats))

########################################
## write and push de_stats and exp_stats

# De destionation -> syn14237651 - parentid = syn17015330, differentialExpressionSummary.tsv 
# exp destination -> syn12514804 - parentid = syn7525089, amp_ad_median_expression.csv
de_parentid <- 'syn17015330'
exp_parentid <- 'syn7525089'

de_annot <- synapser::synGet('syn14237651')$annotations
exp_annot <- de_annot
exp_annot$analysisType <- "summary expression"

#Githubr Code Pull
thisRepo <- githubr::getRepo(
  repository = "Sage-Bionetworks/ampad-rnaseq-reprocessing",
  ref="branch",
  refName='main'
)
thisFile <- githubr::getPermlink(
  repository = thisRepo,
  repositoryPath=paste0('code/',
                        code_file_name
  )
)

#push files to synapse:

# Expression Stats
write.csv(exp_stats, 'amp_ad_median_expression.csv',
            row.names = F,
            quote = F
)

activityName = 'Amp-AD Expression Stats'
activityDescription = 'Expression Summary Statistics for Amp-AD Transcriptomics Data'

ENRICH_OBJ <- synapser::synStore( synapser::File(
  path='amp_ad_median_expression.csv',
  name = 'Amp-AD Expression Stats',
  parentId=exp_parentid ),
  used = synids_used_exp,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = exp_annot)
file.remove("amp_ad_median_expression.csv")

# DE Stats
write.table(de_stats, 'differentialExpressionSummary.tsv',
            row.names = F,
            col.names = T,
            sep='\t',
            quote = F
)

activityName = 'All Differential Expression (Merged)'
activityDescription = 'Merged Differential Expression Statistics fro Amp-AD Transcriptomics Data'

ENRICH_OBJ <- synapser::synStore( synapser::File(
  path='differentialExpressionSummary.tsv',
  name = 'Amp-AD Differential Expression Stats',
  parentId=de_parentid ),
  used = synids_used_de,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = de_annot)
file.remove("differentialExpressionSummary.tsv")
