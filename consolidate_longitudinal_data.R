#Load required dependancies
rPackages <- c("RMySQL", "plyr", "dplyr", "magrittr", "GenomicRanges", "Biostrings", "igraph", "argparse", "sonicLength", "devtools") 
stopifnot(all(sapply(rPackages, require, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)))

#Set up and gather commandline arguments
setArguments <- function(){
  parser <- ArgumentParser(description = "Consolidate all longitudinal data from intSites database.")
  parser$add_argument("--sp_data", default = "hiv_specimen.database", help = "Group to use for specimen data.")
  parser$add_argument("--sites_data", default = "hiv_intsites.database", help = "Group to use for intsite data.")
  parser$add_argument("-s", action = 'store_true', help = "Abundance by sonicLength package (Berry, C. 2012).")
  parser$add_argument("-b", "--bp_corr", action = "store_true", help = "Correct breakpoint variation from pcr amplification.")
  parser$add_argument("-a", "--all_pats", action = "store_true", help = "Run analysis with all patients, not just longitudinal.")
  
  arguments <- parser$parse_args()
  arguments
}

arguments <- setArguments()

specimen.database <- arguments$sp_data
intsites.database <- arguments$sites_data
all.patients <- arguments$all_pats
breakpoint.correction <- arguments$bp_corr
sonic.method <- arguments$s
if(sonic.method){
  sonic.method <- "estAbund"
}else{
  sonic.method <- "fragLen"
}

#Establish connection with databases to query metadata and intsite data
junk <- sapply(dbListConnections(MySQL()), dbDisconnect) 
dbConn1 <- dbConnect(MySQL(), group = specimen.database)
dbConn2 <- dbConnect(MySQL(), group = intsites.database)
stopifnot(dbGetQuery(dbConn1, "SELECT 1")==1) 
stopifnot(dbGetQuery(dbConn2, "SELECT 1")==1) 

sp_query_select <- "SELECT parentAlias, trial, patient, cellType, timepoint, vcn FROM hivsp"
specimen_tb <- dbGetQuery(dbConn1, sp_query_select)
sam_query_select <- "SELECT parentAlias, childAlias, uniqRegion, primerType, linkerNum, bcNum, ngDNAInput FROM hivsam"
sample_tb <- dbGetQuery(dbConn1, sam_query_select)
cont_query_select <- "SELECT parentAlias, childAlias, uniqRegion, primerType, linkerNum, bcNum, ngDNAInput FROM contsam"
cont_tb <- dbGetQuery(dbConn1, cont_query_select)
int_query_select <- "SELECT * FROM samples"
intsites_samples_tb <- dbGetQuery(dbConn2, int_query_select)

sp_sam_tb <- merge(specimen_tb, sample_tb, by = "parentAlias", all = TRUE)
intsite_join_tb <- merge(sp_sam_tb, intsites_samples_tb, by.x = "childAlias", by.y = "sampleName")

if(!all.patients){
  long_pat_tb <- distinct(intsite_join_tb[, c("patient", "timepoint")]) %>%
    group_by(., patient) %>%
    summarize(., tps = length(timepoint))
  long_pats <- long_pat_tb[long_pat_tb$tps >= 2, "patient"]
  intsite_join_tb <- intsite_join_tb[intsite_join_tb$patient %in% long_pats$patient,]
}

sampleIDs <- sprintf("(%s)", paste(unique(intsite_join_tb$sampleID), collapse=",")) 

##Get integration site data from intsites.database
sql <- paste("SELECT DISTINCT * 
             FROM samples JOIN sites 
             ON samples.sampleID = sites.sampleID 
             JOIN pcrbreakpoints 
             ON pcrbreakpoints.siteID = sites.siteID  
             WHERE samples.sampleID in ", sampleIDs ) 
sites.uniq <- suppressWarnings( dbGetQuery(dbConn2, sql) )  
sites.uniq <- sites.uniq[, !duplicated(colnames(sites.uniq))] 

##get multihit sites 
sql <- paste("SELECT * 
             FROM samples JOIN multihitpositions 
             ON samples.sampleID = multihitpositions.sampleID 
             JOIN multihitlengths 
             ON multihitpositions.multihitID = multihitlengths.multihitID ", 
             "WHERE samples.sampleID in ",  sampleIDs ) 
sites.multi <- suppressWarnings( dbGetQuery(dbConn2, sql) )  
sites.multi <- sites.multi[, !duplicated(colnames(sites.multi))] 
sites.multi$breakpoint <- ifelse(sites.multi$strand=="+", 
                                 sites.multi$position+sites.multi$length-1, 
                                 sites.multi$position-sites.multi$length+1)

junk <- sapply(dbListConnections(MySQL()), dbDisconnect) 

#Analyze integration sites
source_url("https://raw.githubusercontent.com/cnobles/cloneTracker/master/cloneTracker.SOURCE_ME.R")
unique.gr <- db_to_granges(sites.uniq, keep.additional.columns = TRUE)
unique.gr$patient <- intsite_join_tb[match(unique.gr$sampleName, intsite_join_tb$childAlias), "patient"]
unique.gr$specimen <- paste0(unique.gr$specimen, "-",
                         sapply(strsplit(unique.gr$sampleName, split = "-"), "[[", 2))
unique.gr$timepoint <- intsite_join_tb[match(unique.gr$sampleName, intsite_join_tb$childAlias), "timepoint"]
multihit.gr <- db_to_granges(sites.multi, keep.additional.columns = TRUE)
multihit.gr$patient <- intsite_join_tb[match(multihit.gr$sampleName, intsite_join_tb$childAlias), "patient"]
multihit.gr$specimen <- paste0(multihit.gr$specimen, "-",
                         sapply(strsplit(multihit.gr$sampleName, split = "-"), "[[", 2))
multihit.gr$timepoint <- intsite_join_tb[match(multihit.gr$sampleName, intsite_join_tb$childAlias), "timepoint"]

#Unique Site Tracking
unique.grl <- split(unique.gr, unique.gr$patient)
unique.std <- unlist(GRangesList(lapply(unique.grl, function(sites){
  standardize_intsites(sites, standardize_breakpoints = breakpoint.correction)
})))
!!!unique.cond <- unlist(GRangesList(lapply(split(unique.std, unique.std$specimen), 
                                         function(sites){
  cond.sites <- condense_intsites(sites, return.abundance = TRUE, 
                                  method = sonic.method, replicates = "sampleName")
  mcols(cond.sites)[, c("sampleName", "sampleID", "siteID", "count")] <- NULL
  cond.sites
})))

patient.list <- split(unique.cond, unique.cond$patient)
longitudinal.list <- lapply(patient.list, function(gr) split(gr, gr$timepoint))
longitudinal.list <- longitudinal.list[sapply(longitudinal.list, length) >= 2]
long.track <- lapply(longitudinal.list, function(grl) track_clones(grl))




