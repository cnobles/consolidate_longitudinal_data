#Load required dependancies
rPackages <- c("RMySQL", "plyr", "dplyr", "magrittr", "GenomicRanges", "hiAnnotator",
               "Biostrings", "igraph", "argparse", "sonicLength", "devtools") 
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
freeze <- unique(intsite_join_tb$refGenome)

annotate_sites <- function(sites, ref_seqs){
  sites <- doAnnotation(annotType = "within", sites, ref_seqs, colnam = "WithinGene",
                        feature.colnam = "name2", asBool = FALSE)
  sites <- doAnnotation(annotType = "nearest", sites, ref_seqs, colnam = "NearestGeneStart",
                        side = "5p", feature.colnam = "name2")
  sites <- doAnnotation(annotType = "nearest", sites, ref_seqs, colnam = "NearestGene",
                        feature.colnam = "name2")
  sites
}

makeUCSCsession(freeze=freeze)
ref.genes <- getUCSCtable("refGene", "RefSeq Genes", bsession=NULL, freeze=freeze)
genes.gr <- makeGRanges(ref.genes)

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

allSites <- unique.std
mcols(allSites) <- mcols(allSites)[,c("sampleName", "specimen", "patient", "timepoint", "count", "refGenome")]
write.table(as.data.frame(allSites, row.names = NULL), file = "allSites.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)

unique.cond <- condense_intsites(unique.std, grouping = "specimen", 
                                 return.abundance = TRUE, method = sonic.method,
                                 replicates = "sampleName")
mcols(unique.cond) <- mcols(unique.cond)[, c("specimen", "refGenome", "patient",
                                             "timepoint", "estAbund", "relAbund",
                                             "relRank", "posID")]
finalSites <- as.data.frame(
  annotate_sites(unique.cond, genes.gr), 
  row.names = NULL)
finalSites$position <- finalSites$start
finalSites <- finalSites[, c("seqnames", "strand", "position", "posID", "specimen", 
                             "patient", "timepoint", "estAbund", "relAbund", 
                             "relRank", "WithinGene", "WithinGeneOrt", "X5pNearestGeneStartDist", 
                             "X5pNearestGeneStart", "X5pNearestGeneStartOrt", "NearestGeneDist", 
                             "NearestGene", "NearestGeneOrt", "refGenome")]
write.table(finalSites, file = "finalSites.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)

patient.list <- split(unique.cond, unique.cond$patient)
longitudinal.list <- lapply(patient.list, function(gr) split(gr, gr$timepoint))
longitudinal.list <- longitudinal.list[sapply(longitudinal.list, length) >= 2]
longitudinal.sites <- lapply(longitudinal.list, track_clones, track.origin = FALSE)

contam.sites <- track_clones(patient.list, track.origin = FALSE)
crossover.sites <- contam.sites[names(contam.sites) %in% unlist(sapply(longitudinal.sites, names))]
distinct.sites <- lapply(longitudinal.sites, function(sites){
  sites[!names(sites) %in% names(contam.sites)]
})

intSites <- list(distinct.sites, crossover.sites, contam.sites)
names(intSites) <- c("distinct.sites", "crossover.sites", "contam.sites")

save(intSites, file = "intSites.RData")

longitudinalSites <- as.data.frame(
  annotate_sites(unlist(GRangesList(lapply(distinct.sites, unlist))), genes.gr), 
  row.names = NULL)
longitudinalSites$position <- longitudinalSites$start
longitudinalSites <- longitudinalSites[, c("seqnames", "strand", "position", "posID", "specimen", 
                                           "patient", "timepoint", "estAbund", "relAbund", 
                                           "relRank", "WithinGene", "WithinGeneOrt", "X5pNearestGeneStartDist", 
                                           "X5pNearestGeneStart", "X5pNearestGeneStartOrt", "NearestGeneDist", 
                                           "NearestGene", "NearestGeneOrt", "refGenome")]
write.table(longitudinalSites, file = "longitudinalSites.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

crossoverSites <- as.data.frame(
  annotate_sites(unlist(contam.sites), genes.gr), 
  row.names = NULL)
crossoverSites$position <- crossoverSites$start
crossoverSites <- crossoverSites[, c("seqnames", "strand", "position", "posID", "specimen", 
                                     "patient", "timepoint", "estAbund", "relAbund", 
                                     "relRank", "WithinGene", "WithinGeneOrt", "X5pNearestGeneStartDist", 
                                     "X5pNearestGeneStart", "X5pNearestGeneStartOrt", "NearestGeneDist", 
                                     "NearestGene", "NearestGeneOrt", "refGenome")]
write.table(crossoverSites, file = "crossoverSites.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(unique(crossoverSites$posID), file = "crossoverPosIDs.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#Analyze breakpoint sharing of fragments
analyze_bps <- function(grl_of_sites){ #grl needs to be split by posid
  bind_rows(lapply(grl_of_sites, function(sites){
    grl <- split(sites, width(sites))
    total.bp <- length(grl)
    patient.crossover.count <- length(grep(TRUE, sapply(grl, function(x) length(unique(x$patient)) > 1)))
    patient.crossover <- patient.crossover.count > 0
    replicate.crossover.count <- length(grep(TRUE, sapply(grl, function(x) length(unique(x$sampleName)) > 1)))
    replicate.crossover <- replicate.crossover.count > 0
    data.frame("site" = unique(generate_posID(sites)), "total.bps" = total.bp, 
               "patient.crossover" = patient.crossover, "patient.co.count" = patient.crossover.count,
               "replicate.crossover" = replicate.crossover, "replicate.co.count" = replicate.crossover.count)
  }))
}

patient.std.list <- split(unique.std, unique.std$patient)
longitudinal.std.list <- lapply(patient.list, function(gr) split(gr, gr$timepoint))
longitudinal.std.list <- longitudinal.std.list[sapply(longitudinal.std.list, length) >= 2]
longitudinal.std.sites <- lapply(longitudinal.std.list, track_clones, track.origin = FALSE)

contam.std.sites <- track_clones(patient.std.list, track.origin = FALSE)
crossover.std.sites <- contam.sites[names(contam.sites) %in% unlist(sapply(longitudinal.sites, names))]
distinct.std.sites <- lapply(longitudinal.sites, function(sites){
  sites[!names(sites) %in% names(contam.std.sites)]
})

contam.shared.bp <- analyze_bps(contam.std.sites)
crossover.shared.bp <- analyze_bps(crossover.std.sites)
distinct.shared.bp <- lapply(distinct.std.sites, analyze_bps)

bp.analysis <- list(distinct.shared.bp, crossover.shared.bp, contam.shared.bp)
names(bp.analysis) <- c("distinct.sites", "crossover.sites", "contam.sites")

save(bp.analysis, file = "bp.analysis.RData")

#Annotate sites
contam.sites.info <- unlist(GRangesList(lapply(contam.sites, function(sites){
  reduce(flank(sites, -1, start = TRUE), min.gapwidth = 0L)
})))
crossover.sites.info <- unlist(GRangesList(lapply(crossover.sites, function(sites){
  reduce(flank(sites, -1, start = TRUE), min.gapwidth = 0L)
})))
distinct.sites.info <- lapply(distinct.sites, function(grl){
  unlist(GRangesList(lapply(grl, function(sites) reduce(flank(sites, -1, start = TRUE), min.gapwidth = 0L))))
})
distinct.sites.info <- distinct.sites.info[sapply(distinct.sites.info, length) > 0]

contam.sites.info <- annotate_sites(contam.sites.info, genes.gr)
crossover.sites.info <- annotate_sites(crossover.sites.info, genes.gr)
distinct.sites.info <- lapply(distinct.sites.info, annotate_sites, ref_seqs = genes.gr)

intSites.info <- list(distinct.sites.info, crossover.sites.info, contam.sites.info)
names(intSites.info) <- c("distinct.sites", "crossover.sites", "contam.sites")

save(intSites.info, file = "intSites.info.RData")

siteInfo <- sort(unique(granges(unique.cond)))
siteInfo$posID <- generate_posID(siteInfo)
siteInfo <- annotate_sites(siteInfo, genes.gr)
siteInfo <- as.data.frame(siteInfo, row.names = NULL)
siteInfo$position <- siteInfo$start
siteInfo <- siteInfo[, c("seqnames", "strand", "position", "posID", "WithinGene", 
                           "WithinGeneOrt", "X5pNearestGeneStartDist", "X5pNearestGeneStart", 
                           "X5pNearestGeneStartOrt", "NearestGeneDist", "NearestGene", "NearestGeneOrt")]
write.table(siteInfo, file = "siteInfo.tsv", quote = FALSE, sep = "\t",
            row.names = FALSE)

#site.info <- lapply(1:length(intSites.info), function(i){
#  sites <- intSites.info[[i]]
#  if(class(sites) == "list"){
#    sites.list <- lapply(1:length(sites), function(j){
#      patient <- names(sites[j])
#      patient.sites <- sites[[j]]
#      mcols(patient.sites)$Patient <- patient
#      patient.sites
#    })
#    sites <- do.call(c, lapply(1:length(sites.list), function(i) sites.list[[i]]))
#  }else{
#    mcols(sites)$Patient <- "NA"
#  }
#  typeTable <- data.frame(row.names = names(intSites.info),
#                          "siteType" = c("distinct", "crossover", "contaminating"))
#  sites$siteType <- typeTable[names(intSites.info[i]), "siteType"]
#  sites
#})
#site.info <- do.call(c, lapply(1:length(site.info), function(i) site.info[[i]]))
#site.info$posid <- names(site.info)
#names(site.info) <- NULL
#mcols(site.info) <- mcols(site.info)[,c("siteType", "posid", "Patient", "WithinGene", "WithinGeneOrt", 
#  "X5pNearestGeneStartDist", "X5pNearestGeneStart", "X5pNearestGeneStartOrt", 
#  "NearestGeneDist", "NearestGene", "NearestGeneOrt")]