#############
# Databases #
#############

#~~~ 
#~ Population Metadata
#~~~

genotypeDB <- read.csv("data/populationKeyWithMetadata.csv", header = T)
key <- read.csv("data/genotypesKey.csv", stringsAsFactors = F)
meta <- read.csv("data/db_sweetcapMetadata_v05.csv")
meta[,49:ncol(meta)] <- scale(meta[,49:ncol(meta)])

#~~~ 
#~ Genomics
#~~~

bcfs_dir <- "/media/mink/purple/mink/1.research/1.databases/2.genomics/1.bcfs/"
VCFS <- paste0(bcfs_dir, "Ia453_sweetcap_v0.4_28M.bcf.gz")
VCFID <- c("28.5M")

#  siteCounts <- readRDS("data/countPlots.RDS")
#  output$counts.sites <-  renderPlot({siteCounts[[3]]})
#  output$counts.snps <-  renderPlot({siteCounts[[1]]})
#  output$counts.indels <-  renderPlot({siteCounts[[2]]})





