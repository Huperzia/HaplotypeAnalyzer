#################
# Load packages #
#################
suppressMessages({
  library(shiny)
  library(bs4Dash)
  library(waiter)
  library(stringi)
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(ggtreeExtra)
  library(tidyr)
  library(dplyr)
  library(ggnewscale)
  library(scales)
  library(stringr)
  library(openxlsx)
  library(reshape2)
  library(cowplot)
  library(log4r)
})


#####################
# Initialize Logger #
#####################

log_file <- paste0("logs/", format(Sys.time(), "%Y-%m-%d_"), "serverlog.txt")
logger <- logger("INFO", appenders = file_appender(log_file, append = F))
info(logger, "Welcome to PopGenBrowser")
logger <- logger("INFO", appenders = file_appender(log_file, append = T))

info(logger, paste("Server Initialization", date()))


################
# Load Modules #
################

info(logger, "Loading Modules")

source("mods/aesthetics.R")
source("mods/database_config.R")

info(logger, "Finished Loading Modules")


########################
# Load Modules Servers #
########################

reac <- reactiveValues()


# interpret Locus input
interpretLocus <- function(range = NA){
  tmp <- unlist(str_split(range, pattern = ":"))
  chr <- tmp[1]
  tmp <- as.numeric(unlist(str_split(tmp[2], pattern = "[.]")))
  RANGE <- as.numeric(format(c(tmp[1], tmp[3]), scientific=F))
  
  return(c(chr, RANGE))
}


getBCF <- function(bed = NA, VCF = NA, VCFS = NA){
  
  FILE <- paste(bed[1], ".", bed[2], "..", bed[3], ".", VCF, ".vcf.gz", sep = "")
  FILES <- list.files("tmp/", pattern = "012")
  
  info(logger, paste("Getting BCF", FILE))    
  info(logger, paste("Getting BCF", FILES))
  
  if(!(paste(FILE, ".012", sep ="") %in% FILES)){
    system.time(system(paste("bcftools view ", VCFS,
                             " -r ", bed[1], ":", bed[2], "-", bed[3],
                             #                                                 " --threads 8",
                             " -Oz | ",
                             "vcftools --gzvcf -",
                             " --012",
                             " --out tmp/", FILE,
                             sep = "")))
  }
  
  tmp <- t(read.table(paste("tmp/", FILE, ".012", sep =""), header=F)[,-1])
  colnames(tmp) <- read.table(paste("tmp/", FILE, ".012.indv", sep =""))[,1]
  rownames(tmp) <- read.table(paste("tmp/", FILE, ".012.pos", sep =""))[,2]
  
  info(logger, paste("Got BCF", FILE))   
  
  return(list(tmp))
  
}


#  observeEvent(list(input$range), {

# Get counts
#    CNTS <- getcnts(bcfs = BCFS, key = key)

# Get pcas
#    tmp <- getpcas(bcfs = BCFS, cnts = CNTS)
#    pcas <- list(tmp[[1]][[1]], tmp[[2]][[1]])
#    CNTS <- list(tmp[[1]][[2]], tmp[[2]][[2]])
#    pcas <- list(tmp[[1]][[1]])
#    CNTS <- list(tmp[[1]][[2]])

# Get dapcs
#    tmp <- getdapcs(bcfs = BCFS, cnts = CNTS)
#    dapcs <- list(tmp[[1]][[1]], tmp[[2]][[1]])
#    CNTS <- list(tmp[[1]][[2]], tmp[[2]][[2]])

#    cnts <- CNTS[[1]]
#    df <- BCFS[[1]]

#    population_overview_samples <- reactive({
#      isolate({

#  })



getHeatMap <- function(df = NULL, trait = NULL){
  df2 <- df %>% select(vcfID, clade, clade.name,
                       clade.color, contains(trait))
  
  traits <- colnames(df2)[-c(1:4)]
  
  df3 <- df2 %>% 
    replace(is.na(.), 0) %>% 
    pivot_longer(traits, names_to = "trait", values_to = "value")
  
  ggplot(df3, aes(x = trait,  y = vcfID)) +
    geom_tile(aes(fill = value)) + 
    V1 + scale_fill_viridis_c()
}



getBarPlots <- function(df = NULL, trait = NULL){
  df2 <- df %>% select(vcfID, clade, clade.name,
                       clade.color, contains(trait))
  
  traits <- colnames(df2)[-c(1:4)]
  
  df3 <- df2 %>% 
    replace(is.na(.), 0) %>% 
    pivot_longer(traits, names_to = "trait", values_to = "value")
  
  df4 <- df3
  df4$conCol <- ifelse(df4$value < 0, "negative","positive")
  
  ggplot(df4, aes(x = vcfID, y = value, fill = conCol)) +
    geom_bar(stat = "identity") + 
    facet_wrap(~trait, nrow = 1) +
    scale_fill_manual(values = drac[c(4,1)], guide = "none") +
    scale_alpha(range = c(0.3, 1), guide = F) +
    V1 + scale_fill_manual(values=c(positive = drac[4], negative = drac[1])) +
    theme(strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.x = element_blank()) +
    coord_flip() 
}









