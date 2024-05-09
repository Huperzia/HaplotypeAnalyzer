#################
# Load packages #
#################
library(log4r)

log_file <- "logs/serverlog.txt"
logger <- logger("INFO", appenders = file_appender(log_file))
info(logger, "Welcome to PopGenBrowser")

info(logger, "Loading libraries")

library(tidyverse)

library(shiny)
library(bs4Dash)
library(waiter)
library(stringi)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

library(scales)
library(stringr)
library(openxlsx)
library(reshape2)

library(cowplot)
library(viridis)
info(logger, "Finished loading libraries")


#############
# Databases #
#############

info(logger, "Loading Databases")

genotypeDB <- read.csv("data/populationKeyWithMetadata.csv", header = T)
#~~~~~~~~~~~~

info(logger, "Finished loading databases")


###########
# Cluster #
###########

#no_cores <- 4
#cl <- makeCluster(no_cores)
#registerDoParallel(cl)

#~~~~~~~~~#


##############
# Aesthetics #
##############

info(logger, "Loading aesthetics")

drac <- c("#50fa7b", "#ffb86c", "#bd93f9", "#ff79c6", 
          "#ff5555", "#f1fa8c", "#6272a4", "#8be9fd", 
          "#282a36", "#44475a", "#44475a", "#f8f8f2")

cols <- drac

V1 <- theme_bw() + 
  theme(plot.background = element_rect(fill = drac[9], color = drac[9]),
        panel.grid = element_blank(),   
        panel.border = element_rect(color = drac[12]),
        panel.background = element_blank(),
        axis.ticks = element_line(color = drac[12]),
        strip.background = element_blank(),
        strip.text = element_text(colour = drac[12], size = 12, hjust = 0)) 


##############
# Functions #
##############

getBarPlots <- function(df = NULL, traits = NULL, levs = NULL){
  df2 <- df %>% select(vcfID, clade, clade.name,
                       clade.color, all_of(traits))
  
  df3 <- df2 %>% 
    replace(is.na(.), 0) %>% 
    pivot_longer(traits, names_to = "trait", values_to = "value")
  
  df4 <- df3
  df4$conCol <- ifelse(df4$value < 0, "negative","positive")
  
  df4$vcfID <- factor(df4$vcfID, levels = rev(levs))

  ggplot(df4, aes(x = vcfID, y = value, fill = conCol, alpha = abs(value))) +
    geom_bar(stat = "identity") + 
    facet_wrap(~trait, nrow = 1, scales = "free") +
    scale_alpha_continuous(range = c(0.75, 1), guide = "none") +
    scale_fill_manual(values=c(positive = drac[1], negative = drac[4]), guide = "none") +
    theme(strip.background = element_blank(),
          strip.text = element_text(colour = drac[12], size = 12, hjust = 0),
          plot.background = element_rect(fill = drac[9], color = NA),   
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank()) +
    coord_flip()
}


getHeatMap <- function(df = NULL, traits = NULL){
  df2 <- df %>% select(vcfID, clade, clade.name,
                       clade.color, traits)
  
  traits <- colnames(df2)[-c(1:4)]
  
  df3 <- df2 %>% 
    replace(is.na(.), 0) %>% 
    pivot_longer(traits, names_to = "trait", values_to = "value")
  
  ggplot(df3, aes(x = trait,  y = vcfID)) +
    geom_tile(aes(fill = value)) + 
    scale_fill_viridis_c()
}





##########
# Server #
##########
info(logger, "Initializing Server")

server <- function(input, output, session) {
  useAutoColor()
  # user menu
  
  output$user <- renderUser({
    dashboardUser(
      name = "Zea mays",
      image = "icons/icon.svg",
      title = "Zea mays",
      subtitle = "User"
    )
  })
  

  ####################
  ####################
  ## OBSERVE INPUTS ##
  ####################
  ####################
  # Start
  #  input <- data.frame(NA)
  #  input$range <- "Chr2:236682000..236687000"
  #  input$range <- "Chr3:223900000..223970000"  
  #  input$VCF <- "intersection"
  #  input$NAME <- "sh2"
	
info(logger, "Observing Inputs")  

  key <- read.csv("data/genotypesKey.csv", stringsAsFactors = F)
  VCFS <- c("/media/mink/purple/mink/1.research/1.databases/2.genomics/1.bcfs/Ia453_sweetcap_v0.4_28M.bcf.gz")
  
  VCFID <- c("Intersection")
  
#  siteCounts <- readRDS("data/countPlots.RDS")
  
#  output$counts.sites <-  renderPlot({siteCounts[[3]]})
#  output$counts.snps <-  renderPlot({siteCounts[[1]]})
#  output$counts.indels <-  renderPlot({siteCounts[[2]]})
  
reac <- reactiveValues(VCFS = VCFS)
  
  observeEvent(list(input$range), {
    
    # Get input range
    tmp <- unlist(str_split(input$range, pattern = ":"))
    chr <- tmp[1]
    tmp <- as.numeric(unlist(str_split(tmp[2], pattern = "[.]")))
    RANGE <- as.numeric(format(c(tmp[1], tmp[3]), scientific=F))
    
info(logger, "Retrieving inputs")

    # Get BCFs
    i = 1
    FILE <- paste(chr, ".", RANGE[1], "..", RANGE[2], ".", VCFID[i], ".vcf.gz", sep = "")

info(logger, paste("Getting BCF", VCFS))
info(logger, paste("Getting BCF", FILE))                      

    FILES <- list.files("tmp/", pattern = "012")

info(logger, paste("Getting BCF", FILES))

     if(!(paste(FILE, ".012", sep ="") %in% FILES)){
     system.time(system(paste("bcftools view ", VCFS[i],
                              " -r ", chr, ":", RANGE[1], "-", RANGE[2],
#                               " --threads 8",
                                " -Oz | ",
                                "vcftools --gzvcf -",
                                " --012",
                                " --out tmp/", FILE,
                                sep = "")))
                      }
			
info(logger, paste("Got BCF", FILE))                      

      tmp <- t(read.table(paste("tmp/", FILE, ".012", sep =""), header=F)[,-1])
      colnames(tmp) <- read.table(paste("tmp/", FILE, ".012.indv", sep =""))[,1]
      rownames(tmp) <- read.table(paste("tmp/", FILE, ".012.pos", sep =""))[,2]
                      
      BCFS <- list(tmp)


info(logger, "Retrieved bcfs")
    
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
        tst <- BCFS[[1]]
        print(dim(tst))
        tree <- as.phylo(hclust(dist(t(tst))))
        print(dim(tst))
        print(tree)
                      
			tst <- BCFS[[1]]
                      reac$tst <- tst
                      reac$tree <- as.phylo(hclust(dist(t(tst))))
                      reac$FILE <- FILE
                      reac$RANGE <- RANGE
		      reac$chr <- chr

  })
  
  
  
  observeEvent(list(input$ID_TYPE, reac$tst, input$populationA, input$populationB, input$trait_textinput), {
	 tst <- isolate(reac$tst)
    tree <- isolate(reac$tree) 
    FILE <- isolate(reac$FILE)
    RANGE <- isolate(reac$RANGE)
	chr <- isolate(reac$chr)

  cat(chr, ":", RANGE, "\n")  

  
  
    ###################
    # Process Samples #
    ###################
    #    input$caption <- "S210260\nS210261\nS210290\nS210291"
    #    input$ID_TYPE <- "vcf_id"
  
	popA <- stri_split_lines(input$populationA)[[1]]
	popB <- stri_split_lines(input$populationB)[[1]]
	TRAITS <- stri_split_lines(input$trait_textinput)[[1]]
	
	TRAITS <- TRAITS[TRAITS %in% colnames(meta)]
	
	print(TRAITS)
	print(popA)
	
	popU <- meta$vcfID[!(meta$vcfID %in% popA) & !(meta$vcfID %in% popB)]
	
	samp.pops <- list(popA = popA, popB = popB, popU = popU)
	
  cat("\n###\n\n###\n")

    #####
    
    
    #################
    # Found/Missing #
    #################
    
    tree2 <- groupOTU(tree, samp.pops)
 
    V2 <-   theme(axis.text = element_blank(), 
            panel.border = element_blank(), 
            panel.background = element_blank(),
            plot.background = element_rect(fill = drac[9], color = NA), 
            axis.ticks = element_blank())        
    
    p <- ggtree(tree2, layout = "rectangular", branch.length = "none", aes(color = group, size = group)) + 
       scale_color_manual(
        values = drac[c(5, 8, 12, 7)],
        guide = "none"
     )  + 
      scale_size_manual(
        values = c(1, 1, 0.2, 0.2),
        guide = "none"
      ) + theme_void() + V2
   
    tst3 <- data.frame(id = row.names(t(tst)), t(tst), check.names = F)
    tst3$id <- factor(tst3$id, levels = rev(get_taxa_name(p)))
    
    tst4 <- tst3 %>% 
      pivot_longer(cols = colnames(tst3)[-c(1)], names_to = "Pos")


    p5 <- getBarPlots(df = meta, traits = TRAITS, levs = get_taxa_name(p))
    
    p8 <-  ggplot(tst4, aes(x = Pos, y = id, fill = factor(value))) +  
      geom_tile() +
      scale_fill_manual(
        values=drac[c(3,10,2,1)],
        guide="none"
      ) + theme_bw() + V2 +
      theme(axis.title = element_blank())        
    
    
    
    p6 <- plot_grid(plotlist = list(p, p8, p5, p + scale_x_reverse()), 
                    ncol = 4, align = 'hv', axis = "tblr", 
                    rel_widths = c(0.1, 0.6, 0.2, 0.1))

    p9 <- p + geom_text(aes(label=label), hjust = 0, color = drac[12], size = 4)  + 
      theme_tree(bgcolor = drac[9], plot.margin = margin(6, 120, 6, 6)) +
      coord_cartesian(clip = 'off') + V2
    
    p10 <- plot_grid(plotlist = list(p9, p8, p5, p + scale_x_reverse()), 
                    ncol = 4, align = 'hv', axis = "tblr", 
                    rel_widths = c(0.2, 0.5, 0.2, 0.1))
    
    output$population_overview_samples4 <- renderPlot( p6 )
    output$population_overview_samples2 <- renderPlot( p10 )
    
    
    ord <- p[[1]] %>% 
      select(label, isTip, y) %>% 
      filter(isTip == T) %>% select(label, y) %>% 
      arrange(desc(y))

    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste(chr, "_", RANGE[1], "-", RANGE[2], "_", "28M", ".csv", sep = "")
      },
      content = function(file) {
        write.csv(tst3[ord$label,], file, row.names = F)
      }
    )


    output$downloadData2 <- downloadHandler(
      filename = function() {
        paste(chr, "_", RANGE[1], "-", RANGE[2], "_", "28M", ".xlsx", sep = "")
      },
      content = function(file) {
          xl.name <- paste(chr, "_", RANGE[1], "-", RANGE[2], "_", "28M", sep = "")
          wb <- createWorkbook()
          addWorksheet(wb, xl.name)
          writeData(wb, xl.name, tst3[ord$label,])
          
          conditionalFormatting(wb, xl.name,
                                cols = 2:ncol(tst3),
                                rows = 2:(nrow(tst3)+1), 
                                rule = "==-1", 
                                style = createStyle(bgFill = drac[3])
          )
          conditionalFormatting(wb, xl.name,
                                cols = 2:ncol(tst3),
                                rows = 2:(nrow(tst3)+1), 
                                rule = "==0", 
                                style = createStyle(bgFill = drac[11])
          )
          conditionalFormatting(wb, xl.name,
                                cols = 2:ncol(tst3),
                                rows = 2:(nrow(tst3)+1), 
                                rule = "==1", 
                                style = createStyle(bgFill = drac[2])
          )
          conditionalFormatting(wb, xl.name,
                                cols = 2:ncol(tst3),
                                rows = 2:(nrow(tst3)+1), 
                                rule = "==2", 
                                style = createStyle(bgFill = drac[1])
          )
          setColWidths(wb, xl.name, cols = 2:ncol(tst3), widths = 1)
          setColWidths(wb, xl.name, cols = 1, widths = 25)
          setRowHeights(wb, xl.name, rows = 1:(nrow(tst3)+1), heights = 15)
          freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
          
          saveWorkbook(wb, file)
      }
    )
    
  })
  
  #########  
  #########  
  
}




























