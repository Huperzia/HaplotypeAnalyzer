############
# getUtils #
############
drac <- c("#50fa7b", "#ffb86c", "#bd93f9", "#ff79c6", 
          "#ff5555", "#f1fa8c", "#6272a4", "#8be9fd", 
          "#282a36", "#44475a", "#44475a", "#f8f8f2")
gradient <- colorRampPalette(drac[c(1,6,2,3)])
#loadfonts(device = "pdf", quiet = TRUE)

font <- list(
  family = "Harding Text Web Regular",
  size = 14,
  color = drac[12])

V <- theme_bw() +
  theme(legend.position = c(0.75, 0.75),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_rect(fill = drac[9], color = drac[9]),
        panel.background = element_rect(fill = drac[9]),
        legend.background = element_rect(fill = drac[9]),
        axis.text = element_text(color = drac[12], size = 18),
        text = element_text(color = drac[12], size = 24,
                            family="Harding Text Web Regular"
        ),
        panel.border = element_rect(color = drac[12]),
        axis.ticks = element_line(color = drac[12]),
        legend.key = element_rect(fill = drac[9], color = drac[9]))

sc <- scale_x_continuous(labels = label_comma(accuracy = 0.1, scale = 1/1000, suffix = "  kb"))


# BCFS
getbcfs <- function(chr = chr, RANGE = RANGE, VCFID = VCFID, VCFS = VCFS){
  
  bcfs <- foreach(i = c(1:1),
                  .combine = "c") %dopar% {
                    FILE <- paste(chr, ".", RANGE[1], "..", RANGE[2], ".", VCFID[i], ".vcf.gz", sep = "")
                    
                    FILES <- list.files("tmp/", pattern = "012")
                    if(!(paste(FILE, ".012", sep ="") %in% FILES)){
                      system.time(system(paste("bcftools view ", VCFS[i],
                                               " -r ", chr, ":", RANGE[1], "-", RANGE[2],
                                               " --threads 8",
                                               " -Oz | ",
                                               "vcftools --gzvcf -",
                                               " --012",
                                               " --out tmp/", FILE,
                                               sep = "")))
                    }
                    
                    tmp <- t(read.table(paste("tmp/", FILE, ".012", sep =""), header=F)[,-1])
                    colnames(tmp) <- read.table(paste("tmp/", FILE, ".012.indv", sep =""))[,1]
                    rownames(tmp) <- read.table(paste("tmp/", FILE, ".012.pos", sep =""))[,2]
                    
                    list(tmp)
                  }
  
  return(bcfs)
}

# CNTS
getcnts <- function(bcfs = bcfs, key = key){
  
#  cnts <- #foreach(i = c(1:1),
           #       .combine = "c") %dopar% {
                    df <- bcfs[[i]]
                    tmp <- data.frame(matrix(NA, nrow = dim(df)[2], ncol = 4))
                    
                    for(ii in 1:dim(df)[2]){
                      tmp[ii,] <- c(sum(df[,ii] == -1), sum(df[,ii] == 0), sum(df[,ii] == 1), sum(df[,ii] == 2))
                    }
                    
                    colnames(tmp) <- c("missing", "ref", "het", "alt")
                    tmp <- cbind(id = colnames(df), tmp, stringsAsFactors = F)
                    tmp <- merge(tmp, key, by = "id")
                    
                    list(tmp)
               #   }
  
  return(cnts)
}

# PCA
getpcas <- function(bcfs = bcfs, cnts = cnts){
  
#  pcas <- foreach(i = c(1:1),
 #                 .combine = "c") %dopar% {
                    df <- bcfs[[i]]
                    X = as.matrix(dist(t(as.matrix(df))))
                    tmp <- prcomp(X)
                    
                    #imp <- percent((summary(pca)$importance/sum(summary(pca)$importance))[1:2], accuracy = 1e-2)
                    tmp <- cbind(id = rownames(X), data.frame(tmp$x), stringsAsFactors= F)
                    tmp <- tmp[, 1:4]
                    tmp.cnts <- merge(cnts[[i]], tmp, by = "id")
                    
                    list(list(tmp, tmp.cnts))
   #               }
  
  #return(pcas)
}

# DAPC
getdapcs <- function(bcfs = bcfs, cnts = cnts, ldas = ldas, nclust = nclust){
  
  dapcs <- foreach(i = c(1:1),
                   .combine = "c",
                   .packages = c("adegenet")) %dopar% {
                     df <- bcfs[[i]]
                     set.seed(1000)
                     x <- new("genind", as.matrix(t(df)), ind.names=colnames(df))
                     
                     grp <- find.clusters(x,
                                          #max.n.clust = 9,
                                          stat = "BIC",
                                          n.clust = nclust,
                                          n.pca = 300, scale = F)
                     
                     tmp <- dapc(x, grp$grp, n.pca = 300, n.da = ldas)
                     tmp <- data.frame(dapc.grp = tmp$grp, dapc.1 = tmp$ind.coord[,1], dapc.2 = tmp$ind.coord[,2], dapc.3 = tmp$ind.coord[,3])
                     tmp$id <- colnames((df))
                     
                     tmp.cnts <- merge(cnts[[i]], tmp, by = "id")
                     list(list(tmp, tmp.cnts))
                   }
  
  return(dapcs)
}

# GENES
getGenes <- function(range){
  
}


#############
# plotUtils #
#############

plotFreqs <- function(cnts = cnts, alleles = alleles, leg = leg, bins = bins, alphas = alphas, cols = cols, contin = contin){
  p <- list()
  
  if(contin == T){
    scm <- scale_color_gradientn(colors = drac[c(1,6,2,3)], guide = F)
    filled <- cut(alleles, breaks = 4)
    alphas <- rep(0.6, length(unique(filled)))
  } else {
    scm <-  scale_color_manual(values = cols, guide = F)
    filled <- alleles
  }

#  alleles <- cnts$PC1    
#  alleles <- cnts$shrunken
  p[[1]] <- ggplot(cnts, aes(x = missing,
                             alpha = filled,
                             fill = filled)) +
    geom_histogram(position="identity", bins = bins) +
    scale_fill_manual(values = cols, guide = F) +
    scale_alpha_manual(values = alphas, guide = F) +
    labs(fill = "Genotype", x = "Missing", y = "Frequency") + V
  
  p[[3]] <- ggplot(cnts, aes(x = missing,
                             y = het + alt,
                             color = alleles)) +
    geom_point(size = 0.6, alpha = 0.8) +
    scm +
    theme(legend.position = "top") + V +
    labs(fill = "", x = "Missing", y = "Het + Alt")
  
  p[[4]] <- ggplot(cnts, aes(x = het + alt,
                             alpha = filled,
                             fill = filled)) +
    geom_histogram(position="identity", bins = bins) +
    scale_fill_manual(values = cols, guide = F) +
    scale_alpha_manual(values = alphas, guide = F) +
    labs(fill = "Genotype", x = "Het + Alt", y = "Frequency") + V + coord_flip()
  
  p[[5]] <- ggplot(cnts, aes(x = missing,
                             y = depth,
                             color = alleles)) +
    geom_point(size = 0.6) +
    scale_alpha_manual(values = 0.7, guide = F) +
#    scale_color_manual(values = cols, guide = F) + 
    scm +
    V +
    theme(axis.ticks.length.y = unit(2, "mm")) +
    labs(fill = "", x = "Missing", y = "Depth")
  
  p[[6]] <- ggplot(cnts, aes(x = het + alt,
                             y = depth,
                             color = alleles)) +
    geom_point(size = 0.6) +
    scale_alpha_manual(values = 0.7, guide = F) +
    #scale_color_manual(values = cols, guide = F) + 
    scm +
    V +
    theme(axis.ticks.length.y = unit(2, "mm")) +
    labs(fill = "", x = "Het + Alt", y = "Depth")
  
  p[[7]] <- cowplot::ggdraw(plot_grid(p[[1]], leg, p[[3]], p[[4]], p[[5]], p[[6]], ncol = 2, align='v', axis="tblr")) +
    theme(plot.background = element_rect(fill=drac[9], color = NA))
  
  
  return(p[[7]])
#  subplot(ggplotly(p[[1]]), ggplotly(leg), 
#          ggplotly(p[[3]]), ggplotly(p[[4]]), 
#          ggplotly(p[[5]]), plot_ly(), nrows = 3, shareX = T)
  
#  return(ggplotly(p[[2]]))
}

plotlyFreqs <- function(cnts = cnts, alleles = alleles, font = font, cols = cols){
  trait <- alleles
  plot_ly(data = cnts, x = ~missing, y = ~het, z = ~alt, color = ~trait, colors = cols,
          marker = list(size = 1),
          text = ~paste("ID: ", id,
                        "<br>Genotype:", trait,
                        sep = " "),
          hovertemplate = paste('%{text}<extra></extra>')) %>%
    layout(title = "",
           font = font,
           titlefont = list(size = 15, face = "bold", color = drac[12]),
           font = list(size = 14, margin = 10, face = "bold", color = drac[12]),
           plot_bgcolor = drac[9],
           paper_bgcolor = drac[9],
           legend = list(orientation = 'h'),
           scene = list(xaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12]),
                        yaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12]),
                        zaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12])))
}

plotPCAs <- function(cnts = cnts, alleles = alleles, leg = leg, bins = bins, alphas = alphas, cols = cols, contin = contin){
  p <- list()
  
  if(contin == T){
    scm <- scale_color_gradientn(colors = drac[c(1,6,2,3)], guide = F)
    filled <- cut(alleles, breaks = 4)
    alphas <- rep(0.6, length(unique(filled)))
  } else {
    scm <-  scale_color_manual(values = cols, guide = F)
    filled <- alleles
  }
  
  p[[7]] <- ggplot(cnts, aes(x = PC1,
                             alpha = filled,
                             fill = filled)) +
    geom_histogram(position="identity", bins = bins) +
    scale_fill_manual(values = cols, guide = F) +
    scale_alpha_manual(values = alphas, guide = F) +
    labs(fill = "Genotype", x = "PC1", y = "Frequency") + V
  
  p[[8]] <- ggplot(cnts, aes(x = PC1,
                             y = PC2,
                             color = alleles)) +
    geom_point(size = 0.6) +
    scale_alpha_manual(values = 0.7, guide = F) +
    scm +
    V +
    labs(fill = "", x = "PC1", y = "PC2")
  
  p[[9]] <- ggplot(cnts, aes(x = PC2,
                             alpha = filled,
                             fill = filled)) +
    geom_histogram(position="identity", bins = bins) +
    scale_fill_manual(values = cols, guide = F) +
    scale_alpha_manual(values = alphas, guide = F) +
    labs(fill = "Genotype", x = "PC2", y = "Frequency") + V + coord_flip()
  
  p[[10]] <- ggplot(cnts, aes(x = PC1,
                              y = depth,
                              color = alleles)) +
    geom_point(size = 0.6) +
    scale_alpha_manual(values = 0.7, guide = F) +
    scm + V +
    labs(fill = "", x = "PC1", y = "Depth")
  
  p[[11]] <- ggplot(cnts, aes(x = PC2,
                              y = depth,
                              color = alleles)) +
    geom_point(size = 0.6) +
    scale_alpha_manual(values = 0.7, guide = F) +
    scm + V +
    labs(fill = "", x = "PC2", y = "Depth")
  
  p[[12]] <- cowplot::ggdraw(plot_grid(p[[7]], leg, p[[8]], p[[9]], p[[10]], p[[11]], ncol = 2, align='v', axis="tblr")) +
    theme(plot.background = element_rect(fill=drac[9], color = NA))
  
  return(p[[12]])
}

plotlyPCAs <- function(cnts = cnts, alleles = alleles, font = font, cols = cols){
  trait <- alleles
  plot_ly(data = cnts, x = ~PC1, y = ~PC2, z = ~PC3, color = ~trait, colors = cols,
          marker = list(size = 1),
          text = ~paste("ID: ", id,
                        "<br>Genotype:", trait,
                        sep = " "),
          hovertemplate = paste('%{text}<extra></extra>')) %>%
    layout(title = "",
           legend = list(orientation = 'h'),
           font = font,
           titlefont = list(size = 24, face = "bold", color = drac[12]),
           #                       font = list(size = 14, margin = 10, face = "bold", color = drac[12]),
           plot_bgcolor = drac[9],
           paper_bgcolor = drac[9],
           scene = list(xaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12]),
                        yaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12]),
                        zaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12])))
  
}


plotDAPCs <- function(cnts = cnts, alleles = alleles, leg = leg, bins = bins, alphas = alphas, cols = cols, contin = contin){
  p <- list()
  cols <- drac[c(1,4,3,2,6)]
  
  if(contin == T){
    scm <- scale_color_gradientn(colors = drac[c(1,6,2,3)], guide = F)
    filled <- cut(alleles, breaks = 4)
    alphas <- rep(0.6, length(unique(filled)))
  } else {
    scm <-  scale_color_manual(values = cols, guide = F)
    filled <- alleles
  }
  
  p[[7]] <- ggplot(cnts, aes(x = dapc.1,
                             alpha = filled,
                             fill = filled)) +
    geom_histogram(position="identity", bins = bins) +
    scale_fill_manual(values = cols, guide = F) +
    scale_alpha_manual(values = alphas, guide = F) +
    labs(fill = "Genotype", x = "DAPC1", y = "Frequency") + V
  
  p[[8]] <- ggplot(cnts, aes(x = dapc.1,
                             y = dapc.2,
                             color = alleles)) +
    geom_point(size = 0.6, alpha = alphas[1]) +
    scale_alpha_manual(values = 0.7, guide = F) +
    scm +  V +
    labs(fill = "", x = "DAPC1", y = "DAPC2")
  
  p[[9]] <- ggplot(cnts, aes(x = dapc.2,
                             alpha = filled,
                             fill = filled)) +
    geom_histogram(position="identity", bins = bins) +
    scale_fill_manual(values = cols, guide = F) +
    scale_alpha_manual(values = alphas, guide = F) +
    labs(fill = "Genotype", x = "DAPC2", y = "Frequency") + V + coord_flip()
  
  p[[10]] <- ggplot(cnts, aes(x = dapc.1,
                              y = depth,
                              color = alleles)) +
    geom_point(size = 0.6, alpha = alphas[1]) +
    scale_alpha_manual(values = 0.7, guide = F) +
    scm +  V +
    labs(fill = "", x = "DAPC1", y = "Depth")
  
  p[[11]] <- ggplot(cnts, aes(x = dapc.2,
                              y = depth,
                              color = alleles)) +
    geom_point(size = 0.6, alpha = alphas[1]) +
    scale_alpha_manual(values = 0.7, guide = F) +
    scm + V +
    labs(fill = "", x = "DAPC2", y = "Depth")
  
  p[[12]] <- cowplot::ggdraw(plot_grid(p[[7]], leg, p[[8]], p[[9]], p[[10]], p[[11]], ncol = 2, align='v', axis="tblr")) +
    theme(plot.background = element_rect(fill=drac[9], color = NA))
  
  return(p[[12]])
}

plotlyDAPCs <- function(cnts = cnts, alleles = alleles, font = font, cols = cols){
  trait <- alleles
  ply <- plot_ly(data = cnts, x = ~dapc.1, y = ~dapc.2, z = ~dapc.3, color = ~trait, colors = cols,
                 marker = list(size = 1),
                 text = ~paste("ID: ", id,
                               "<br>Genotype:", trait,
                               "<br>DAPC Group:", dapc.grp,
                               sep = " "),
                 hovertemplate = paste('%{text}<extra></extra>')) %>%
    layout(title = "",
           legend = list(orientation = 'h'),
           font = font,
           titlefont = list(size = 15, face = "bold", color = drac[12]),
           #                       font = list(size = 14, margin = 10, face = "bold", color = drac[12]),
           plot_bgcolor = drac[9],
           paper_bgcolor = drac[9],
           scene = list(xaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12]),
                        yaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12]),
                        zaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12])))
  return(ply)
}


plotlyTraits <- function(cnts = cnts, alleles = alleles, font = font, cols = cols){
#  alleles2 <- (40 - alleles)/4
  trait <- alleles
  plot_ly(data = cnts, x = ~PC1, y = ~PC2, z = ~PC3, color = ~trait, colors = gradient(30),
          marker = list(size = ~alleles),
          text = ~paste("ID: ", id,
                        "<br>Trait:", trait,
                        sep = " "),
          hovertemplate = paste('%{text}<extra></extra>')) %>%
    layout(title = "",
           legend = list(title=list(text='<b> Trait </b>'), orientation = 'h'),
           font = font,
           titlefont = list(size = 24, face = "bold", color = drac[12]),
           #                       font = list(size = 14, margin = 10, face = "bold", color = drac[12]),
           plot_bgcolor = drac[9],
           paper_bgcolor = drac[9],
           scene = list(xaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12]),
                        yaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12]),
                        zaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12])))
  
}

plotlyTraitsDAPC  <- function(cnts = cnts, alleles = alleles, font = font, cols = cols, nclust = nclust){
  trait <- alleles
  ply <- plot_ly(data = cnts, x = ~dapc.1, y = ~dapc.2, z = ~dapc.3, color = ~dapc.grp, colors = drac[1:nclust],
                 marker = list(size = ~trait),

                 text = ~paste("ID: ", id,
                               "<br>Trait:", trait,
                               "<br>DAPC Group:", dapc.grp,
                               sep = " "),
                 hovertemplate = paste('%{text}<extra></extra>')) %>%
    layout(title = "",
           legend = list(orientation = 'h'),
           font = font,
           titlefont = list(size = 15, face = "bold", color = drac[12]),
           #                       font = list(size = 14, margin = 10, face = "bold", color = drac[12]),
           plot_bgcolor = drac[9],
           paper_bgcolor = drac[9],
           scene = list(xaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12]),
                        yaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12]),
                        zaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12])))
  return(ply)
}


plotlyTraitsDAPCdist  <- function(cnts = cnts, alleles = alleles, font = font, cols = cols, nclust = nclust){
  trait <- alleles
  dapc.group <- jitter(as.numeric(cnts$dapc.grp))
  ply <- plot_ly(data = cnts, 
                 x = dapc.group, 
                 y = ~trait, 
                 color = ~dapc.grp, 
                 colors = drac[1:nclust],
                 type = "scatter",
                size = ~trait,
                 text = ~paste("ID: ", id,
                               "<br>Trait:", trait,
                               sep = " "),
                 hoverinfo = "text",
                 hovertemplate = paste('%{text}<extra></extra>')) %>%
    layout(title = "",
           legend = list(orientation = 'h'),
           font = font,
           titlefont = list(size = 15, face = "bold", color = drac[12]),
           plot_bgcolor = drac[9],
           paper_bgcolor = drac[9],
           scene = list(xaxis = list(title = "DAPC Cluster", showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', dtick = 1, showticklabels = TRUE, color = drac[12]),
                        yaxis = list(title = "Trait", showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12]),
                        zaxis = list(showgrid = TRUE, showline = TRUE, autotick = TRUE, ticks = '', showticklabels = TRUE, color = drac[12])))
  return(ply)

}
#~~~~~~~~~~~#
