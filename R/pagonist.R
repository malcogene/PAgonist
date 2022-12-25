
options(stringsAsFactors = FALSE)

AppUI <- function() {
  currentfilePath <- dirname(rstudioapi::getSourceEditorContext()$path)
  pkgname <-(pkg <- strsplit(currentfilePath, '/'))[[1]][length(pkg[[1]])-1]; pkgname
  shiny::runApp(system.file("shiny", package=pkgname))} #
appU <- function() { loc <- gsub('.*:', '', getAnywhere("AppUI")$where[1])
shiny::runApp(system.file("shiny", package=loc))  }




#' @name cdcMap
#' @title ...
#' @description  ...
#' @param sfPower numerical value of ... for
#' @return A user interface will be shown in users' default web browser.
#' @import shiny ggplot2 ggrepel
#' @export
cdcMap <-
  function(rankedNamedVec = {
    rankedNamedVec = c(13:-12)
    names(rankedNamedVec) = letters
    rankedNamedVec
  },
  subsetOfVec = sample(letters, 5),
  concordantDirected = T,
  scaleFreeConstant = 1,
  dnPenalizedPowerConstant = 2,
  plot = T,
  theme_c = F,
  offTargetLabel = T,
  offTargetMaxN = 2,
  offTargetEnrichmentAnalysis = F,
  offTargetE.pw.in = NULL,
  legend.position = "top",
  verbose = T)
  {

    par(lend = 2, ljoin = 1)
    positiveToNegativeInx = which(sign(rankedNamedVec <-
                                         sort(rankedNamedVec, decreasing = T)) < 0)[1]
    subsetOfVec <- intersect(subsetOfVec, names(rankedNamedVec))
    PDF.Hit <- PDF.Miss <- rep(0, length(rankedNamedVec))
    Hits <- names(rankedNamedVec) %in% subsetOfVec
    # if( sum(rankedNamedVec[Hits]) == 0 ) return("sum(rankedNamedVec[Hits]) == 0")
    PDF.Hit[Hits] <- abs(rankedNamedVec[Hits]) ^ scaleFreeConstant
    PDF.Hit <-  PDF.Hit / sum(PDF.Hit)
    CDF.Hit <- cumsum(PDF.Hit)
    PDF.Miss[!Hits] <-
      1 / (length(rankedNamedVec) - length(subsetOfVec))
    CDF.Miss <- cumsum(PDF.Miss)
    deltaD.df = data.frame(
      Rank = 1:length(CDF.Miss),
      Name = names(rankedNamedVec),
      Value = rankedNamedVec,
      Hit = as.integer(Hits),
      CDF.miss = CDF.Miss,
      CDF.Hit = CDF.Hit,
      deltaD = (CDF.Hit - CDF.Miss)
    )
    deltaD.df$maxD = ifelse(abs(deltaD.df$deltaD) == max(abs(deltaD.df$deltaD)), 1, 0)
    Max <- deltaD.df$maxD == 1

    if (concordantDirected) {
      if (!is.na(positiveToNegativeInx)) {
        positiveHit.CDF = sum(PDF.Hit[1:(positiveToNegativeInx - 1)])
        negativeHit.CDF = sum(PDF.Hit[positiveToNegativeInx:length(PDF.Hit)])
        discordanceConstant <-
          (1 + (positiveHit.CDF - negativeHit.CDF) / (positiveHit.CDF + negativeHit.CDF)) /
          2
        if (verbose)
          cat(
            sprintf(
              "\n------------------------------\npreDiscordance Constant: %s",
              round(discordanceConstant, 3)
            )
          )

        # 2*atan(dnPenalizedPowerConstant)/pi +.5
        if (deltaD.df$deltaD[Max] > 0) {
          discordanceConstant = discordanceConstant ^ dnPenalizedPowerConstant
          if (verbose)
            cat(
              sprintf(
                "\n------------------------------\nDown-penalized Constant: %s",
                round(discordanceConstant, 3)
              )
            )
        } else {
          discordanceConstant = 1 - discordanceConstant
        }
        if (verbose)
          cat(
            sprintf(
              "\n------------------------------\npostDiscordance Constant: %s",
              round(discordanceConstant, 3)
            )
          )
      } else {
        discordanceConstant <- 1
      }
    } else {
      discordanceConstant <- 1
    }
    deltaD.df$adj.deltaD <- discordanceConstant * deltaD.df$deltaD
    if (verbose)
      cat(
        sprintf(
          "\n-------------------------------\nEnrichment Score: %s\nDirection penalized Enrichment Score: %s\n\n\n",
          round(deltaD.df$deltaD[Max], 3),
          round(deltaD.df$adj.deltaD[Max], 3)
        )
      )


    if (deltaD.df$adj.deltaD[Max] > 0) {
      dataS = subset(deltaD.df, (Rank > Rank[positiveToNegativeInx - 1]) &
                       (Hit == 1))
    } else {
      dataS = subset(deltaD.df, (Rank < Rank[positiveToNegativeInx + 1]) &
                       (Hit == 1))
    }

    if (verbose)
      print(dataS$Name)
    if (offTargetEnrichmentAnalysis) {
      if (!is.null(offTargetE.pw.in)) {
        offTargetERes <- ORA(dataS$Name, pw.in = offTargetE.pw.in)
      } else {
        offTargetERes <- ORA(dataS$Name)
      }
      if (verbose)
        print(offTargetERes)
    } else {
      offTargetERes <- NULL
    }

    if (plot) {
      suppressWarnings({
        p <-
          ggplot(deltaD.df, aes(x = Rank, y = CDF.miss))  # + geom_area(aes(x=Rank, y=deltaD, fill="Up"))
        if (!is.na(positiveToNegativeInx)) {
          p <-
            p + geom_ribbon(
              data = subset(deltaD.df, Rank < Rank[positiveToNegativeInx + 1]),
              aes(ymax = deltaD, fill = "Up"),
              ymin = 0,
              alpha = .3
            ) + geom_ribbon(
              data = subset(deltaD.df, Rank < Rank[positiveToNegativeInx + 1]),
              aes(ymax = adj.deltaD, fill = "Up"),
              ymin = 0,
              alpha = .7
            ) + geom_ribbon(
              data = subset(deltaD.df, Rank > Rank[positiveToNegativeInx - 1]),
              aes(ymax = deltaD, fill = "Down"),
              ymin = 0,
              alpha = .3
            ) + geom_ribbon(
              data = subset(deltaD.df, Rank > Rank[positiveToNegativeInx - 1]),
              aes(ymax = adj.deltaD, fill = "Down"),
              ymin = 0,
              alpha = .5
            )
        }
        if (all(names(rankedNamedVec) == letters)) {
          label = substitute(paste(italic(adj.D ~ statistic), " (", italic(KS), "): ", d), list(d =
                                                                                                  round(deltaD.df$adj.deltaD[Max], 3)))

          p <-
            p + geom_segment(
              aes(
                x = Rank,
                y = CDF.Hit,
                xend = Rank,
                yend = CDF.Miss
              ),
              color = ifelse(Hits, "#FFFFFF00", "#00000030"),
              size = .4,
              linetype = ifelse(Hits, 1, 2)
            ) + geom_point(color = "grey50", size = ifelse(Hits, 0, 1.5)) + geom_point(
              aes(y = CDF.Hit),
              color = ifelse(Hits, "#3F814A", "grey50"),
              size = ifelse(Hits, 1.5, 0)
            ) + geom_text(
              aes(
                x = Rank,
                y = CDF.Hit,
                label = ifelse(Hit, Name, "")
              ),
              size = 3.2,
              hjust = -.4,
              vjust = -.5,
              color = "#3F814A"
            ) + geom_text(
              aes(
                x = Rank,
                y = CDF.Miss,
                label = ifelse(Hit, "", Name)
              ),
              size = 3.2,
              hjust = -.4,
              vjust = -.5,
              color = "grey50"
            ) + annotate(
              'text',
              x = deltaD.df$Rank[Max],
              y = -(ymax / 20),
              label = label,
              size = 2.5,
              color = "#101221"
            )
        }
        ymax <-
          ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2]
        p <-
          p + geom_segment(
            aes(
              x = Rank,
              y = CDF.Hit,
              xend = Rank,
              yend = CDF.Miss
            ),
            color = ifelse(Hits, "#3F814A", "#FFFFFF00"),
            size = .4,
            linetype = ifelse(Hits, 1, 2)
          ) + geom_segment(
            aes(
              x = Rank[Max],
              y = CDF.Hit[Max],
              xend = Rank[Max],
              yend = CDF.Miss[Max]
            ),
            color = "#244337",
            linetype = 1,
            size = 1.2
          ) + geom_line(aes(y = deltaD, color = "ES"),
                        linetype = 1,
                        size = .3) + geom_line(aes(y = CDF.miss, color = "refCDF")) + geom_line(aes(y =
                                                                                                      CDF.Hit, color = "ECDF")) + geom_line(aes(y = adj.deltaD, color = "CES"), size =
                                                                                                                                              1.15) + labs(y = "Probability", color = "") + scale_discrete_manual(
                                                                                                                                                aesthetics = c("color"),
                                                                                                                                                labels = c("GS eCDF", "non-GS eCDF", "Running Score (ES)", "DES"),
                                                                                                                                                values = c(
                                                                                                                                                  ECDF = "#53996B",
                                                                                                                                                  refCDF = "#00000030",
                                                                                                                                                  ES = "#72F94C",
                                                                                                                                                  CES = "#72F94C"
                                                                                                                                                )
                                                                                                                                              ) + scale_discrete_manual(
                                                                                                                                                aesthetics = c("fill"),
                                                                                                                                                labels = c("Up", "Down"),
                                                                                                                                                values = c(Up = "#F29BA1", Down = "#3382BA")
                                                                                                                                              ) + guides(fill = guide_legend("Direction")) + theme_grey() + theme(legend.position =
                                                                                                                                                                                                                    legend.position,
                                                                                                                                                                                                                  legend.key.size = unit(.35, 'cm'))

        if (offTargetLabel) {
          if (dim(dataS)[1] > offTargetMaxN) {
            if (deltaD.df$adj.deltaD[Max] > 0) {
              dataS2 <-
                dataS[dim(dataS)[1]:(dim(dataS)[1] - (offTargetMaxN - 1)),]
            } else {
              dataS2 <- dataS[1:offTargetMaxN,]
            }
          } else {
            dataS2 <- dataS
          }
          if (all(names(rankedNamedVec) != letters))
            dataS2$Name <- mapg(dataS2$Name)   # <-- pkg에서 삭제요!
          p <-
            p + geom_label_repel(
              data = dataS2,
              aes(y = CDF.Hit, label = Name),
              nudge_y = 2,
              fontface = 4,
              color = "white",
              fill = "red",
              box.padding = unit(0.25, "lines"),
              point.padding = unit(0.3, "lines"),
              size = 2.5,
              segment.size  = 0.1,
              segment.color = "grey40",
              direction     = "x"
            )
        }
        if (theme_c)
          p <-
          p + theme_centro(grid = T) + theme(legend.position = legend.position,
                                             legend.key.size = unit(.35, 'cm'))
        print(p)
      })
    }
    list(
      deltaD.df = deltaD.df,
      maxD = deltaD.df[deltaD.df$maxD == 1, "adj.deltaD"],
      plot = p,
      offTargetERes = offTargetERes,
      offTarget = dataS
    )
  }






#' @name cdcScore
#' @title ...
#' @description  ...
#' @param sfPower numerical value of ... for
#' @return A user interface will be shown in users' default web browser.
#' @import shiny ggplot2 ggrepel
#' @export
cdcScore <-
  function(pert.matrix,
           upSet = NULL,
           dnSet = NULL,
           concordantDirected = F,
           dnPenalizedPowerConstant = 1,
           plot = F,
           theme_c = F,
           tar.OffTargetE = NULL,
           offTargetE.pw.in = NULL,
           offTargetMaxN = 2) {
    if (!is.null(upSet))
      tmpUp <-
        apply(pert.matrix, 2, function(x)
          cdcMap(
            sort(x, decreasing = T),
            upSet,
            plot = plot,
            concordantDirected = concordantDirected,
            dnPenalizedPowerConstant = dnPenalizedPowerConstant,
            theme_c = theme_c,
            offTargetMaxN = offTargetMaxN,
            offTargetEnrichmentAnalysis = F
          )$maxD)
    if (!is.null(dnSet))
      tmpDn <-
        apply(pert.matrix, 2, function(x)
          cdcMap(
            sort(x, decreasing = T),
            dnSet,
            plot = plot,
            concordantDirected = concordantDirected,
            dnPenalizedPowerConstant = dnPenalizedPowerConstant,
            theme_c = theme_c,
            offTargetMaxN = offTargetMaxN,
            offTargetEnrichmentAnalysis = F
          )$maxD)

    if (is.null(offTargetE.pw.in)) {
      offTargetERes.Up <- offTargetERes.Dn <- NULL
    } else {
      if (is.null(tar.OffTargetE))
        tar.OffTargetE = colnames(pert.matrix)[1]
      if (!is.null(upSet))
        offTargetERes.Up <-
          apply(pert.matrix[, tar.OffTargetE, drop = F], 2, function(x)
            cdcMap(
              sort(x, decreasing = T),
              upSet,
              plot = F,
              concordantDirected = concordantDirected,
              dnPenalizedPowerConstant = dnPenalizedPowerConstant,
              theme_c = theme_c,
              offTargetEnrichmentAnalysis = T,
              offTargetE.pw.in = offTargetE.pw.in
            ))
      if (!is.null(dnSet))
        offTargetERes.Dn <-
          apply(pert.matrix[, tar.OffTargetE, drop = F], 2, function(x)
            cdcMap(
              sort(x, decreasing = T),
              dnSet,
              plot = F,
              concordantDirected = concordantDirected,
              dnPenalizedPowerConstant = dnPenalizedPowerConstant,
              theme_c = theme_c,
              offTargetEnrichmentAnalysis = T,
              offTargetE.pw.in = offTargetE.pw.in
            ))
    }


    if (!is.null(upSet) & !is.null(upSet)) {
      cdcScore = ifelse(sign(tmpUp) != sign(tmpDn), (tmpUp - tmpDn) / 2, 0)
    } else if (!is.null(upSet) & is.null(upSet)) {
      cdcScore = tmpUp
    } else {
      cdcScore = tmpDn
    }
    cdcScore = sort(cdcScore)
    # cat(gsub("__.*", "", names(cdcScore[1])))
    list(
      cdcScore = cdcScore,
      offTargetEResUp = offTargetERes.Up,
      offTargetEResDn = offTargetERes.Dn
    )
  }






















#' @export
is.installed <- function(RequiredPackages) {
  RequiredPackages = c('ggplot2', 'ggrepel')

  pinx <- which(RequiredPackages %in% installed.packages()[,1])
  if(length(pinx) !=0) {installPackages<- RequiredPackages[-pinx] };
  if(length(installPackages) !=0) {
    Inx <- readline(prompt= sprintf("\nThis function needs %s package(s). Whould you like to install?\n\nEnter Y or an empty line to skip install and return:\n\n", installPackages) );
    if( Inx == 'y' || Inx == 'Y' ) {
      for (i in installPackages) { # Installs packages if not yet installed
        # if (!is.element(i, installed.packages()[,1]))
        install.packages(i)
        require(i, character.only = T)
        # }
      } } else { stop() }
  } else {
    for (i in RequiredPackages) {
      require(i, character.only = T)
    }
  }
}









#' @export
is.installed.bioconductor <- function(RequiredPackages) {
  pinx <- which(RequiredPackages %in% installed.packages()[,1])
  pinx <- which(RequiredPackages %in% installed.packages()[,1])
  if(length(pinx) !=0) {installPackages<- RequiredPackages[-pinx] };
  if(length(installPackages) !=0) {
    Inx <- readline(prompt= sprintf("\nThis function needs %s bioconductor package(s). Whould you like to install?\n\nEnter Y or an empty line to skip install and return", installPackages) );
    if( Inx == 'y' || Inx == 'Y' ) {
      for (i in installPackages) { # Installs packages if not yet installed
        # if (!is.element(i, installed.packages()[,1])) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install(i)
        require(i, character.only = T)
        # }
      } } else { stop() }
  } else {
    for (i in RequiredPackages) {
      require(i, character.only = T)
    }
  }
}



#' @export
loadUrl <- function(url, downloadPath = NA, sep=c("RData"," ", "," , "\t", ";", "xls", "gsheet"), ...) {
  cat('onedrive: copy link\n googlesheet: share-> Anyone with the link\n sep: "RData", ..."xls", "gsheet"\n')
  if(!is.na(downloadPath))  { tmpFile <- downloadPath

  } else { tmpFile <- file.path(getwd(), paste0(substr(Sys.time(), 1, 10), '.rda' ))  }
  url2 <- gsub("e=.*", "download=1", url)
  download.file(url2, tmpFile, mode="wb") # For windows user: mode = "wb" cf. binary
  sep <- match.arg(sep)
  if(sep == "RData") {
    print(tmpFile)
    tmpFile <-  gsub("\\\\", "/", tmpFile)
    justLoaded <- try(load(tmpFile), silent = T);
    try(assign(justLoaded, eval(as.symbol(justLoaded)),.GlobalEnv ), silent = T);
    if(class(justLoaded)=="try-error"){ justLoaded <- try(read.delim(tmpFile, ...), silent = T); message("Need 'sep' argument, is it txt file?")  }
  } else if(sep == "xls") {
    is.installed('readxl')
    justLoaded <- try(read_excel(tmpFile,...), silent = T)

  } else if(sep == "gsheet") {
    is.installed('gsheet')
    cat('gsheet should be public, click share-> Anyone with the link')
    justLoaded <- gsheet2tbl(url,...)
  } else {
    justLoaded <- try(read.delim(tmpFile, sep=sep, ...), silent = T)
  }
  justLoaded
}



#' @export
browseEntrez <- function(entrezIDs) {
  for(i in entrezIDs) {
    browseURL(paste0("https://www.ncbi.nlm.nih.gov/gene/", i))
  }
}


#' @export
peep <- function(x, boxplot = F ) {
  if(is.null(dim(x))) { if(length(x) > 10)  { print(x[1:10]) } else { print(x) }  } else if (dim(x)[1] >=10 && dim(x)[2]>=5 ){ print(x[1:5, 1:3]); boxplot(x[1:5, 1:3]) } else {print(head(x)); boxplot(x)} }


#' @export
normalize.q <- function(x= data.frame(matrix(sample(12, replace = T), 4)), filter.sd.quantile = 0.1, tied = c("average", "min", "max"), verbose = T ) {
  # compare to normalize.quantiles, 1. accept data.frame 2. tie control option:"average", "min", "max"  3. sd.filter 4. peep & plot & verbose...

  x <- x[rowSums(x)>0, ]
  x <- x[apply(x,1,sd) >= quantile(apply(x,1,sd), filter.sd.quantile), ]
  cat(sprintf("\nrowSums(x)==0, =<quantile(sd(row), %s) were filtered\n\n", filter.sd.quantile))

  tied <- match.arg(tied)
  rank <- apply(x, 2, rank,ties.method="min");
  if(any(tied %in% c("average", "max"))) rank.max <- apply(x, 2, rank,ties.method="max");
  sorted <- apply(x, 2, sort)
  sorted.row.mean <- apply(sorted, 1, mean);
  x2 <- apply(rank, 2, function(x) sorted.row.mean[x])
  if(any(tied %in% c("average", "max"))) x2.max <- apply(rank.max, 2, function(x) sorted.row.mean[x])
  if(tied=="average") { x2 <- (x2+x2.max)/2 } else if (tied=="max"){x2 <- x2.max } else { }

  if( class(x) == "data.frame") { x2 <- as.data.frame(x2); rownames(x2) <- rownames(x) }
  if(verbose) {
    op <- par(no.readonly = T); par(mfrow=c(1,2), mar=c(3,3,1,1))
    cat('Original matrix or data.frame\n'); peep(x, T)
    cat('Sort the original matrix from lowest to highest\n'); peep(rank)
    cat('Determine the ranks of original matix\n');peep(sorted)
    cat('\nCalculate the means\n\n'); peep(sorted.row.mean)
    cat('\n\nFinally substitute the means into our ranked matrix\n'); peep(x2, T)
    cat(sprintf('If the values were tied, %s is used\n\n', tied))
    par(op)
    'In the example on Wikipedia, if the values were tied, the min value is used but in the normalize.quantiles() function, the average is used'
  }
  x2
}




#' @export
DEGs <- function(Exp, cl, adj.pval = 0.1,  logFC = 2, geomTextN=5, heatmapUpN = 25, plotDEG =T, multipleRegression=F, rowCounts=F, meanFilter=10, PDF=T, cname='temp' ) {
  try(dev.off(), silent = T)

  is.installed(c('ggplot2', 'ggrepel'))
  is.installed.bioconductor(c('limma', 'ComplexHeatmap'))


  if(rowCounts) { Exp <- Exp[apply(Exp, 1, mean) > meanFilter, ]; Exp <- voom(Exp, plot = T) }

  if(multipleRegression) { fit <-eBayes(lmFit(Exp, model.matrix(~ .,cl))); print(topTable(fit, 2))                                     } else {
    fit <-eBayes(lmFit(Exp, model.matrix(~ .,cl[, 1, drop=F]))); print(topTable(fit, 2)) }

  tT <- topTable(fit, number = dim(Exp)[1])
  tT$Gene <- rownames(tT)
  tT.up <- tT[order(tT$logFC, decreasing = T ),]; tT.down<- tT[order(tT$logFC),]
  tT.filter <- data.frame(tT[inx<-(tT$adj.P.Val<adj.pval) & (abs(tT$logFC) > logFC ), ]); print(tT.filter)

  if(PDF) {
    pdf(file = file.path(getwd(),sprintf("%s.pdf", cname )), width=5, height=5)  }
  if(plotDEG) {
    if(any(colnames(tT) == "logFC" && dim(tT.filter)[1] != 0) ) {
      require(ggplot2); require(ggrepel)
      tT$Cutoff_value <- c("Not Sig", sprintf("FDR < %s & logFC > %s", adj.pval, logFC))[as.numeric(inx)+1]
      gplot <- ggplot(tT, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point(aes(color = Cutoff_value)) + labs(title ="c") + scale_color_manual(values = c("red", "grey")) +  theme_bw(base_size = 12) + theme(legend.position = "bottom") + geom_hline(yintercept= -log10(adj.pval), linetype="dashed", color = "#FF000050") + geom_vline(xintercept= c(logFC, -logFC), linetype="dashed", color = "#FF000050")

      g <- gplot + geom_text_repel(
        data = dd <- rbind(tT.up[1:geomTextN, ], tT.down[1:geomTextN, ]),
        aes(label = Gene),
        size = 3,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      ) }
    print(g)

    if(dim(Exp)[1] >= heatmapUpN*2 ) {
      bluered <- colorRampPalette(c("blue", "white", "red"))(256)
      # stats::heatmap(Exp[rbind(tT.down[1:heatmapUpN, ],tT.up[heatmapUpN:1, ])$Gene,], col = bluered, scale = "row", main = sprintf("top%s", heatmapUpN*2), Colv = NA, Rowv=NA    )

      # if(HUGO) { colnames(d)[1:dim(matGS)[1]] <- .mapid(colnames(d)[1:dim(matGS)[1]]) }

      colse=c("#00000020", "#000000", "#0000BF10", "#0000BF30", "#0000BF50", "#0000BF70","#0000BF90","#0000BF")
      colTemp <- colse[as.numeric(as.factor(cl[,1]))]
      names(colTemp ) <- cl[,1]
      colTemp<-list(colTemp); names(colTemp) <- colnames(cl)[1]


      h <- Heatmap( t(scale(t(d<-Exp[rbind(tT.up[1:heatmapUpN, ], tT.down[heatmapUpN:1, ])$Gene,]))),  col = bluered, name="Exprs", rect_gp = gpar(col = "gray12", lty = 1, lwd = 0.2), cluster_rows = T, cluster_columns = T, show_row_names = T, row_names_gp =gpar(fontsize = 5),   split = data.frame(cyl = factor(c(rep('UP', heatmapUpN), rep('DOWN', heatmapUpN)), levels=c('UP','DOWN' ))),gap = unit(1.5, "mm"), top_annotation = HeatmapAnnotation(df=cl, col=colTemp, show_annotation_name = T, annotation_name_gp = gpar(fontsize = 7), annotation_name_side ='left', annotation_height=c(1.5) ) )

      pre.h<<-h
      draw(h)
    }
  }
  if(PDF) { dev.off()
  }
  return(list(fit = fit, tT.filter= tT.filter, tT.up=tT.up, tT.down=tT.down))
}


