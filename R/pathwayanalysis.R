
### to do items

# fix dependcy calls...

# add dependencies to description:
# tm, wordcloud, clusterprofiler, enrichplot, ggnewscale(?); others?





# survival analysis: forest plots

# also, in deplots.r, make sure to check heatmapplot:
# added suport for sig vs nonsig genes --> test it
# need to document some of the new params in that function


#### pathway analysis, clusterprofiler wrapper ####

#' Run Overrepresentation Analysis via the ClusterProfiler workflow
#'
#' @param sigres data.frame of significant genes. first column is gene names; second column is logfoldchange or other effect size.
#' @param term2gene data.frame of pathways to check. first column is pathways, second column is genes.
#'
#' @return a list object: first element is a ClusterProfiler object; second element is an "emap plot", a graph of pathways connected by jaccard similarity; third element is a data.frame with the connected nodes from the emap plot; 4th element is a list of pathways in each connected cluster; 5th element is genes from the pathways in each cluster
#' @export
#'
#' @examples
pathwayanalysis <- function(sigres,
                            term2gene){


  message('\nInitiating analysis')


  message('\nRunning clusterProfiler::enricher...')

  msigdb_ora <- clusterProfiler::enricher(sigres[,1],
                                          TERM2GENE = term2gene,
                                          # TERM2NAME = term2name,
                                          # OrgDb = org.Mm.eg.db,
                                          # keyType = 'ENSEMBL',
                                          pvalueCutoff = 0.05
  )

  message('\nPrepping emapplots...')

  #run termsimilarity
  # use all significant?
  # sigtab <- table(msigdb_ora@result$p.adjust < 0.05)
  # msigdb_ora <- enrichplot::pairwise_termsim(msigdb_ora, showCategory = sigtab['TRUE'] )

  #go with default, use onyl top 200
  msigdb_ora <- enrichplot::pairwise_termsim(msigdb_ora)


  #try to plot color = num genes in overexp gene list?
  msigdb_ora@result$Percent_of_DEGs <- (msigdb_ora@result$Count / nrow(sigres))*100


  set.seed(2021)
  options(ggrepel.max.overlaps = 2000)

  suppressMessages(
    emap_total <- enrichplot::emapplot(msigdb_ora, color='Percent_of_DEGs', repel=T, showCategory = 200,
                           layout='graphopt', cex_label_category=0.4, cex_category = 0.3,cex_line = 0.3)+
      scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent of DEGs')
  )



  #get termsim and min_dist
  mat <- msigdb_ora@termsim
  min_dist = 0.2 #default from package


  #keep only termsim matrix elements that are in final plot
  # the top 200 by pvalue
  dat <- emap_total$data


  #try to recapture clusters...
  message('\nRecovering connected nodes...')

  conectlist <- list()
  for(i in c(rownames(mat)) ){

    #i don't know why the cols and rows are not the same
    # just check both
    x = na.omit( mat[i,] )
    xnames <- names(x[x>=min_dist])

    y = na.omit( mat[,i] )
    ynames <- names(y[y>=min_dist])

    conectlist[[i]] <- unique(c(i,xnames, ynames))

  }



  #select single=node clusters, remove from list
  # singlenodes <- names( conectlist[which(sapply(conectlist, length)==0)] )
  #
  # conectlist <- conectlist[!(names(conectlist) %in% singlenodes)]


  #biggest first
  # sort by number of pathways in each cluster
  conectlist <- conectlist[order(sapply(conectlist, length), decreasing = T)]


  #pairwise overlap?
  # start with first cluster
  reslist <- conectlist[1]

  for(pway in names(conectlist)[-1] ){


    # message(pway)
    cpws <- conectlist[[pway]]

    #thru cpws and reslist
    # if it's in, add it
    # if not, add new to reslist
    for(residx in c(1:length(reslist)) ){


      respways <- reslist[[residx]]

      # check if any of the cpws is in res for this residx and break from this subloop
      if( any(cpws %in% respways) ){

        # message(' - matches with res ', residx)
        reslist[[residx]] <- unique(c(respways, cpws))

        break()

      }


      #if not, check if there is another reslist index...
      # if not, add
      if(residx == length(reslist)){

        # message(' - - - no match with res ', residx, ', adding another res index')
        reslist[[pway]] <- cpws

      }


    }


  }




  #it seems this can sometimes fail... not sure why...


  #try to collapse lists wiht overlap
  # maybe we could have done it this way above too...?

  backup <- reslist
  reslist2 <- list()
  skiplist <- c()

  for(residx in c(1:length(reslist)) ){

    # first, we check if we have previously already added the list...
    if( names(reslist)[residx] %in% skiplist ){next()}

    cpws <- reslist[[residx]]

    #check the lenght of intersects between this list nd all others
    # if there are, check the intersects of that list and all others and collapse recursively...
    # if not, just add to list

    intersectlens <- sapply(reslist[-residx], function(x){length( intersect(x, cpws) )})

    if( any(intersectlens > 0) ){


      # message('Intersect found for reslist index : ', residx, '; intersecting indices:')

      which_have_intersects <- names(intersectlens[intersectlens>0])

      for(j in which_have_intersects){

        intersectingidx <- which(names(reslist)==j)
        # message(' - reslist index: ',  intersectingidx)

        #add to skiplist
        skiplist <- c(skiplist, j)

        #update the reslist to include the intersecting pathway...
        cpws <- unique( c(cpws, reslist[[j]]) )
      }


    }

    # add to reslist2
    # only the intersects which were added get ignored...
    reslist2[[residx]] <- cpws


  }


  reslist <- reslist2
  rm(reslist2)

  #finally, we have the result!
  # totally agnostic to cluster number


  #sort by number of pathways in each cluster
  reslist <- reslist[order(sapply(reslist, length), decreasing = T)]

  #rename them
  names(reslist) <- paste0('Cluster_', seq(1:length(reslist)))





  #finally, do some actual analysis...
  # get total genes in cluster
  # get percent of DEGs from cluster

  clustpercs <- list()
  cum_perc <- 0
  clustgeneslist <- list()

  for(clustdex in 1:length(reslist) ){
    clust <- reslist[[clustdex]]

    #get the genes for each of these...
    clustergenes <- c()
    for(pw in clust){
      clustergenes <- c(clustergenes,
                        intersect(pull(term2gene[term2gene$gs_name ==pw,2] ) ,
                                  sigres[,1] )
      )

    }

    # get unique genes for this cluster
    # add them to list...
    clustergenes <- unique(clustergenes)
    clustgeneslist[[names(reslist)[clustdex]]] <- clustergenes

    #calculaate proportion:
    # length of unique genes in cluster / total significant DEGs
    perc <- (length(clustergenes) / nrow(sigres) * 100)

    #add to running perc
    cum_perc <- cum_perc + perc

    clustpercs[[clustdex]] <- data.frame(clust = names(reslist)[clustdex],
                                         perc = perc,
                                         num_genes = length(clustergenes),
                                         num_pways = length(reslist[[clustdex]]),
                                         cum_perc = cum_perc)



  }

  clustpercs <- bind_rows(clustpercs)

  #remove empty rows... not sure why these are a thing
  clustpercs <- clustpercs[clustpercs$num_genes>0,]


  #add cluster info to obj metadata
  #names(reslist) <- clustpercs$relabel
  pwayclust <- unlist2(reslist)
  pwayclust <- data.frame(Cluster=names(pwayclust),
                          ID = pwayclust)

  #add in names and %
  pwayclust$PercDEGs <- clustpercs[match(pwayclust$Cluster, clustpercs$clust ),'perc']

  #add clustering info to object res
  pwayres <- msigdb_ora@result
  pwayres <- pwayres[pwayres$ID %in% pwayclust$ID,]
  pwayres$Cluster <- pwayclust[match(pwayres$ID, pwayclust$ID), 'Cluster']
  pwayres$Percent_of_DEGs_CLUSTER <- pwayclust[match(pwayclust$ID, pwayres$ID), "PercDEGs"]



  #match to add into the msigdb_ora res
  msigdb_ora@result$Cluster <- NA ; msigdb_ora@result$Percent_of_DEGs_CLUSTER <- NA
  msigdb_ora@result[match(pwayres$ID, msigdb_ora@result$ID),] <- pwayres


  #add some extra plots of overall res
  #dot plot
  dotplot_all <- enrichplot::dotplot(msigdb_ora, showCategory=30,font.size=6)+
    ggtitle('Overall Pathways Dotplot') +
    theme(plot.title = element_text(hjust = 0.5))


  #wordcloud plot
  wc <- pathwayanalysis_wordcloud(unlist(reslist),
                                  scale=c(3,0.3), random.order=F, random.color=F,
                                  colors= RColorBrewer::brewer.pal(8, "Dark2"))
  title('Overall Pathways Wordcloud')
  wc <- recordPlot()


  plots <- list(emap_noclusters = emap_total,
                dotplot_all = dotplot_all,
                wordcloud_allpways = wc)


  output <- list(msigdb_ora_object = msigdb_ora,
                 plots = plots,
                 cluster_percentages = clustpercs,
                 pathways_in_each_cluster = reslist,
                 genes_in_each_cluster = clustgeneslist)

  output

}



#' A fast way to make a wordcloud from a character vector of pathway names
#'
#' This will make a wordcloud based on a character vector of pathway names.
#' If using with the rest of the pathway analysis pipeline, see the example below.
#'
#'
#' @param pways character vector. pathways to make a wordcloud with. or, more generally, a vector of "documents" with words separated by some delimiter.
#' @param delimiter a single character string. delimiter separating words in each of the elements of the vector `pways`
#' @param excludeWords words to exlude; for example, "GOBP"
#' @param ... other arguments passed on to `wordcloud::wordcloud()`
#'
#' @return a wordcloud  as a `recordedplot` plot object.
#' @export
#'
#' @examples
#'
#'
pathwayanalysis_wordcloud <- function(pways, delimiter, excludeWords, ...){

  #credit to http://www.sthda.com/english/wiki/word-cloud-generator-in-r-one-killer-function-to-do-everything-you-need
  # for basic source code and to help understand corpus object class

  if(missing(excludeWords)){excludeWords <- NULL}
  if(missing(delimiter)){delimiter <- '_'}

  #replace underscores with spaces
  pways <- gsub(delimiter,' ', pways)

  #format as corpus
  pways <- tm::Corpus(tm::VectorSource(pways))

  #remove unneeded words
  pways <- tm::tm_map(pways, tm::removeWords, excludeWords)

  #remove whitespaces?
  pways <- tm::tm_map(pways, tm::stripWhitespace)

  # wordcloud::wordcloud(pways
  #                      , scale=c(5,0.5)     # Set min and max scale
  #                      , max.words=100      # Set top n words
  #                      , random.order=FALSE # Words in decreasing freq
  #                      , rot.per=0.35       # % of vertical words
  #                      , use.r.layout=FALSE # Use C++ collision detection
  #                      , colors=brewer.pal(8, "Dark2"))

  wordcloud::wordcloud(pways, ...)

  recordPlot()

}



#' Declare significant clusters in the output of pathway analysis function.
#'
#' Set the signficant clusters in the pathwayanalysis_out object. Will return an updated version of pathwayanalysis_out[[3]] object with significance (T/F)
#'
#' @param pathwayanalysis_out the output of pathwayanalysis function
#' @param sigclusts vector of integers. the rows in pathwayanalysis_out[[3]] to use as significant clusters. If not provided, will declare clusters with > 10% of DEGs as significant.
#'
#' @return an updated pathwayanalysis_out object.
#' @export
#'
#' @examples
pathwayanalysis_declare_significant_clusters <- function(pathwayanalysis_out, sigclusts){

  clustpercs <- pathwayanalysis_out[[3]]

  if(missing(sigclusts)){
    #just check any above 10%?
    # get percents for the top clusters

    message('Setting significant clusters as clusters with > 10% of DEGs')

    clustpercs$significant <- F
    clustpercs[clustpercs$perc>=10, "significant"] <- T

  } else{

    message('Using manually set significant clusters')

    clustpercs$significant <- F
    clustpercs[sigclusts,"significant"] <- T
  }




  pathwayanalysis_out[[3]] <- clustpercs

  pathwayanalysis_out


}



#' Relabel some pathway clusters with biologically relevant names
#'
#' Pathways often share overlap in their genes, which are identified in this pipeline as clusters. Typically, pathways will group together as a "functional module" that can be identified based on literature and biological knowledge, so this function labels the clusters as such. This is usually kind of difficult and takes some time inspecting genes and pathway names/descriptions in each cluster.
#'
#' @param pathwayanalysis_out output from the `pathwayanalysis()` function.
#' @param clusterlabels character vector; labels to give to the clusters of pathways. If null, will use generic names.
#'
#' @return will return a list, similar in format to the output of `pathwayanalysis` function, but with updated object in element 1 and cluster labels in element 3.
#' @export
#'
#' @examples
pathwayanalysis_relabel_significant_clusters <- function(pathwayanalysis_out,
                                                         clusterlabels){




  msigdb_ora <- pathwayanalysis_out[[1]]
  clustpercs <- pathwayanalysis_out[[3]]
  reslist <- pathwayanalysis_out[[4]]
  clustgeneslist <- pathwayanalysis_out[[5]]

  #make sure clustpercs has significant clusters
  if( !('significant' %in% colnames(clustpercs)) ){

    stop('No significant clusters found, please run pathwayanalysis_declare_significant_clusters() function')
  }


  if(missing(clusterlabels) ){
    warning('No cluster labels provided;\nWill default to generic labels from significant clusters')
    clusterlabels <- clustpercs[clustpercs$significant==T,'clust']
  }

  #relabel significant with the names
  clustpercs$relabel <- clustpercs$clust
  clustpercs[clustpercs$significant==T,'relabel'] <- clusterlabels

  # add percent of degs to each cluster name
  clustpercs$relabel <- paste0(clustpercs$relabel, '\n', round(clustpercs$perc, 1), ' % of DEGs')




  #add to obj
  sigclusts <- clustpercs[clustpercs$significant == T,]

  # get pathways in the top clusters
  sigclustpways <- reslist[names(reslist) %in% sigclusts$clust]
  names(sigclustpways) <- sigclusts$relabel

  # get pathwyas in top clusters as a df
  sigclustpwaysdf <- bind_rows( lapply(1:length(sigclustpways), function(x){
    i=sigclustpways[[x]]
    name=names(sigclustpways)[x]
    data.frame(pways=i, clust=name)
  })
  )


  msigdb_ora@result$ClusterMain <- NA
  msigdb_ora@result[match(sigclustpwaysdf$pways, msigdb_ora@result$ID),'ClusterMain'] <- sigclustpwaysdf$clust


  pathwayanalysis_out[[1]] <- msigdb_ora
  pathwayanalysis_out[[3]] <- clustpercs


  pathwayanalysis_out

}







#' Make a Venn Diagram with genes from the significant clusters.
#'
#' Please run `pathwayanalysis_declare_significant_clusters` first.
#'
#'
#' @param pathwayanalysis_out the output of the pathwayanalysis function.
#'
#' @return
#' @export
#'
#' @examples
pathwayanalysis_venn <- function(pathwayanalysis_out){


  clustpercs <- pathwayanalysis_out[[3]]
  reslist <- pathwayanalysis_out[[4]]
  clustgeneslist <- pathwayanalysis_out[[5]]

  #make sure clustpercs has significant clusters
  if( !('significant' %in% colnames(clustpercs)) ){

    stop('No significant clusters found, please run pathwayanalysis_declare_significant_clusters() function')
  }

  #make sure clustpercs has significant clusters
  if( !('relabel' %in% colnames(clustpercs)) ){

    stop('Informative cluster names not found, please run pathwayanalysis_relabel_significant_clusters() function')
  }



  sigclusts <- clustpercs[clustpercs$significant == T,]

  # get actuaal genes of the top clusters
  sigclustgenes <- clustgeneslist[names(clustgeneslist) %in% sigclusts$clust]
  names(sigclustgenes) <- sigclusts$relabel

  # venn of significant clusters
  sig_venn <- ggVennDiagram::ggVennDiagram(sigclustgenes)+
    labs(title = 'Venn Diagram of genes from top clusters of enriched pathways')+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(expand = expansion(mult = .2))

  sig_venn


}




#' Add more emapplots (graph plots with edges between similar pathways) containing cluster info.
#'
#' @param pathwayanalysis_out output from the `pathwayanalysis()` function.
#' @param palette character vector of colors to use (hex codes)
#'
#' @return updated `pathwayanalysis_out` object with plots added to `pathwayanalysis_out[[2]]`
#' @export
#'
#' @examples
pathwayanalysis_finalize_emap_plots <- function(pathwayanalysis_out, palette){

  if(missing(palette)){palette <- RColorBrewer::brewer.pal(Inf, 'Accent')}

  msigdb_ora <- pathwayanalysis_out[[1]]
  plots <- pathwayanalysis_out[[2]]
  emap_total <- pathwayanalysis_out[[2]]$emap_noclusters
  clustpercs <- pathwayanalysis_out[[3]]
  reslist <- pathwayanalysis_out[[4]]
  clustgeneslist <- pathwayanalysis_out[[5]]

  #make sure clustpercs has significant clusters
  if( !('significant' %in% colnames(clustpercs)) ){

    stop('No significant clusters found, please run pathwayanalysis_declare_significant_clusters() function')
  }

  #make sure clustpercs has significant clusters
  if( !('relabel' %in% colnames(clustpercs)) ){

    stop('Informative cluster names not found, please run pathwayanalysis_relabel_significant_clusters() function')
  }



  sigclusts <- clustpercs[clustpercs$significant == T,]

  # get pathways in the top clusters
  sigclustpways <- reslist[names(reslist) %in% sigclusts$clust]
  names(sigclustpways) <- sigclusts$relabel

  # get pathwyas in top clusters as a df
  sigclustpwaysdf <- bind_rows( lapply(1:length(sigclustpways), function(x){
    i=sigclustpways[[x]]
    name=names(sigclustpways)[x]
    data.frame(pways=i, clust=name)
  })
  )


  #colors are defined above
  pal <- palette
  levs <- nrow(sigclusts)
  pal <- pal[1:levs]
  pal <- c(pal, 'grey50')

  #get coords from existing emap
  dat <- emap_total$data
  dat$Cluster <- NA
  dat[match(sigclustpwaysdf$pways, dat$name),'Cluster'] <- sigclustpwaysdf$clust


  repelaggr <- aggregate(x ~ Cluster, dat, mean)
  repelaggr$y <- aggregate(y ~ Cluster, dat, mean)[,2]
  repelaggr$Color <- pal[1:nrow(repelaggr)]


  emap_total_withclust <- emap_total +  ggrepel::geom_label_repel(inherit.aes = F,
                                                                  data = repelaggr, aes(x=x,y=y,label=Cluster, fill=Cluster),
                                                                  fill= repelaggr$Color,
                                                                  box.padding = 5, max.overlaps = 200,
                                                                  # direction = 'x',
                                                                  min.segment.length = Inf)


  rownames(dat) <- dat$name #rownames of coords must match names...
  emap_final <-  enrichplot::emapplot(msigdb_ora, color='ClusterMain', repel=F, showCategory = 200, # layout='graphopt',
                          node_label='None', #cex_label_category=0.4,
                          coords = dat[,1:2],
                          cex_category = 0.3,cex_line = 0.3 )+
    scale_fill_manual(values = pal, name = 'Cluster')



  # plots <- list(emap_noclusters = emap_total,
  #               emap_total_withclust = emap_total_withclust,
  #               emap_onlyclusters = emap_final)
  # pathwayanalysis_out[[2]] <- plots

  #instead of replacing plots list element, just append to it

  pathwayanalysis_out[[2]]['emap_total_withclust'] <- emap_total_withclust
  pathwayanalysis_out[[2]]['emap_onlyclusters'] <- emap_final

  pathwayanalysis_out


}





#' Make plots for pathways in each cluster.
#'
#' @param pathwayanalysis_out output from the `pathwayanalysis()` function.
#' @param sigres
#'
#' @return
#' @export
#'
#' @examples
pathwayanalysis_sigcluster_subplots <- function(pathwayanalysis_out, sigres){



  msigdb_ora <- pathwayanalysis_out[[1]]
  clustpercs <- pathwayanalysis_out[[3]]
  reslist <- pathwayanalysis_out[[4]]
  clustgeneslist <- pathwayanalysis_out[[5]]

  #make sure clustpercs has significant clusters
  if( !('significant' %in% colnames(clustpercs)) ){

    stop('No significant clusters found, please run pathwayanalysis_declare_significant_clusters() function')
  }

  #make sure clustpercs has significant clusters
  if( !('relabel' %in% colnames(clustpercs)) ){

    stop('Informative cluster names not found, please run pathwayanalysis_relabel_significant_clusters() function')
  }



  sigclusts <- clustpercs[clustpercs$significant == T,]

  #get the pwayres object from the msigdb dataframe
  pwayres <- msigdb_ora@result

  # get pathways in the top clusters
  sigclustpways <- reslist[names(reslist) %in% sigclusts$clust]
  names(sigclustpways) <- sigclusts$relabel


  namedfc <- sigres[,2]
  names(namedfc) <- sigres[,1]

  subclustlist <- list()
  for(sigclustidx in c(1:nrow(sigclusts)) ){

    origlabel <- sigclusts[sigclustidx, 'clust']
    relabel <- sigclusts[sigclustidx, 'relabel']

    message(relabel)

    #get the pathways
    sigpways <- sigclustpways[[sigclustidx]]

    #recompute...

    ### suba emap skip it for now, THIS IS BROKEN CURRENTLY
    # https://github.com/YuLab-SMU/clusterProfiler/issues/488

    # subemap <- emapplot(msigdb_ora, color='Percent_of_DEGs', repel=T,
    #                     #coords = subdat[,1:2],
    #                     showCategory = sigpways,
    #                     cex_label_category=0.4,
    #                     cex_category = 0.3,cex_line = 0.3
    # )+
    #   scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent of DEGs')+
    #   ggtitle(relabel) +
    #   theme(plot.title = element_text(hjust = 0.5))



    #selct top pathways from cluster to show...
    pwayres_inthisclust <- pwayres[pwayres$ID %in% sigpways,]

    #select top 3 based n num DEGs
    num_to_pick <- ifelse(nrow(pwayres_inthisclust) >= 3, 3, nrow(pwayres_inthisclust))
    pwayres_inthisclust <- pwayres_inthisclust[order(pwayres_inthisclust$Percent_of_DEGs, decreasing = T),]
    pways_to_show <- pwayres_inthisclust$ID[1:num_to_pick]

    subcnet <- enrichplot::cnetplot(msigdb_ora, showCategory = pways_to_show, foldChange = namedfc,
                        shadowtext='category',
                        cex_gene = 0.5, cex_label_gene = 0.5, cex_label_category = 0.7)+
      ggtitle(relabel) +
      theme(plot.title = element_text(hjust = 0.5))




    # subora@result$qvalue <- subora@result$Percent_of_DEGs
    subddotplot <-  enrichplot::dotplot(msigdb_ora, showCategory=head(pwayres_inthisclust$ID,30),font.size=6)+
      ggtitle(relabel) +
      theme(plot.title = element_text(hjust = 0.5))



    wc <- pathwayanalysis_wordcloud(unlist(reslist),
                                    scale=c(3,0.3), random.order=F, random.color=F,
                                    colors= RColorBrewer::brewer.pal(8, "Dark2"))
    title('Overall Pathways Wordcloud')
    wc <- recordPlot()



    # pdf(paste0(clustoutdir, '/subplots_', origlabel, '.pdf'), width = 8, height = 8)
    # print(subemap)
    # print(subcnet)
    # print(subddotplot)
    #
    # dev.off()


    subclustlist[[sigclustidx]] <- list(
                                        #subemap,
                                        subcnet,
                                        subddotplot,
                                        wc)

  }



  pathwayanalysis_out$cluster_subplots <- subclustlist


  pathwayanalysis_out

}








#' Unlist2 from AnnotationDBI package
#'
#' See `?AnnotationDbi::unlist2()` by Hervé Pagès
#'
#' @param x
#' @param recursive
#' @param use.names
#' @param what.names
#'
#' @return
#' @export
#'
#' @examples
unlist2 <- function (x, recursive = TRUE, use.names = TRUE, what.names = "inherited")
{
  ans <- unlist(x, recursive, FALSE)
  if (!use.names)
    return(ans)
  if (!is.character(what.names) || length(what.names) != 1)
    stop("'what.names' must be a single string")
  what.names <- match.arg(what.names, c("inherited", "full"))
  names(ans) <- unlist(make.name.tree(x, recursive, what.names),
                       recursive, FALSE)
  ans
}
