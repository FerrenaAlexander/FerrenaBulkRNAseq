
### to do items

# finish off:
# left off from venn diagram
# make the wordcloud save as a list object? ie, https://stackoverflow.com/questions/29583849/save-a-plot-in-an-object



#### pathway analysis, clusterprofiler wrapper ####

#' Run Overrepresentation Analysis via the ClusterProfiler workflow
#'
#' @param sigres data.frame of significant genes. first column is gene names; second column is logfoldchange or other effect size; third column is Pvalue (adusted is best)
#' @param term2gene data.frame of pathways to check. first column is pathways, second column is genes.
#'
#' @return a list object: first element is a ClusterProfiler object; second element is an "emap plot", a graph of pathways connected by jaccard similarity; third element is a data.frame with the connected nodes from the emap plot; 4th element is a list of pathways in each connected cluster.
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
  msigdb_ora <- enrichplot::pairwise_termsim(msigdb_ora)

  #try to plot color = num genes in overexp gene list?
  msigdb_ora@result$Percent_of_DEGs <- (msigdb_ora@result$Count / nrow(sigres))*100


  set.seed(2021)
  options(ggrepel.max.overlaps = 2000)

  suppressMessages(
    emap_total <- emapplot(msigdb_ora, color='Percent_of_DEGs', repel=T, showCategory = 200,
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


  output <- list(msigdb_ora, emap_total, clustpercs, reslist)

}



#' A fast way to make a wordcloud from a character vector of pathway names
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
wordcloud_pathways <- function(pways, delimiter, excludeWords, ...){

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




#' Relabel some pathway clusters with biologically relevant names
#'
#' @param pathwayanalysis_out output of the `pathwayanalysis` function.
#' @param clusterlabels character vector; labels to give to the clusters of pathways
#' @param clusters_to_label vector of integers; which clusters to re-name.
#'
#' @return will return a list, similar in format to the output of `pathwayanalysis` function, but with updated cluster labels in element 3.
#' @export
#'
#' @examples
pathwayanalysis_relabel_significant_clusters <- function(pathwayanalysis_out,
                                                         clusterlabels,
                                                         clusters_to_label){



  clustpercs <- pathwayanalysis_out[[3]]

  clustpercs$relabel <- clustpercs$clust
  clustpercs[clusters_to_label,'relabel'] <- clusterlabels
  clustpercs$relabel <- paste0(clustpercs$relabel, '\n', round(clustpercs$perc, 1), ' % of DEGs')

  pathwayanalysis_out[[3]] <- clustpercs


  pathwayanalysis_out

}

