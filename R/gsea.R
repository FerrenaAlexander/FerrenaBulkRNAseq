

### gsea functions ###


### gsea results and leading edge
#' GSEA results.
#'
#' Using FGSEA, DESeq2 results, and input pathways, compute GSEA results.
#'
#' @param deg dataframe with column names similar to the output of DESeq2 "results()" function
#' @param pathways list of character vectors containing gene names in same format as deg; each list element's name should be the name of the pathway
#' @param nperm how many permutations to run in FGSEA (may be deprecated). Default is 10000
#' @param weightmethod one of "pvalue" or "foldchange". how to weight the DEG results. Default is "Pvalue"
#' @param onlypos T/F. whether to use only positive log fold change genes. Default is False.
#'
#' @return a list with two elements; 1, a dataframe with FGSEA results; 2, a list of character vectors showing the leading edge genes for each pathway
#' @export
#'
#' @examples
gsea.results <- function(deg,
                         pathways=NULL,
                         nperm=NULL,
                         weightmethod=NULL,
                         onlypos=NULL
){


  if(is.null( pathways )) {pathways <- tamlabscpipeline::hallmark}
  if(is.null( nperm )) {nperm <- 10000}
  if(is.null( weightmethod )) {weightmethod <- 'pvalue'}
  if(is.null( onlypos )) {onlypos <- F}



  tmp <- deg

  if(onlypos == T){
    tmp[tmp$log2FoldChange >0,]
  }

  if(weightmethod == 'pvalue') {
    scores <- log(tmp$pvalue)

    #fix underflow
    if(any(scores == -Inf)){

      numuf = length(scores[scores==-Inf])
      adder = rev(1:numuf)
      for(uf in rev(1:numuf)){
        scores[uf] <- scores[numuf + 1] + (adder[uf] * -1)
      }

    }

    scores <- (-1 * scores) * sign(tmp$log2FoldChange)
    names(scores) <- tmp$Gene_name
    rm(tmp)
    scores <- sort(scores, decreasing = T)
  }

  if(weightmethod == 'foldchange'){
    scores <- 10^tmp$log2FoldChange
    names(scores) <- tmp$Gene_name
    rm(tmp)
    scores <- sort(scores, decreasing = T)
  }

  if(weightmethod != 'pvalue' & weightmethod != 'foldchange'){
    stop('weightmethod must be one of either "pvalue" or "foldchange"; default (null) is pvalue')
  }


  suppressWarnings( fgseaRes <- fgsea(pathways=pathways, stats=scores, nperm=nperm) )
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    as.data.frame()


  res <- fgseaResTidy

  rm(fgseaRes)
  rm(fgseaResTidy)
  rm(scores)




  #got some NAs at some point in the past, not surw why
  res$NES[is.na(res$NES)] <- 0

  #deal with leading edge
  leadingedge <- res$leadingEdge
  res <- res[,-8]

  names(leadingedge) <- res$pathway

  #return list: first element iis GSEA results DF, second element of leading edge genes
  list(res, leadingedge)

}


### gsea dotplot
gsea.dotplot.onecol <- function(deg,
                                pathways=NULL,
                                nperm=NULL,
                                weightmethod=NULL,
                                onlypos=NULL,
                                filter_nonsig_pathways=NULL,
                                top_pathways = NULL,
                                ntop = NULL,
                                fix_long_pathway_names = NULL,
                                pathwayfontsize=NULL,
                                gsubpattern=NULL
){


  if(is.null( pathways )) {pathways <- tamlabscpipeline::hallmark}
  if(is.null( nperm )) {nperm <- 10000}
  if(is.null( weightmethod )) {weightmethod <- 'pvalue'}
  if(is.null( onlypos )) {onlypos <- F}
  if(is.null( filter_nonsig_pathways )) {filter_nonsig_pathways <- F}
  if(is.null( top_pathways  )) {top_pathways <- T}
  if(is.null( ntop  )) {ntop <- 25}
  if(is.null( fix_long_pathway_names  )) {fix_long_pathway_names <- T}
  if(is.null( pathwayfontsize  )) {pathwayfontsize <- 10}
  if(is.null( gsubpattern  )) {gsubpattern <- 'HALLMARK_'}



  tmp <- deg

  if(onlypos == T){
    tmp[tmp$log2FoldChange >0,]
  }

  if(weightmethod == 'pvalue') {
    scores <- log(tmp$pvalue)

    #fix underflow
    if(any(scores == -Inf)){

      numuf = length(scores[scores==-Inf])
      adder = rev(1:numuf)
      for(uf in rev(1:numuf)){
        scores[uf] <- scores[numuf + 1] + (adder[uf] * -1)
      }

    }

    scores <- (-1 * scores) * sign(tmp$log2FoldChange)
    names(scores) <- tmp$Gene_name
    rm(tmp)
    scores <- sort(scores, decreasing = T)
  }

  if(weightmethod == 'foldchange'){
    scores <- 10^tmp$log2FoldChange
    names(scores) <- tmp$Gene_name
    rm(tmp)
    scores <- sort(scores, decreasing = T)
  }

  if(weightmethod != 'pvalue' & weightmethod != 'foldchange'){
    stop('weightmethod must be one of either "pvalue" or "foldchange"; default (null) is pvalue')
  }


  suppressWarnings( fgseaRes <- fgsea(pathways=pathways, stats=scores, nperm=nperm) )
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    as.data.frame()


  res <- fgseaResTidy

  rm(fgseaRes)
  rm(fgseaResTidy)
  rm(scores)



  #remove repetitive pattern
  if(!is.null(gsubpattern)){
    res$pathway <- gsub(gsubpattern, '', res$pathway)
  }

  #got some NAs at some point in the past, not surw why
  res$NES[is.na(res$NES)] <- 0

  #get rid of leading edge...
  res <- res[,-8]


  #fix long names
  if(fix_long_pathway_names == T){
    pnames <- res$pathway
    pnewnames <- c()
    for(p in pnames){


      if(str_length(p) > 50){

        #try to find and replace underscore closest to halfway...
        halfway <- round(str_length(p)/2)
        underscorepos <- str_locate_all(p,'_')[[1]][,1]

        distances <- abs(halfway-underscorepos)

        which_und_is_closest <- which.min(distances)

        split <- str_split_fixed(p,'_',Inf)

        newp <- paste0(paste(split[1,1:which_und_is_closest], collapse = '_'),
                       '        ', '\n',
                       paste(split[1,(which_und_is_closest+1):ncol(split)], collapse = '_')
        )



      } else{newp <- p}

      pnewnames <- c(pnewnames, newp)

    }

    res$pathway <- pnewnames
  }


  #level sorting: find the highest/lowest pathways; l2FC; only 2-way...
  res <- res[order(res$NES, decreasing = F),]
  pwaylevs <- res[order(res$NES, decreasing = F),'pathway']




  #for plotting... only keep top N (ntop) from each from each.
  #if yes, adjust pathways if ntop < total num;
  # first check if all sign NES are equal (rare), if not then take top and bottom

  if(top_pathways == T){

    if(ntop < nrow(res) ){

      if( all(sign(res$NES) == sign(res$NES)[1]) ){
        pwaylevs <- tail(pwaylevs, ntop)

      } else {


        pwaylevs <- unique(
          c(
            head(pwaylevs, ntop ),
            tail(pwaylevs, ntop)
          )
        )
      }
    }
  }

  res <- res[res$pathway %in% pwaylevs,]



  #levels set levels as defined above
  res$pathway <- factor(res$pathway, levels = pwaylevs)


  #to remove nonsignificant rows...
  if(filter_nonsig_pathways == T){
    res2 <-  data.frame()
    for(set in unique(res$pathway)){
      res_tmp <- res[res$pathway == set,]
      if(any(res_tmp$padj < 0.25)){
        res2 <- rbind(res2, res_tmp)
      } else {NULL}
    }


    if( nrow(res2) > 0 ) {
      res <- res2
      rm(res2) } else { message('No significant results for ', j, ' detected.') }

  }

  #plot

  #significance marks
  resx <- res[,c('pathway', 'padj')]

  levs <- levels(resx$pathway)
  resx$pathway <- as.character(resx$pathway)
  resx <- resx[match(levs, resx$pathway),]

  resx$pathway <- ifelse(resx$padj < 0.05,
                         yes = paste0(resx$pathway, '    ***'),
                         no = resx$pathway)


  resx$pathway <- ifelse(resx$padj >= 0.05 & resx$padj < 0.25,
                         yes = paste0(resx$pathway, '     * '),
                         no = resx$pathway)


  resx$pathway <- ifelse(resx$padj >= 0.25,
                         yes = paste0(resx$pathway, '        '),
                         no = resx$pathway)


  res$pathway <- plyr::mapvalues(x = res$pathway,
                                 from = levels(res$pathway),
                                 to = resx$pathway)



  gseaout <- ggplot(res, aes(x='', y = pathway , col = NES, size = size) )+
    geom_point()+
    scale_color_distiller(palette = 'RdBu', name="Normalized\nEnrichment\nScore")+
    theme_minimal()+
    theme(legend.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = pathwayfontsize, face='bold'),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
    )+
    guides(size=guide_legend(title="Number\nof Genes"))+
    labs(caption = paste0('*** padj < 0.05\n',
                          ' *  padj < 0.25')
    )

  gseaout


}



