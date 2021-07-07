# https://r-pkgs.org/whole-game.html



#' Make a VolcanoPlot
#'
#' Using DESeq2 results, make a volcanoplot
#'
#' color can be set with colors or passed using metadata. the latter is encouraged, in keeping with setting a color scheme at the start and keeping it consistent in all analysis.
#'
#' @param results dataframe, with column names similar to the output of DESeq2::results() function
#' @param metadata dataframe, for coloring samples using a consistently set color scheme. metadata is a dataframe similar to coldata used in DESEQ2; include a contrast column and colors column
#' @param contrast character string. name of contrast column in metadata. default will search for "Conditions"
#' @param condition1 character string. will be used to find the corresponding contrast-color pairing for upregulated genes
#' @param condition2 character string. will be used to find the corresponding contrast-color pairing for downregulated genes
#' @param colors character vector of length two with string names or codes of colors. if not using metadata-based coloring, will use this input or default to red/blue.
#' @param pval_thres numeric, alpha value cutoff for significant genes. default = 0.05
#' @param lfc_thres numeric, logFC cutoff for significant genes, will be symmetric on both sides. default = 0
#' @param lfclimit numeric, lofFC cutoff for outlier genes, will be symmetric on both sides. default = 15
#' @param pointsize numeric, point size of graph, passed to geom_point, default = 0.5
#' @param outliergene_pointsize numeric, point size of logFC outlier genes, default = 2
#' @param outliergene_shape numeric, code for shape symbol of logFC outlier genes, default = 6, which is a triangle
#' @param repeltextsize numeric, size of repel text gene labels, default is 2
#'
#' @return a ggplot2 object.
#' @export
#'
#' @examples
#'
#'
#' With default colors:
#' volcanoplot(results)
#'
#' With set colors:
#' volcanoplot(results, colors = c('Purple', 'Pink'))
#'
#'
#' To use metadata for coloration, we need a metadata dataframe
#'  that looks like this:
#'
#' Samples   Condition    Color
#' Sample1   KO           Red
#' Sample2   KO           Red
#' Sample3   WT           Blue
#' Sample4   WT           Blue
#'
#' contrast in this case can be input as contrast = "Condition"
#' metadata needs a column called Color, set it at the outset,
#' it will make all your plotting easier!
#' condition1 and condition2 should be set,
#' ie condition1 = "KO" and condition2 = "WT"
#'
#' If we have all that, we can do this:
#' volcanoplot(results, metadata = metadata,
#'             condition1 = 'KO', condition2 = 'WT')
#'
#'
volcanoplot <- function(results,
                        metadata,
                        contrast,
                        condition1, condition2,
                        colors,
                        pval_thres,
                        lfc_thres,
                        lfclimit,

                        pointsize,
                        outliergene_pointsize,
                        outliergene_shape,
                        repeltextsize
                        ){

  if( missing(results) ){stop('Provide results df in format of DESeq2::results() output')}
  if( missing(metadata) ){message('No metadata provided, will rever to provided or default colors')}
  if( missing(contrast) ){contrast <- 'Condition'}
  if( missing(pval_thres) ) {pval_thres <- 0.05}
  if( missing(lfc_thres) ) {lfc_thres <- 0}
  if( missing(lfclimit) ) {lfclimit <- 15}

  if( missing(pointsize) ) {pointsize <- 0.5}
  if( missing(outliergene_pointsize) ) {outliergene_pointsize <- 2}
  if( missing(outliergene_shape) ) {outliergene_shape <- 6}
  if( missing(repeltextsize) ) {repeltextsize <- 2}



  #set as DF and order by Padj
  restmp <- as.data.frame(results)

  #order by pvalue
  restmp <- restmp[order(restmp$pvalue),]

  #remove outlier genes with p = NA (cooks distance)
  restmp <- restmp[!is.na(restmp$pvalue),]
  restmp <- restmp[!is.na(restmp$padj),]

  #make df that's eaiser to deal with:
  colskeep <- c('baseMean', 'log2FoldChange',
                'pvalue', 'padj')
  restmp <- restmp[,c(colskeep)]

  #set up gene_name col
  restmp$Gene_name <- rownames(restmp)

  #set up repel label by applying cutoffs ; this essentially selects the significant genes
  # significance is defined by lfc_thres and pval_thres
  restmp$repel <- NA
  restmp[abs(restmp$log2FoldChange) > lfc_thres &restmp$pvalue < pval_thres, 'repel'] <-
    restmp[abs(restmp$log2FoldChange) > lfc_thres & restmp$pvalue < pval_thres, "Gene_name"]


  #check if colors are provided in metadata. If not, use provided or preset red/blue
  if(missing(metadata) ){

    if(missing(colors)){
      colors <- c("#E41A1C", "#377EB8") #red and blue preset
    }

    #color repel by condition
    restmp$repelcol <- 'grey'
    restmp[!is.na(restmp$repel) & restmp$log2FoldChange > 0,'repelcol'] <- colors[1]
    restmp[!is.na(restmp$repel) & restmp$log2FoldChange < 0,'repelcol'] <- colors[2]



  } else{

    #check if color column is set
    colorcol <- grep('color', colnames(metadata), ignore.case = T)
    if(length(colorcol) != 1) {stop ('Metadata provided but no unambigious "Color" column found. Make sure to have a single column labelled "Color"')}

    #check if contrast column is set
    contrastcol <- grep(contrast, colnames(metadata), ignore.case = T)
    if(length(contrastcol) != 1) {stop ('Metadata provided but no unambigious contrast column found. Make sure to set contrast.')}

    #check if condition labels are provided
    if(missing(condition1) | missing(condition2)){ stop('Metadata provided but contrast conditions are not labelled. Make sure to set Condition1 and Condition2')}

    #use predefined colors in metadata to call colors
    restmp$repelcol <- 'grey'
    restmp[!is.na(restmp$repel) & restmp$log2FoldChange > 0,'repelcol'] <- metadata[metadata[,contrastcol] == condition1,colorcol][1]
    restmp[!is.na(restmp$repel) & restmp$log2FoldChange < 0,'repelcol'] <- metadata[metadata[,contrastcol] == condition2,colorcol][1]

  }

  #make genes that don't pass thres grey color
  restmp$gencolscheme <- 'grey'
  #restmp[!is.na(restmp$repel),'gencolscheme'] <- 'black'


  #get neg log p
  restmp$neglogtenp <- -log10(restmp$pvalue)


  #color the significant DEGs
  # these two lines are equivalent at this point, but the second is more readable
  #significantres <- restmp[!is.na(restmp$repel),]
  significantres <- restmp[abs(restmp$log2FoldChange) > lfc_thres & restmp$pvalue < pval_thres,]


  ### repel labelling ###
  #if any extreme genes exist that fall out of LFC limit, we mark and label them after
  if( any(abs(significantres$log2FoldChange) > lfclimit) ){

    ext <- restmp[abs(restmp$log2FoldChange) > lfclimit,]
    ext <- ext[ext$pvalue < pval_thres,]

    ext$log2FoldChange <- ifelse(sign(ext$log2FoldChange) == -1,
                                 yes = -lfclimit, no = lfclimit)

  }

  #get significant DEGs which do not fall outside of LFC boundary, label 50 per side
  keeprepel <- significantres[abs(restmp$log2FoldChange) <= lfclimit,]

  #for repel labelling: get the first get top 25 pval genes per LFC side
  keeprepel <- rbind( head ( keeprepel[keeprepel$log2FoldChange > lfc_thres & keeprepel$log2FoldChange < lfclimit, ], 25),
                      head ( keeprepel[keeprepel$log2FoldChange < -lfc_thres & keeprepel$log2FoldChange > -lfclimit, ], 25)
  )

  #next get top 25 LFC genes per LFC side
  # after removing genes so as not to have double rows
  significantres <- significantres[!(significantres$Gene_name %in% keeprepel$Gene_name),]
  significantres <- significantres[order(abs(significantres$log2FoldChange), decreasing = T),]

  keeprepel <- rbind(keeprepel,
                     head ( significantres[significantres$log2FoldChange > lfc_thres, ], 25),
                     head ( significantres[significantres$log2FoldChange < -lfc_thres, ], 25)
  )







  #if any repel genes are out of bounds, make sure they are labelled
  if( any(abs(significantres$log2FoldChange) > lfclimit) ){

    keeprepel <- rbind(keeprepel, ext)

  }




  #volcano plot

  suppressWarnings(
    vesuvius <- ggplot(restmp, aes(x = log2FoldChange, y = neglogtenp))+
      geom_point(col = restmp$repelcol, size = pointsize)+
      geom_point(data= ext, col = ext$repelcol, size = outliergene_pointsize, shape = outliergene_shape)+
      geom_vline(xintercept = c(lfc_thres, -1*lfc_thres), linetype = 'dotted', col = 'firebrick', alpha=0.7)+
      geom_hline(yintercept = -log10(pval_thres),  linetype = 'dotted', col = 'firebrick', alpha = 0.7)+
      ggrepel::geom_text_repel(data = keeprepel, aes(label = repel), col = keeprepel$repelcol, size = repeltextsize)+
      theme_light()+
      scale_x_continuous(breaks = seq(-lfclimit,lfclimit,by=5), limits = c(-lfclimit,lfclimit))+
      labs(y = '-Log10 P-value',
           x = paste0('Log2 Fold Change')
      )

  )

 vesuvius


}
