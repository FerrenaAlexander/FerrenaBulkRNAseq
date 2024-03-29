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
#' @param use_padj T/f, default = T. Use adjustec p-value (passed as df[,'padj']) for y axis and coloring of significant genes.
#' @param pointsize numeric, point size of graph, passed to geom_point, default = 0.5
#' @param outliergene_pointsize numeric, point size of logFC outlier genes, default = 2
#' @param outliergene_shape numeric, code for shape symbol of logFC outlier genes, default = 6, which is a triangle
#' @param repeltextsize numeric, size of repel text gene labels, default is 2
#' @param change_gene_label T/F, whether to change gene labels, defautl F
#' @param gene_label_equivalency data.frame, if change_gene_label is set to T, need to provide a data.frame with two columns, first column with current gene labels (rownames of res) and second column with gene labels you want to plot as labels. Useful for when you have res with IDs but want to show gene symbols.
#'
#' @return a ggplot2 object.
#' @export
#'
#' @examples
#'
#'
#' With default colors (red up, blue down):
#' volcanoplot(results)
#'
#' With set colors:
#' volcanoplot(results, colors = c('Purple', 'Pink'))
#'
#' With different labels, for example if your results have EnsemblIDs but you want to show gene symbols:
#' volcanoplot(results, change_gene_label = T, gene_label_equivalency = gene_name_df)
#'
#' where gene_name_df is a data.frame with first two columns like this:
#' GeneLabelsNum  Gene SymbolsLetter
#' 1              A
#' 2              B
#' 3              C
#'
#' Your results gene names should be like the first column and the labelled points will be like the second column.
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

                        use_padj,
                        pointsize,
                        outliergene_pointsize,
                        outliergene_shape,
                        repeltextsize,
                        change_gene_label,
                        gene_label_equivalency
){

  if( missing(results) ){stop('Provide results df in format of DESeq2::results() output')}
  if( missing(metadata) ){message('No metadata provided, will rever to provided or default colors')}
  if( missing(contrast) ){contrast <- 'Condition'}
  if( missing(pval_thres) ) {pval_thres <- 0.05}
  if( missing(lfc_thres) ) {lfc_thres <- 0}
  if( missing(lfclimit) ) {lfclimit <- 15}

  if( missing(use_padj) ) {use_padj <- T}
  if( missing(pointsize) ) {pointsize <- 0.5}
  if( missing(outliergene_pointsize) ) {outliergene_pointsize <- 2}
  if( missing(outliergene_shape) ) {outliergene_shape <- 6}
  if( missing(repeltextsize) ) {repeltextsize <- 2}

  if( missing(change_gene_label) ) {change_gene_label <- F}
  if( change_gene_label ==T & missing(gene_label_equivalency) ) { stop('If trying to change gene labels, please provide a data.frame mapping current labels with desired labels')}



  #set as DF and order by Pvalue
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


  # switch to default using padj...
  if(use_padj==T){
    restmp$pvalue <- restmp$padj
  }


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


  if(change_gene_label==T){
    gene_label_equivalency <- gene_label_equivalency[match(keeprepel$repel, gene_label_equivalency[,1]),]
    keeprepel$repel <- gene_label_equivalency[,2]
  }


  #volcano plot

  suppressWarnings(
    vesuvius <- ggplot(restmp, aes(x = log2FoldChange, y = neglogtenp))+
      geom_point(col = restmp$repelcol, size = pointsize)+
      geom_vline(xintercept = c(lfc_thres, -1*lfc_thres), linetype = 'dotted', col = 'firebrick', alpha=0.7)+
      geom_hline(yintercept = -log10(pval_thres),  linetype = 'dotted', col = 'firebrick', alpha = 0.7)+
      ggrepel::geom_text_repel(data = keeprepel, aes(label = repel), col = keeprepel$repelcol, size = repeltextsize)+
      theme_light()+
      scale_x_continuous(breaks = seq(-lfclimit,lfclimit,by=5), limits = c(-lfclimit,lfclimit))+
      labs(y = '-Log10 P-value',
           x = paste0('Log2 Fold Change')
      )

  )

  #add extreme genes if they exist
  if( any(abs(significantres$log2FoldChange) > lfclimit) ){

    vesuvius <- vesuvius + geom_point(data= ext, col = ext$repelcol, size = outliergene_pointsize, shape = outliergene_shape)

  }



  if(use_padj==T){
    vesuvius <- vesuvius+labs(y="-Log10 Adjusted P-value")
  }

  vesuvius

}










#' Make an expression heatmap with hierarchical clustering
#'
#' assumes a normalized counts matrix, should work with TPM too. Currently, color depends on addding metadata. Later, will add support for manually defined colors.
#'
#'
#'
#'
#' @param expmatrix matrix or data.frame with rows = genes and columns = samples
#' @param genes character vector of gene names (rownames) in the expmatrix
#' @param metadata data.frame in DESeq2 "coldata" format. Needs a "Sample" column with sample IDs matching colnames of expmat, a "Condition" column with the associated conditon, and a "Color" column that has the sample colors
#' @param heatmaptitle string, name of heatmap on top
#' @param legendtitle string, what to call the values in the gene expression matrix. Will default to "Zscaled Normalized Counts"
#' @param do.log2 T/F, whether to log2(x+1) transform gene expression values (rows) in the matrix
#' @param do.scale T/F, whether to scale (z-transform, (x - mean(x) / sd(x))) the gene expression values (rows) of the matrix
#' @param change_gene_label T/F, whether to change gene labels, defautl F
#' @param gene_label_equivalency data.frame, if change_gene_label is set to T, need to provide a data.frame with two columns, first column with current gene labels (rownames of res) and second column with gene labels you want to plot as labels. Useful for when you have res with IDs but want to show gene symbols.
#' @param clustering_distance_columns distance to use for ComplexHeatmap::Heatmap clustering of rows, default is pearson
#' @param ... Additionaly arguments to pass to ComplexHeatmap::Heatmap
#'
#'
#' @return
#' @export
#'
#' @examples
heatmapplot <- function(expmatrix,
                        genes,
                        metadata,

                        heatmaptitle,
                        legendtitle,

                        do.log2,
                        do.scale,

                        change_gene_label,
                        gene_label_equivalency,

                        clustering_distance_columns,

                        deres,
                        pval_thres,
                        lfc_thres,
                        use_padj,

                        ...

){

  if( missing(expmatrix) ){stop("Please provide Expression matrix")}
  if( missing(genes) ){stop("Please provide input genes")}
  if( missing(metadata) ){stop("Please provide metadata dataframe")}

  if( missing(legendtitle) ){legendtitle <- 'Expression'}
  if( missing(heatmaptitle) ){heatmaptitle <- ''}

  if( missing(do.log2) ){do.log2 <- T}
  if( missing(do.scale) ){do.scale <- T}

  if( missing(change_gene_label) ) {change_gene_label <- F}
  if( change_gene_label ==T & missing(gene_label_equivalency) ) { stop('If trying to change gene labels, please provide a data.frame mapping current labels with desired labels')}

  if( missing(clustering_distance_columns) ){ clustering_distance_columns <- 'pearson' }

  if( missing(pval_thres) ){pval_thres <- 0.05}
  if( missing(lfc_thres) ){lfc_thres <- 1}
  if( missing(use_padj) ){use_padj <- T}

  #subset geneex pmatrix
  tmpgem <- expmatrix[rownames(expmatrix) %in% genes,,drop=F]

  # LATER EDIT - adding drop=F should fix that...
  # #if only one gene, above forces tmpgem to be vector, but it needs to stay a one-row matrix
  # goodgenes <- rownames(expmatrix)[rownames(expmatrix) %in% genes]
  # if( length( goodgenes ) == 0 ){stop('no input genes in matrix :(')}
  # if( length( goodgenes ) == 1 ){
  #   tmpgem <- as.matrix(t(tmpgem))
  #   rownames(tmpgem) <- goodgenes
  # }

  #some genes have all 0s, remove those...
  if(any(rowSums(tmpgem) == 0)){
    zeros <- names(which(rowSums(tmpgem) == 0))
    warning("removing genes with all zeros:\n",
            paste(zeros, collapse = ', '))

    tmpgem <- tmpgem[!(rownames(tmpgem) %in% zeros),]

  }

  #log2 transform rows (gene exp  values)
  if(do.log2==T){
    tmpgem <- t( apply(tmpgem, 1, function(x){log2(x+1)}) )
  }

  #scale the rows; ie, scale gene expression value for each gene across samples
  if(do.scale==T){
    tmpgem <- t(scale(t(tmpgem)))
  }

  #set up the column annotation
  #annotation df, has sample (column) name and color
  annotdf <- data.frame(sample = colnames(tmpgem),
                        sampcheck = metadata[match(colnames(tmpgem), metadata$Sample), "Sample"],
                        cond = metadata[match(colnames(tmpgem), metadata$Sample), "Condition"],
                        color = metadata[match(colnames(tmpgem), metadata$Sample), "Color"])

  #define colors, need to use this list thing for complexheatmap
  hacol <- list(Condition = annotdf[,4]) ; names(hacol[[1]]) <- annotdf[,3]

  #create the annotation object
  ha <- ComplexHeatmap::HeatmapAnnotation(Condition = annotdf[,3], col = hacol, name = 'Condition')



  # if deres is given, then try to add asterisks to sig DEGs
  if(!missing(deres)){
    message('adding significance values')

    deres <- deres[complete.cases(deres,)]

    if(use_padj==T){
      sig <- deres[,deres$padj < 0.05]
    } else{
      sig <- deres[,deres$pvalue < 0.05]
    }

    sig <- sig[abs(sig$log2FoldChange) > lfc_thres,]

  sig <- sig[rownames(sig) %in% rownames(tmpgem)]

  rownames( tmpgem[match(rownames(sig), rownames(tmpgem)), ] ) <- rownames(sig)


  }





  #rowlabs, if less than 75 genes then label them, if not leave blank
  if (nrow(tmpgem) >75){  row_labels =  (rep('', nrow(tmpgem))) }else{ row_labels <- rownames(tmpgem)}

  #optionally chnage gene labels
  if( change_gene_label==T & length(row_labels) <= 75 ){

    gene_label_equivalency <- gene_label_equivalency[match(row_labels, gene_label_equivalency[,1]),]

    row_labels <- gene_label_equivalency[,2]

  }



  #heatmap
  hm <- ComplexHeatmap::Heatmap(tmpgem,
                                name = legendtitle,
                                cluster_rows = T, cluster_columns = T,
                                col = rev(colorRampPalette( brewer.pal(9, "RdYlBu") )(255)) ,
                                row_labels = row_labels,
                                row_names_gp = gpar(fontsize = 7.5),
                                column_names_gp = gpar(fontsize = 7.5), column_names_rot = 45,
                                column_title = heatmaptitle,
                                top_annotation = ha,
                                ...)

  hm

}








#' Marker annotated heatmap
#'
#' Heatmaps with annotations for samples (columns) and rows (genes), useful to plot celltype markers
#'
#' @param res a data.frame. one row for each gene. first three columns must be 1) gene symbol or id; 2) fold change; and 3) pvalue or padj. rest of the columns are normalized counts.
#' @param geneannots a data.frame. one row for each gene. first column must be gene, second column must be celltype, third column must be color for each celltype
#' @param metadata a data.frame. sample metadata, each row s a sample. must have rownames = colnames of samples in res (ie, colnames of `res[,-c(1:3)]` ); must also have a column called "Condition" for each sample condition, and a column called "Color" for each sample color; if two conditions, must have two colors
#' @param lfc_thres numeric, a cutoff for lfc to be considered significant, default is 1
#' @param pval_thres numeric, a cutoff for pval/padj to be considered significant, default is 0.05
#' @param do.scale T/F, whether to scale genes, default is T
#' @param do.log T/F, whether to log1p genes, default is T
#' @param legend_caption string, caption to use above figure legend, by default will guess based on scale/log options
#' @param drop.dup.genes T/F whether to drop duplicated genes, by default T, if set to F and dup genes are detected will throw an error
#' @param ... additional options passed on to `ComplexHeatmap::Heatmap()`, see `?ComplexHeatmap::Heatmap()`
#'
#' @return
#' @export
#'
#' @examples
markerheatmap <- function(res,
                          geneannots,
                          metadata,
                          lfc_thres,
                          pval_thres,
                          do.scale,
                          do.log,
                          legend_caption,
                          drop.dup.genes,
                          ...){

  require(ComplexHeatmap)

  if( missing(lfc_thres) ){lfc_thres <- 1}
  if( missing(pval_thres) ){pval_thres <- 0.05}
  if( missing(do.scale) ){do.scale <- T}
  if( missing(do.log) ){do.log <- T}

  if( missing(drop.dup.genes) ){drop.dup.genes <- T}


  ## parse inputs ##
  #genes, first column of gene annots
  genes <- geneannots[,1]

  #if no genes are in the res first column, stop here
  if( all( !(genes %in% res[,1]) ) ){stop('No input genes found in first column of res')}

  #if any genes are missing, warn and move on
  if( any( !(genes %in% res[,1]) ) ){
    missinggenes <- genes[!(genes %in% res[,1])]
    missinggenewarn <- paste0("Some input genes are missing from the matrix, including:\n",
                              paste(missinggenes, collapse = ', '))
    warning(missinggenewarn)
  }

  #if any genes are duplicated, error
  if( any( duplicated(genes) )  ){
    dup_input_genes <- genes[duplicated(genes)]
    dup_input_msg <- paste0("Some input genes duplicated, including:\n",
                            paste(dup_input_genes, collapse = ', '))
    stop(dup_input_msg)
  }



  #make sure there are no duplicate genes:
  dupgenes_gem <- res[,1]
  dupgenes_gem <- dupgenes_gem[table(dupgenes_gem)>1]

  # if the duplicates are in the gene list, raise an error,
  # if not, just drop them from res
  if( any(dupgenes_gem %in% genes) ){
    dupgenes_gem_IN_INPUT_GENES <- dupgenes_gem[dupgenes_gem %in% genes]
    dupmsg <- paste0("Some input genes are duplicated in the matrix, including:\n",
                     paste(dupgenes_gem_IN_INPUT_GENES, collapse = ', '))

    # by default, we just drop the input duplicated genes
    if(drop.dup.genes == T){
      warning(dupmsg)
      genes <- genes[!(genes %in% dupgenes_gem)]
      res <- res[!(res$gene_symbol %in% dupgenes_gem),]
    } else{

      stop(dupmsg)

    }


    #if none of the input genes are dup genes, just drop them from the matrix
  } else{
    #using res, subset res to remove the duplicate genes
    res <- res[!(res$gene_symbol %in% dupgenes_gem),]
  }


  #subset the stats and gem
  res <- res[match(genes, res[,1]),]

  #gene stats: first three columns
  # gene name; LFC; pvalue or padj
  stats <- res[,1:3]

  #gem, everything except first 3 columns of res
  gem <- res[,-c(1:3)]

  #rownames of gem are first column of res
  rownames(gem) <- res[,1]


  #log or not

  if(do.log==T){
    gem <- log1p(gem)
  }

  #scale or not

  if(do.scale==T){
    gem <- t(scale(t(gem)))
  } else{
    gem <- as.matrix(gem)
  }



  # set up label, given logged / scaled gem
  # use user-provided if given, if not try to guess from scale/log options

  if(missing(legend_caption)){


    if(do.log==T & do.scale==T){

      legend_caption <- paste0('Scaled\nLogged\nNormalized\nCounts')

    }

    if(do.log==T & do.scale==F){

      legend_caption <- paste0('Logged\nNormalized\nCounts')

    }

    if(do.log==F & do.scale==T){

      legend_caption <- paste0('Scaled\nNormalized\nCounts')

    }

    if(do.log==F & do.scale==F){

      legend_caption <- paste0('Normalized\nCounts')

    }

  }

  #set up the column annotation
  # this is kind of complicated...



  #annotation df, has sample (column) name and color
  # metadata needs:
  # rownames of df need to be the column names of the samples
  # it needs a column called "Condition"
  # it needs a column called "Color"
  annotdf <- data.frame(Sample = colnames(gem),
                        # sampcheck = metadata[match(colnames(gem), rownames(metadata)), "Sample"],
                        Condition = metadata[match(colnames(gem), rownames(metadata)), "Condition"],
                        Color = metadata[match(colnames(gem), rownames(metadata)), "Color"])

  #define colors, need to use this list thing for complexheatmap
  # named vector: vector elements are colors, names are the values
  # not factor, each element needs it...?
  hacol <- list(Condition = annotdf$Color) ; names(hacol[[1]]) <- annotdf$Condition

  #create the annotation object
  ha <- ComplexHeatmap::HeatmapAnnotation(Condition = annotdf[,'Condition'], col = hacol, name = 'Condition')


  #set up the marker annotation
  # each gene is annotated according to celltype of origin
  #geneannots[,1] --> gene
  #geneannots[,2] --> celltype
  #geneannots[,3] --> color

  #define colors, need to use this list thing for complexheatmap
  hacol <- list(Celltype = geneannots[,3]) ; names(hacol[[1]]) <- geneannots[,2]

  #create the annotation object
  ha_gene <- ComplexHeatmap::rowAnnotation(Celltype = geneannots[,2], col = hacol)




  #sig marks
  # using stats df, apply thesholds of LFC (second column) and pvalue (third column)
  sigres <- stats[abs(stats[,2]) > lfc_thres,]
  sigres <- sigres[sigres[,3] < pval_thres,]

  #if the genes are in the significant res, we can add the significant marks
  rownames(gem) <- ifelse(rownames(gem) %in% sigres[,1], yes = paste0('* ', rownames(gem), ' *'), no = rownames(gem))


  ComplexHeatmap::Heatmap(gem,
                          top_annotation = ha,
                          right_annotation = ha_gene,
                          name = legend_caption,
                          ...
  )
}



