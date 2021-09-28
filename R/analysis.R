# https://r-pkgs.org/whole-game.html



### module score ###

#' Calculate a sample-wise module score given a list of genes
#'
#' this is inspired by Tirosh et al, Science (2016) and the AddModuleScore function in Seurat
#'
#' @param gem matrix, a gene expresssion matrix with rows=genes and columns=samples
#' @param genelist character vector, genes for which to calculate module score
#' @param numbins integer, how many bins to split the input genelist into, for selection of control genes
#' @param numcontrolgenesperbin integer, how many control genes to select per bin
#'
#' @return a vector of numerics, sample-wise module scores
#' @export
#'
#' @examples
modulescore <- function(gem, genelist, numbins=NULL, numcontrolgenesperbin=NULL){
  if(is.null(numbins)){numbins=24}
  if(is.null(numcontrolgenesperbin)){numcontrolgenesperbin=100}

  testgenes <- genelist
  controlgenes <- rownames(gem)
  controlgenes <- controlgenes[!(controlgenes %in% testgenes)]

  #testmeans <- Matrix::rowMeans(gem[rownames(gem) %in% genes,])

  ### bin the genes by avg expression ###

  #get average gene expression
  avgs <- Matrix::rowMeans(gem)
  avgs <- avgs[order(avgs)]

  #cut; get the gene bins
  bins <- cut_number(avgs, n=numbins)

  bindf <- data.frame(genenames = names(avgs),
                      bins = bins)

  #select control genes from same expression bins
  controlgenes <- c()
  for(genedex in 1:length(testgenes)){

    #get gene and bin
    gene <- testgenes[genedex]
    genebin <- bindf[bindf$genenames == gene,2]

    #select all genes from that bin, to sample from
    tmpcontrols <- bindf[bindf$bins == genebin,]

    #exclude the actual test genes
    tmpcontrols <- tmpcontrols[!(tmpcontrols$genenames %in% testgenes),]

    #if num controls exceeds bin size, take all genes in bin...
    numtotake <- ifelse(nrow(tmpcontrols) < numcontrolgenesperbin,
                        yes = nrow(tmpcontrols),
                        no = numcontrolgenesperbin)

    #sample the genes
    tmpcontrols <- sample(tmpcontrols$genenames, size = numtotake, replace = F)

    controlgenes <- unique(c(controlgenes,
                             tmpcontrols
    ))
  }



  #get control gene mean expression for each sample
  samplemeans_controls <- Matrix::colMeans( gem[rownames(gem) %in% controlgenes,] )

  #get test genes mean expression for each samoke
  samplemeans_test <- Matrix::colMeans( gem[rownames(gem) %in% testgenes,] )

  #subtract them to get module score
  modulescore <- samplemeans_test - samplemeans_controls

  return(modulescore)
}










### filtering outlier genes ###

# we have observed a strange phenomenon in which each individual replicate
# expresses certain outlier genes at a very high level.
# Thus we will try to define and filter out outlier genes.
# Outlier genes will be defined as genes with 5x higher expression
# than ALL other samples in the same genotype.
# if the gene has this characteristic, it is removed from both genotypes
#
#
#
#
# #normalize with DESeq2
# gemx <- t( t(gem) / DESeq2::estimateSizeFactorsForMatrix(gem) )
#
#
#
# badgenes <- c()
# for(genotype in unique(md$Condition) ){
#
#   #get the gem for each genotype
#   mdtmp <- md[md$Condition == genotype,]
#   tmpgem <- gemx[,colnames(gemx) %in% mdtmp$SampleID]
#
#
#   message('\nStarting ', genotype)
#
#
#   #progress bar
#   total = nrow(tmpgem) - 1
#   pb <- txtProgressBar(min = 0, max = total, style = 3)
#
#
#   genotype_badgenes <- c()
#   for(genedex in c(1:nrow(tmpgem)) ){
#
#     #get genename
#     genename <- rownames(tmpgem)[genedex]
#
#     #get gem for genotype
#     generow <- tmpgem[genedex,]
#
#     ## test if the max is higher than 5x any of the others... ##
#     #get the max gene
#     maxsamp <- max(generow)
#
#     #test if ALL of the others are hgiher than maxsamp / 5
#     maxdivfive <- maxsamp / 5
#
#     # actual test:
#     outliergene <- ifelse(test = all(test = generow[-which.max(generow)] < maxdivfive),
#                           yes = T, no = F)
#
#
#     if(outliergene){
#       genotype_badgenes <- c(genotype_badgenes, genename)
#     }
#
#
#     setTxtProgressBar(pb, genedex)
#
#   } #end gem sweep
#
#
#
#   badgenes <- c(badgenes, genotype_badgenes)
#
# }
#
# badgenes <- unique(badgenes)
#
# #remove these from the gem and test...
#
# gem <- gem[!(rownames(gem) %in% badgenes),]
#
# rm(gemx, badgenes, tmpgem, mdtmp, pb, genotype_badgenes)
# rm(genedex,genename,generow,maxdivfive,maxsamp, outliergene, total, genotype)







