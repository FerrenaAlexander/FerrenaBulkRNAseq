% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEplots.R
\name{markerheatmap}
\alias{markerheatmap}
\title{Marker annotated heatmap}
\usage{
markerheatmap(
  res,
  geneannots,
  metadata,
  lfc_thres,
  pval_thres,
  do.scale,
  do.log,
  legend_caption,
  drop.dup.genes,
  ...
)
}
\arguments{
\item{res}{a data.frame. one row for each gene. first three columns must be 1) gene symbol or id; 2) fold change; and 3) pvalue or padj. rest of the columns are normalized counts.}

\item{geneannots}{a data.frame. one row for each gene. first column must be gene, second column must be celltype, third column must be color for each celltype}

\item{metadata}{a data.frame. sample metadata, each row s a sample. must have rownames = colnames of samples in res (ie, colnames of \code{res[,-c(1:3)]} ); must also have a column called "Condition" for each sample condition, and a column called "Color" for each sample color; if two conditions, must have two colors}

\item{lfc_thres}{numeric, a cutoff for lfc to be considered significant, default is 1}

\item{pval_thres}{numeric, a cutoff for pval/padj to be considered significant, default is 0.05}

\item{do.scale}{T/F, whether to scale genes, default is T}

\item{do.log}{T/F, whether to log1p genes, default is T}

\item{legend_caption}{string, caption to use above figure legend, by default will guess based on scale/log options}

\item{drop.dup.genes}{T/F whether to drop duplicated genes, by default T, if set to F and dup genes are detected will throw an error}

\item{...}{additional options passed on to \code{ComplexHeatmap::Heatmap()}, see \code{?ComplexHeatmap::Heatmap()}}
}
\description{
Heatmaps with annotations for samples (columns) and rows (genes), useful to plot celltype markers
}
