#!/usr/bin/env Rscript

# =============================================================
# analyze_24G_population_structure.R
#
# Runs PCA, distance analysis, and UPGMA clustering for:
#   1. all 24 genomes
#   2. Illumina-only genomes
#   3. PacBio-only genomes
#
# Color = country
# Shape = sequencing technology
# =============================================================

required_packages <- c("ggplot2", "ggrepel")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required R package not installed: %s", pkg))
  }
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript analyze_24G_population_structure.R MATRIX METADATA LABEL OUTPUT_DIR")
}

MATRIX_FILE <- args[1]
METADATA_FILE <- args[2]
LABEL <- args[3]
OUTPUT_DIR <- args[4]

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

matrix_df <- read.delim(
  MATRIX_FILE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata <- read.delim(
  METADATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required <- c("sample_id","country","region_state","year","genome_size_mb","technology")
missing <- setdiff(required, colnames(metadata))
if (length(missing) > 0) {
  stop(sprintf("Metadata missing columns: %s", paste(missing, collapse=", ")))
}
if (!all(metadata$sample_id %in% colnames(matrix_df))) {
  stop("At least one metadata sample is absent from the SNP matrix.")
}
matrix_df <- matrix_df[, metadata$sample_id, drop=FALSE]

run_one <- function(sample_ids, analysis_name) {
  analysis_dir <- file.path(OUTPUT_DIR, analysis_name)
  dir.create(analysis_dir, recursive=TRUE, showWarnings=FALSE)

  x <- t(as.matrix(matrix_df[, sample_ids, drop=FALSE]))
  keep <- apply(x, 2, var) > 0
  x <- x[, keep, drop=FALSE]

  if (nrow(x) < 3 || ncol(x) < 2) {
    warning(sprintf("Skipping %s: insufficient samples or variable SNPs.", analysis_name))
    return(NULL)
  }

  pca <- prcomp(x, center=TRUE, scale.=FALSE)
  variance <- summary(pca)$importance[2,] * 100
  n_axes <- min(5, ncol(pca$x))
  scores <- as.data.frame(pca$x[,seq_len(n_axes),drop=FALSE])
  scores$sample_id <- rownames(scores)
  scores <- merge(scores, metadata, by="sample_id", all.x=TRUE, sort=FALSE)
  scores <- scores[match(sample_ids, scores$sample_id),,drop=FALSE]

  var_df <- data.frame(
    PC=paste0("PC",seq_along(variance)),
    pct_variance=variance,
    cumulative_variance=cumsum(variance)
  )

  write.table(
    scores,
    file.path(analysis_dir,paste0(LABEL,"_",analysis_name,"_pca_coordinates.tsv")),
    sep="\t",quote=FALSE,row.names=FALSE
  )
  write.table(
    var_df,
    file.path(analysis_dir,paste0(LABEL,"_",analysis_name,"_variance_explained.tsv")),
    sep="\t",quote=FALSE,row.names=FALSE
  )

  pc1 <- sprintf("PC1 (%.1f%%)",variance[1])
  pc2 <- sprintf("PC2 (%.1f%%)",variance[2])

  p12 <- ggplot(
    scores,
    aes(x=PC1,y=PC2,color=country,shape=technology,label=sample_id)
  ) +
    geom_point(size=3.8,alpha=0.9) +
    geom_text_repel(size=3,box.padding=0.35,max.overlaps=Inf,show.legend=FALSE) +
    labs(
      title=sprintf("%s — %s",LABEL,analysis_name),
      subtitle=sprintf("%d samples; %s variable SNPs",
                       nrow(x),format(ncol(x),big.mark=",")),
      x=pc1,y=pc2,color="Country",shape="Technology"
    ) +
    theme_classic(base_size=12) +
    theme(
      plot.title=element_text(face="bold"),
      panel.grid.major=element_line(color="grey92",linewidth=0.4),
      legend.position="right"
    )

  ggsave(
    file.path(analysis_dir,paste0(LABEL,"_",analysis_name,"_PC1vsPC2.pdf")),
    p12,width=8.2,height=6.3,dpi=300,device="pdf"
  )

  if ("PC3" %in% colnames(scores)) {
    pc3 <- sprintf("PC3 (%.1f%%)",variance[3])
    p13 <- ggplot(
      scores,
      aes(x=PC1,y=PC3,color=country,shape=technology,label=sample_id)
    ) +
      geom_point(size=3.8,alpha=0.9) +
      geom_text_repel(size=3,box.padding=0.35,max.overlaps=Inf,show.legend=FALSE) +
      labs(
        title=sprintf("%s — %s",LABEL,analysis_name),
        x=pc1,y=pc3,color="Country",shape="Technology"
      ) +
      theme_classic(base_size=12) +
      theme(
        plot.title=element_text(face="bold"),
        panel.grid.major=element_line(color="grey92",linewidth=0.4)
      )
    ggsave(
      file.path(analysis_dir,paste0(LABEL,"_",analysis_name,"_PC1vsPC3.pdf")),
      p13,width=8.2,height=6.3,dpi=300,device="pdf"
    )
  }

  n_scree <- min(15,length(variance))
  scree <- data.frame(
    PC=factor(paste0("PC",seq_len(n_scree)),levels=paste0("PC",seq_len(n_scree))),
    pct_variance=variance[seq_len(n_scree)]
  )
  ps <- ggplot(scree,aes(x=PC,y=pct_variance)) +
    geom_col(width=0.65) +
    geom_line(aes(group=1),linewidth=0.7) +
    geom_point(size=2.2) +
    labs(
      title=sprintf("Scree plot — %s — %s",LABEL,analysis_name),
      x="Principal component",y="Variance explained (%)"
    ) +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold"),
          axis.text.x=element_text(angle=45,hjust=1))
  ggsave(
    file.path(analysis_dir,paste0(LABEL,"_",analysis_name,"_screeplot.pdf")),
    ps,width=7,height=4.5,dpi=300,device="pdf"
  )

  dmat <- as.matrix(dist(x,method="manhattan")) / ncol(x)
  write.table(
    dmat,
    file.path(analysis_dir,paste0(LABEL,"_",analysis_name,"_distance_matrix.tsv")),
    sep="\t",quote=FALSE,col.names=NA
  )

  dlong <- as.data.frame(as.table(dmat))
  colnames(dlong) <- c("sample_1","sample_2","distance")
  ph <- ggplot(dlong,aes(x=sample_1,y=sample_2,fill=distance)) +
    geom_tile() +
    geom_text(aes(label=sprintf("%.3f",distance)),size=2.3) +
    labs(
      title=sprintf("Genetic-distance heatmap — %s — %s",LABEL,analysis_name),
      x=NULL,y=NULL,fill="Distance"
    ) +
    theme_minimal(base_size=10) +
    theme(
      plot.title=element_text(face="bold"),
      axis.text.x=element_text(angle=60,hjust=1),
      panel.grid=element_blank()
    )
  ggsave(
    file.path(analysis_dir,paste0(LABEL,"_",analysis_name,"_distance_heatmap.pdf")),
    ph,width=9,height=8,dpi=300,device="pdf"
  )

  hc <- hclust(as.dist(dmat),method="average")
  pdf(
    file.path(analysis_dir,paste0(LABEL,"_",analysis_name,"_UPGMA_dendrogram.pdf")),
    width=8,height=6
  )
  plot(hc,main=sprintf("UPGMA — %s — %s",LABEL,analysis_name),
       xlab="",sub="",cex=0.8)
  dev.off()

  invisible(TRUE)
}

run_one(metadata$sample_id,"all_24_genomes")
run_one(metadata$sample_id[metadata$technology=="Illumina"],"illumina_only")
run_one(metadata$sample_id[metadata$technology=="PacBio"],"pacbio_only")

cat("Population-structure analyses completed.\n")
