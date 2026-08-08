#!/usr/bin/env Rscript
# =============================================================
# run_balanced_ibd_robustness.R
#
# Tests whether the genetic-geographic relationship persists
# after sampling equal numbers of SNPs per PGGB community.
#
# For each replicate:
#   - sample N SNPs per community
#   - calculate normalized Manhattan genetic distances
#   - compare with geographic distances
#   - record Pearson, Spearman, and descriptive slope
#
# Usage:
# Rscript run_balanced_ibd_robustness.R \
#   MATRIX FEATURE_METADATA GEO_MATRIX PROFILE OUTPUT_DIR \
#   REPLICATES SEED N_PER_COMMUNITY
#
# N_PER_COMMUNITY=0 means use the smallest available community.
# =============================================================

required_pkgs <- c("ggplot2")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg,quietly=TRUE)) {
    stop(sprintf("Required R package not installed: %s",pkg))
  }
}
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=8) {
  stop("Usage: Rscript run_balanced_ibd_robustness.R MATRIX FEATURES GEO PROFILE OUTPUT REPS SEED N_PER_COMMUNITY")
}

MATRIX_FILE <- args[1]
FEATURE_FILE <- args[2]
GEO_FILE <- args[3]
PROFILE <- args[4]
OUTPUT_DIR <- args[5]
REPLICATES <- as.integer(args[6])
BASE_SEED <- as.integer(args[7])
N_REQUESTED <- as.integer(args[8])

dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)

matrix_df <- read.delim(
  MATRIX_FILE,row.names=1,check.names=FALSE,stringsAsFactors=FALSE
)
features <- read.delim(
  FEATURE_FILE,row.names=1,check.names=FALSE,stringsAsFactors=FALSE
)
geo <- as.matrix(read.delim(
  GEO_FILE,row.names=1,check.names=FALSE,stringsAsFactors=FALSE
))

common_features <- intersect(rownames(matrix_df),rownames(features))
matrix_df <- matrix_df[common_features,,drop=FALSE]
features <- features[common_features,,drop=FALSE]

samples <- intersect(colnames(matrix_df),rownames(geo))
matrix_df <- matrix_df[,samples,drop=FALSE]
geo <- geo[samples,samples,drop=FALSE]

counts <- table(features$community)
n_per <- if (N_REQUESTED>0) min(N_REQUESTED,min(counts)) else min(counts)

sampling_table <- data.frame(
  community=names(counts),
  available_snps=as.integer(counts),
  selected_per_replicate=n_per,
  selected_pct=round(n_per/as.integer(counts)*100,4)
)
write.table(
  sampling_table,
  file.path(OUTPUT_DIR,paste0(PROFILE,"_balanced_sampling_plan.tsv")),
  sep="\t",quote=FALSE,row.names=FALSE
)

geo_vec <- as.vector(geo[upper.tri(geo)])
log_geo <- log10(geo_vec+1)

results <- vector("list",REPLICATES)

for (rep in seq_len(REPLICATES)) {
  set.seed(BASE_SEED+rep-1)

  selected <- unlist(lapply(
    split(rownames(features),features$community),
    function(ids) sample(ids,n_per,replace=FALSE)
  ),use.names=FALSE)

  x <- t(as.matrix(matrix_df[selected,,drop=FALSE]))
  variable <- apply(x,2,var)>0
  x <- x[,variable,drop=FALSE]

  d <- as.matrix(dist(x,method="manhattan"))/ncol(x)
  dvec <- as.vector(d[upper.tri(d)])

  fit <- lm(dvec ~ log_geo)

  results[[rep]] <- data.frame(
    profile=PROFILE,
    replicate=rep,
    seed=BASE_SEED+rep-1,
    communities=length(counts),
    snps_per_community=n_per,
    snps_selected=length(selected),
    variable_snps=ncol(x),
    pearson=cor(dvec,geo_vec,method="pearson"),
    spearman=cor(dvec,geo_vec,method="spearman"),
    slope_log10_km=coef(fit)[2],
    r_squared=summary(fit)$r.squared
  )
}

result_df <- do.call(rbind,results)
write.table(
  result_df,
  file.path(OUTPUT_DIR,paste0(PROFILE,"_balanced_ibd_replicates.tsv")),
  sep="\t",quote=FALSE,row.names=FALSE
)

summarize_metric <- function(x) {
  c(
    median=median(x),
    q025=quantile(x,0.025),
    q975=quantile(x,0.975),
    mean=mean(x),
    sd=sd(x),
    positive_fraction=mean(x>0)
  )
}

metrics <- c("pearson","spearman","slope_log10_km","r_squared")
summary_rows <- do.call(rbind,lapply(metrics,function(metric) {
  values <- summarize_metric(result_df[[metric]])
  data.frame(
    profile=PROFILE,
    metric=metric,
    statistic=names(values),
    value=as.numeric(values),
    row.names=NULL
  )
}))
write.table(
  summary_rows,
  file.path(OUTPUT_DIR,paste0(PROFILE,"_balanced_ibd_summary.tsv")),
  sep="\t",quote=FALSE,row.names=FALSE
)

long <- rbind(
  data.frame(replicate=result_df$replicate,metric="Pearson",value=result_df$pearson),
  data.frame(replicate=result_df$replicate,metric="Spearman",value=result_df$spearman),
  data.frame(replicate=result_df$replicate,metric="Slope",value=result_df$slope_log10_km)
)

p <- ggplot(long,aes(x=metric,y=value,fill=metric)) +
  geom_violin(trim=FALSE,alpha=0.65) +
  geom_boxplot(width=0.15,outlier.shape=NA) +
  geom_hline(yintercept=0,linetype="dashed",linewidth=0.5) +
  labs(
    title=paste("Balanced isolation-by-distance robustness —",PROFILE),
    subtitle=sprintf("%d replicates; %d SNPs per community",REPLICATES,n_per),
    x=NULL,y="Statistic"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title=element_text(face="bold"),
    legend.position="none"
  )

ggsave(
  file.path(OUTPUT_DIR,paste0(PROFILE,"_balanced_ibd_robustness.pdf")),
  p,width=7,height=5.5,dpi=300,device="pdf"
)

cat(sprintf(
  "Balanced IBD robustness completed: %s, %d replicates, %d SNPs/community\n",
  PROFILE,REPLICATES,n_per
))
