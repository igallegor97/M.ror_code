#!/usr/bin/env Rscript
# =============================================================
# build_and_analyze_spatial_genetics.R
#
# Integrates PGGB genetic distances with curated coordinates.
# Produces geographic matrices, pairwise tables, Mantel tests,
# MRM models, descriptive regressions, distance classes, and
# discrepant-pair rankings.
# =============================================================

required_pkgs <- c("geosphere","vegan","ecodist")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg,quietly=TRUE)) {
    stop(sprintf("Required R package not installed: %s",pkg))
  }
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop("Usage: Rscript build_and_analyze_spatial_genetics.R COORDS METADATA ALL11 CONS9 OUTPUT")
}

COORD_FILE <- args[1]
METADATA_FILE <- args[2]
ALL11_FILE <- args[3]
CONS9_FILE <- args[4]
OUTPUT_DIR <- args[5]
dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)

coords <- read.delim(COORD_FILE,stringsAsFactors=FALSE,check.names=FALSE)
meta <- read.delim(METADATA_FILE,stringsAsFactors=FALSE,check.names=FALSE)

read_distance <- function(path) {
  x <- read.delim(path,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
  as.matrix(x)
}

all11 <- read_distance(ALL11_FILE)
cons9 <- read_distance(CONS9_FILE)

valid <- coords[
  !is.na(coords$latitude) & !is.na(coords$longitude) &
  coords$latitude!="" & coords$longitude!="",
  ,drop=FALSE
]
valid$latitude <- as.numeric(valid$latitude)
valid$longitude <- as.numeric(valid$longitude)

samples <- valid$sample_id
if (!all(samples %in% rownames(all11)) || !all(samples %in% rownames(cons9))) {
  stop("At least one coordinate sample is absent from a genetic distance matrix.")
}

all11 <- all11[samples,samples,drop=FALSE]
cons9 <- cons9[samples,samples,drop=FALSE]

xy <- as.matrix(valid[,c("longitude","latitude")])
geo_m <- geosphere::distm(xy,fun=geosphere::distGeo)
geo_km <- geo_m/1000
rownames(geo_km) <- samples
colnames(geo_km) <- samples

write.table(
  geo_km,file.path(OUTPUT_DIR,"geographic_distance_matrix_km.tsv"),
  sep="\t",quote=FALSE,col.names=NA
)

upper_pairs <- function(mat) {
  idx <- which(upper.tri(mat),arr.ind=TRUE)
  data.frame(
    sample_1=rownames(mat)[idx[,1]],
    sample_2=colnames(mat)[idx[,2]],
    value=mat[idx],
    stringsAsFactors=FALSE
  )
}

geo_long <- upper_pairs(geo_km)
all_long <- upper_pairs(all11)
cons_long <- upper_pairs(cons9)

pairs <- geo_long
names(pairs)[3] <- "geographic_distance_km"
pairs$genetic_distance_all11 <- all_long$value
pairs$genetic_distance_conservative9 <- cons_long$value

m1 <- meta[match(pairs$sample_1,meta$sample_id),]
m2 <- meta[match(pairs$sample_2,meta$sample_id),]

pairs$country_1 <- m1$country
pairs$country_2 <- m2$country
pairs$region_1 <- m1$region_state
pairs$region_2 <- m2$region_state
pairs$technology_1 <- m1$technology
pairs$technology_2 <- m2$technology
pairs$same_country <- ifelse(pairs$country_1==pairs$country_2,"same country","different countries")
pairs$same_region <- ifelse(
  pairs$country_1==pairs$country_2 & pairs$region_1==pairs$region_2,
  "same region","different regions"
)
pairs$technology_pair <- ifelse(
  pairs$technology_1==pairs$technology_2,
  paste0(pairs$technology_1,"-",pairs$technology_2),
  "mixed"
)

pairs$geo_z <- as.numeric(scale(log10(pairs$geographic_distance_km+1)))
pairs$all11_z <- as.numeric(scale(pairs$genetic_distance_all11))
pairs$cons9_z <- as.numeric(scale(pairs$genetic_distance_conservative9))
pairs$discrepancy_all11 <- pairs$all11_z-pairs$geo_z
pairs$discrepancy_conservative9 <- pairs$cons9_z-pairs$geo_z

# Quantile-based distance classes
breaks <- unique(quantile(
  pairs$geographic_distance_km,
  probs=seq(0,1,length.out=6),na.rm=TRUE
))
if (length(breaks)>=3) {
  pairs$distance_class <- cut(
    pairs$geographic_distance_km,breaks=breaks,
    include.lowest=TRUE,dig.lab=6
  )
}

write.table(
  pairs,file.path(OUTPUT_DIR,"geographic_genetic_pair_table.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE
)

# Mantel tests for the complete coordinate set and sensitivity subsets.
subset_definitions <- list(
  all_valid_coordinates = samples,
  high_precision = valid$sample_id[
    valid$coordinate_precision %in% c("locality","municipality")
  ],
  Colombia_only = valid$sample_id[valid$country=="Colombia"],
  Illumina_only = valid$sample_id[
    meta$technology[match(valid$sample_id,meta$sample_id)]=="Illumina"
  ],
  excluding_Colombia = valid$sample_id[valid$country!="Colombia"]
)

mantel_rows <- list()
subset_summary <- list()

for (subset_name in names(subset_definitions)) {
  subset_samples <- unique(subset_definitions[[subset_name]])
  subset_samples <- subset_samples[subset_samples %in% samples]

  subset_summary[[length(subset_summary)+1]] <- data.frame(
    subset=subset_name,
    n_samples=length(subset_samples),
    samples=paste(subset_samples,collapse=",")
  )

  if (length(subset_samples)<4) {
    next
  }

  geo_subset <- geo_km[subset_samples,subset_samples,drop=FALSE]

  for (profile in c("all11","conservative9")) {
    g_full <- if (profile=="all11") all11 else cons9
    g_subset <- g_full[subset_samples,subset_samples,drop=FALSE]

    for (method in c("pearson","spearman")) {
      x <- vegan::mantel(
        as.dist(g_subset),as.dist(geo_subset),
        method=method,permutations=9999,na.rm=TRUE
      )
      mantel_rows[[length(mantel_rows)+1]] <- data.frame(
        profile=profile,
        subset=subset_name,
        method=method,
        n_samples=length(subset_samples),
        statistic=unname(x$statistic),
        p_value=x$signif,
        permutations=x$permutations
      )
    }
  }
}

write.table(
  do.call(rbind,subset_summary),
  file.path(OUTPUT_DIR,"spatial_subset_summary.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE
)

mantel_df <- do.call(rbind,mantel_rows)
write.table(
  mantel_df,file.path(OUTPUT_DIR,"mantel_results.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE
)

# MRM: genetic distance ~ geographic distance + technology difference
technology <- meta$technology[match(samples,meta$sample_id)]
tech_mat <- outer(technology,technology,FUN=function(a,b) as.numeric(a!=b))
diag(tech_mat) <- 0
rownames(tech_mat) <- samples
colnames(tech_mat) <- samples

run_mrm <- function(gmat,profile) {
  fit <- ecodist::MRM(
    as.dist(gmat) ~ as.dist(log10(geo_km+1)) + as.dist(tech_mat),
    nperm=9999
  )
  coef <- as.data.frame(fit$coef)
  coef$term <- rownames(coef)
  rownames(coef) <- NULL
  coef$profile <- profile
  coef
}
mrm <- rbind(
  run_mrm(all11,"all11"),
  run_mrm(cons9,"conservative9")
)
write.table(
  mrm,file.path(OUTPUT_DIR,"mrm_results.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE
)

# Descriptive ordinary regressions (effect-size summaries only)
lm_rows <- list()
for (profile in c("all11","conservative9")) {
  response <- if (profile=="all11") "genetic_distance_all11" else "genetic_distance_conservative9"
  fit <- lm(
    as.formula(paste(response,"~ log10(geographic_distance_km+1)")),
    data=pairs
  )
  s <- summary(fit)
  lm_rows[[length(lm_rows)+1]] <- data.frame(
    profile=profile,
    intercept=coef(fit)[1],
    slope_log10_km=coef(fit)[2],
    r_squared=s$r.squared,
    adjusted_r_squared=s$adj.r.squared
  )
}
write.table(
  do.call(rbind,lm_rows),
  file.path(OUTPUT_DIR,"descriptive_regression_results.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE
)

if ("distance_class" %in% colnames(pairs)) {
  class_summary <- aggregate(
    cbind(genetic_distance_all11,genetic_distance_conservative9) ~ distance_class,
    data=pairs,
    FUN=function(x) c(n=length(x),median=median(x),mean=mean(x),sd=sd(x))
  )
  write.table(
    class_summary,file.path(OUTPUT_DIR,"distance_class_summary.tsv"),
    sep="\t",quote=FALSE,row.names=FALSE
  )
}

write.table(
  pairs[order(pairs$discrepancy_all11),],
  file.path(OUTPUT_DIR,"discrepant_pairs_all11.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE
)
write.table(
  pairs[order(pairs$discrepancy_conservative9),],
  file.path(OUTPUT_DIR,"discrepant_pairs_conservative9.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE
)

cat("Spatial-genetic statistics completed.\n")
