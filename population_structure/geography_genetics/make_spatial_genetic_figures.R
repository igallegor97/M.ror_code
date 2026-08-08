#!/usr/bin/env Rscript
# =============================================================
# make_spatial_genetic_figures.R
#
# Creates publication-oriented maps and plots:
#   1. Latin America sampling map
#   2. Colombia/Ecuador regional zoom
#   3. Isolation-by-distance scatter plots
#   4. Nearest-genetic-neighbor maps
#   5. Top-similarity network maps
#   6. Genetic/geographic heatmaps
#   7. Discrepant-pair plots
#   8. Genetic distance by geographic-distance class
# =============================================================

required_pkgs <- c(
  "ggplot2","ggrepel","sf","rnaturalearth",
  "rnaturalearthdata","ggspatial"
)
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg,quietly=TRUE)) {
    stop(sprintf("Required R package not installed: %s",pkg))
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(sf)
  library(rnaturalearth)
  library(ggspatial)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Usage: Rscript make_spatial_genetic_figures.R COORDS METADATA PAIRS OUTPUT")
}

COORD_FILE <- args[1]
METADATA_FILE <- args[2]
PAIR_FILE <- args[3]
OUTPUT_DIR <- args[4]
dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)

coords <- read.delim(COORD_FILE,stringsAsFactors=FALSE,check.names=FALSE)
meta <- read.delim(METADATA_FILE,stringsAsFactors=FALSE,check.names=FALSE)
pairs <- read.delim(PAIR_FILE,stringsAsFactors=FALSE,check.names=FALSE)

valid <- coords[
  !is.na(coords$latitude) & !is.na(coords$longitude) &
  coords$latitude!="" & coords$longitude!="",
  ,drop=FALSE
]
valid$latitude <- as.numeric(valid$latitude)
valid$longitude <- as.numeric(valid$longitude)
# Add only metadata columns that are not already present in the
# coordinate table. Both files contain country and region_state;
# merging all columns would create country.x/country.y.
metadata_to_add <- meta[, c(
  "sample_id",
  setdiff(colnames(meta), colnames(valid))
), drop = FALSE]

original_order <- valid$sample_id

valid <- merge(
  valid,
  metadata_to_add,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

# Restore the original coordinate-table order after merge.
valid <- valid[match(original_order, valid$sample_id), , drop = FALSE]

points_sf <- st_as_sf(
  valid,coords=c("longitude","latitude"),crs=4326,remove=FALSE
)
world <- ne_countries(scale="medium",returnclass="sf")
latin <- world[world$continent %in% c("South America","North America"),]

base_map <- ggplot() +
  geom_sf(data=latin,fill="grey96",color="grey65",linewidth=0.25) +
  geom_sf(
    data=points_sf,
    aes(color=country,shape=technology),
    size=3.2,alpha=0.95
  ) +
  geom_text_repel(
    data=valid,
    aes(x=longitude,y=latitude,label=sample_id,color=country),
    size=3,max.overlaps=Inf,show.legend=FALSE,
    box.padding=0.35
  ) +
  annotation_north_arrow(
    location="bl",which_north="true",
    pad_x=grid::unit(0.15,"in"),pad_y=grid::unit(0.45,"in"),
    style=north_arrow_fancy_orienteering
  ) +
  coord_sf(xlim=c(-105,-55),ylim=c(-25,25),expand=FALSE) +
  labs(
    title="Geographic distribution of Moniliophthora roreri isolates",
    subtitle="Color indicates country; shape indicates sequencing technology",
    color="Country",shape="Technology"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title=element_text(face="bold"),
    legend.position="right"
  )

ggsave(
  file.path(OUTPUT_DIR,"Fig1_sampling_map_Latin_America.pdf"),
  base_map,width=10,height=7,dpi=300,device="pdf"
)
ggsave(
  file.path(OUTPUT_DIR,"Fig1_sampling_map_Latin_America.png"),
  base_map,width=10,height=7,dpi=300
)

zoom <- ggplot() +
  geom_sf(data=latin,fill="grey96",color="grey65",linewidth=0.3) +
  geom_sf(
    data=points_sf,
    aes(color=country,shape=technology),
    size=3.5,alpha=0.95
  ) +
  geom_text_repel(
    data=valid,
    aes(x=longitude,y=latitude,label=sample_id,color=country),
    size=3,max.overlaps=Inf,show.legend=FALSE
  ) +
  coord_sf(xlim=c(-80.5,-72),ylim=c(-5.5,13),expand=FALSE) +
  labs(
    title="Regional view: Colombia and Ecuador",
    color="Country",shape="Technology"
  ) +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold"))

ggsave(
  file.path(OUTPUT_DIR,"Fig2_sampling_map_Colombia_Ecuador_zoom.pdf"),
  zoom,width=8,height=8,dpi=300,device="pdf"
)

# Isolation by distance plots
make_ibd <- function(response,profile_label,file_stub) {
  p <- ggplot(
    pairs,
    aes(
      x=geographic_distance_km,
      y=.data[[response]],
      color=same_country
    )
  ) +
    geom_point(alpha=0.65,size=2.2) +
    geom_smooth(
      aes(group=1),
      method="lm",formula=y~log10(x+1),
      se=TRUE,color="black",linewidth=0.8
    ) +
    scale_x_continuous(trans="log1p") +
    labs(
      title=paste("Isolation by distance —",profile_label),
      subtitle="Each point represents one pair of isolates",
      x="Geographic distance (km; log1p scale)",
      y="Genetic distance",
      color="Pair category"
    ) +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold"))
  ggsave(
    file.path(OUTPUT_DIR,paste0(file_stub,".pdf")),
    p,width=7.5,height=5.8,dpi=300,device="pdf"
  )
}
make_ibd("genetic_distance_all11","all 11 communities","Fig3A_IBD_all11")
make_ibd("genetic_distance_conservative9","conservative 9 communities","Fig3B_IBD_conservative9")

# Build geographic segments from pair rows
coord_lookup <- valid[,c("sample_id","longitude","latitude")]

add_coords <- function(df) {
  a <- coord_lookup[match(df$sample_1,coord_lookup$sample_id),]
  b <- coord_lookup[match(df$sample_2,coord_lookup$sample_id),]
  df$lon1 <- a$longitude
  df$lat1 <- a$latitude
  df$lon2 <- b$longitude
  df$lat2 <- b$latitude
  df
}

# Nearest genetic neighbor, one directed choice per sample, deduplicated.
nearest_edges <- function(response) {
  edges <- list()
  for (sid in valid$sample_id) {
    x <- pairs[pairs$sample_1==sid | pairs$sample_2==sid,,drop=FALSE]
    x <- x[order(x[[response]]),,drop=FALSE]
    if (nrow(x)>0) edges[[length(edges)+1]] <- x[1,,drop=FALSE]
  }
  e <- unique(do.call(rbind,edges))
  add_coords(e)
}

plot_network <- function(edges,response,title,file_name) {
  # geom_curve cannot draw a curve when both endpoints are identical.
  # This occurs for isolates represented by the same municipality-level
  # coordinates (for example Co44-CO8, Co52-Co82, or E7-Mr030).
  # Keep these pairs in the exported edge tables, but exclude them from
  # curved-line rendering and show only the overlapping sample points.
  curve_edges <- edges[
    !(is.na(edges$lon1) | is.na(edges$lat1) |
      is.na(edges$lon2) | is.na(edges$lat2)) &
    !(edges$lon1 == edges$lon2 & edges$lat1 == edges$lat2),
    ,drop=FALSE
  ]

  colocated_edges <- edges[
    !(is.na(edges$lon1) | is.na(edges$lat1) |
      is.na(edges$lon2) | is.na(edges$lat2)) &
    (edges$lon1 == edges$lon2 & edges$lat1 == edges$lat2),
    ,drop=FALSE
  ]

  p <- ggplot() +
    geom_sf(data=latin,fill="grey97",color="grey70",linewidth=0.25)

  if (nrow(curve_edges) > 0) {
    p <- p + geom_curve(
      data=curve_edges,
      aes(
        x=lon1,y=lat1,xend=lon2,yend=lat2,
        linewidth=1-.data[[response]],
        alpha=1-.data[[response]]
      ),
      curvature=0.12,color="grey20"
    )
  }

  # Mark municipality-level colocated genetic connections with a ring.
  if (nrow(colocated_edges) > 0) {
    colocated_points <- unique(colocated_edges[,c("lon1","lat1")])
    p <- p + geom_point(
      data=colocated_points,
      aes(x=lon1,y=lat1),
      shape=21,fill=NA,color="black",stroke=1.1,size=5,
      inherit.aes=FALSE
    )
  }

  p <- p +
    geom_sf(
      data=points_sf,
      aes(color=country,shape=technology),
      size=3.2
    ) +
    geom_text_repel(
      data=valid,
      aes(x=longitude,y=latitude,label=sample_id,color=country),
      size=2.8,max.overlaps=Inf,show.legend=FALSE
    ) +
    scale_linewidth_continuous(name="Genetic similarity") +
    scale_alpha_continuous(name="Genetic similarity",range=c(0.25,0.9)) +
    coord_sf(xlim=c(-105,-55),ylim=c(-25,25),expand=FALSE) +
    labs(title=title,color="Country",shape="Technology") +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold"))
  ggsave(file.path(OUTPUT_DIR,file_name),p,width=10,height=7,dpi=300,device="pdf")
}

near_all <- nearest_edges("genetic_distance_all11")
near_cons <- nearest_edges("genetic_distance_conservative9")
write.table(near_all,file.path(OUTPUT_DIR,"nearest_neighbor_edges_all11.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
write.table(near_cons,file.path(OUTPUT_DIR,"nearest_neighbor_edges_conservative9.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)

plot_network(
  near_all,"genetic_distance_all11",
  "Nearest genetic neighbor of each isolate - all 11 communities",
  "Fig4A_nearest_neighbor_map_all11.pdf"
)
plot_network(
  near_cons,"genetic_distance_conservative9",
  "Nearest genetic neighbor of each isolate - conservative 9 communities",
  "Fig4B_nearest_neighbor_map_conservative9.pdf"
)

# Top 10% most similar pairs
top_network <- function(response) {
  threshold <- quantile(pairs[[response]],0.10,na.rm=TRUE)
  add_coords(pairs[pairs[[response]]<=threshold,,drop=FALSE])
}
top_all <- top_network("genetic_distance_all11")
top_cons <- top_network("genetic_distance_conservative9")
write.table(top_all,file.path(OUTPUT_DIR,"top10pct_similarity_edges_all11.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
write.table(top_cons,file.path(OUTPUT_DIR,"top10pct_similarity_edges_conservative9.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)

plot_network(
  top_all,"genetic_distance_all11",
  "Top 10% most genetically similar pairs - all 11 communities",
  "Fig5A_top_similarity_network_all11.pdf"
)
plot_network(
  top_cons,"genetic_distance_conservative9",
  "Top 10% most genetically similar pairs - conservative 9 communities",
  "Fig5B_top_similarity_network_conservative9.pdf"
)

# Paired heatmaps
matrix_from_pairs <- function(value_col) {
  s <- valid$sample_id
  m <- matrix(0,nrow=length(s),ncol=length(s),dimnames=list(s,s))
  for (i in seq_len(nrow(pairs))) {
    a <- pairs$sample_1[i]; b <- pairs$sample_2[i]
    m[a,b] <- pairs[[value_col]][i]
    m[b,a] <- pairs[[value_col]][i]
  }
  m
}

plot_heatmap <- function(mat,title,file_name,digits=3) {
  long <- as.data.frame(as.table(mat))
  colnames(long) <- c("sample_1","sample_2","value")
  p <- ggplot(long,aes(x=sample_1,y=sample_2,fill=value)) +
    geom_tile() +
    geom_text(aes(label=sprintf(paste0("%.",digits,"f"),value)),size=2) +
    labs(title=title,x=NULL,y=NULL,fill="Distance") +
    theme_minimal(base_size=9) +
    theme(
      plot.title=element_text(face="bold"),
      axis.text.x=element_text(angle=60,hjust=1),
      panel.grid=element_blank()
    )
  ggsave(file.path(OUTPUT_DIR,file_name),p,width=9,height=8,dpi=300,device="pdf")
}

plot_heatmap(
  matrix_from_pairs("geographic_distance_km"),
  "Geographic distance among isolates (km)",
  "Fig6A_geographic_distance_heatmap.pdf",0
)
plot_heatmap(
  matrix_from_pairs("genetic_distance_all11"),
  "Genetic distance - all 11 communities",
  "Fig6B_genetic_distance_heatmap_all11.pdf",3
)
plot_heatmap(
  matrix_from_pairs("genetic_distance_conservative9"),
  "Genetic distance - conservative 9 communities",
  "Fig6C_genetic_distance_heatmap_conservative9.pdf",3
)

# Discrepant pairs: five low and five high for each profile
plot_discrepancy <- function(column,title,file_name) {
  ordered <- pairs[order(pairs[[column]]),]
  selected <- rbind(head(ordered,5),tail(ordered,5))
  selected$pair <- paste(selected$sample_1,selected$sample_2,sep=" - ")
  selected$interpretation <- ifelse(
    selected[[column]]<0,
    "genetically closer than geographic distance predicts",
    "genetically farther than geographic distance predicts"
  )
  selected$pair <- factor(selected$pair,levels=selected$pair[order(selected[[column]])])
  p <- ggplot(
    selected,
    aes(x=pair,y=.data[[column]],fill=interpretation)
  ) +
    geom_col(width=0.7) +
    coord_flip() +
    geom_hline(yintercept=0,linewidth=0.5) +
    labs(
      title=title,
      x=NULL,
      y="Standardized genetic–geographic discrepancy",
      fill="Interpretation"
    ) +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold"))
  ggsave(file.path(OUTPUT_DIR,file_name),p,width=9,height=6,dpi=300,device="pdf")
}

plot_discrepancy(
  "discrepancy_all11",
  "Pairs with strongest genetic–geographic discrepancy - all 11",
  "Fig7A_discrepant_pairs_all11.pdf"
)
plot_discrepancy(
  "discrepancy_conservative9",
  "Pairs with strongest genetic–geographic discrepancy - conservative 9",
  "Fig7B_discrepant_pairs_conservative9.pdf"
)

if ("distance_class" %in% colnames(pairs)) {
  long <- rbind(
    data.frame(distance_class=pairs$distance_class,
               profile="all11",genetic_distance=pairs$genetic_distance_all11),
    data.frame(distance_class=pairs$distance_class,
               profile="conservative9",genetic_distance=pairs$genetic_distance_conservative9)
  )
  pclass <- ggplot(
    long,aes(x=distance_class,y=genetic_distance,fill=profile)
  ) +
    geom_boxplot(outlier.alpha=0.4) +
    labs(
      title="Genetic distance across geographic-distance classes",
      x="Geographic-distance quantile class",
      y="Genetic distance",
      fill="Profile"
    ) +
    theme_classic(base_size=11) +
    theme(
      plot.title=element_text(face="bold"),
      axis.text.x=element_text(angle=35,hjust=1)
    )
  ggsave(
    file.path(OUTPUT_DIR,"Fig8_distance_class_boxplots.pdf"),
    pclass,width=9,height=5.5,dpi=300,device="pdf"
  )
}

cat("Spatial-genetic figures completed.\n")
