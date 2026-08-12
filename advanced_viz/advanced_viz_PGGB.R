#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[[i + 1L]]
}
master_path <- arg("--master")
output_dir <- arg("--output-dir")
top_n <- as.integer(arg("--top-n", "15"))
png_dpi <- as.integer(arg("--png-dpi", "300"))
go_obo <- arg("--go-obo", "")
if (is.null(master_path) || is.null(output_dir)) stop("Required: --master TABLE --output-dir DIR")

plot_dir <- file.path(output_dir, "plots")
data_dir <- file.path(output_dir, "plot_data")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("MrorB3", "MrorC26", "MrorCO8", "MrorCO84", "MrorE7")
pal <- c(blue="#0072B2", sky="#56B4E9", green="#009E73", orange="#E69F00",
         vermillion="#D55E00", purple="#CC79A7", yellow="#F0E442", grey="#7A7A7A")
sample_pal <- setNames(unname(pal[c("blue","orange","green","purple","vermillion")]), samples)

theme_set(theme_minimal(base_size = 10, base_family = "sans") +
  theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
        plot.title=element_text(face="bold", size=rel(1.25)),
        plot.subtitle=element_text(colour="grey35"), axis.title=element_text(face="bold"),
        strip.text=element_text(face="bold"), legend.position="bottom",
        legend.title=element_text(face="bold"), plot.caption=element_text(colour="grey40", hjust=0)))

save_plot <- function(p, stem, width, height) {
  ggsave(file.path(plot_dir, paste0(stem, ".pdf")), p, width=width, height=height, units="in", device="pdf", bg="white")
  ggsave(file.path(plot_dir, paste0(stem, ".svg")), p, width=width, height=height, units="in", device=svglite::svglite, bg="white")
  ggsave(file.path(plot_dir, paste0(stem, ".png")), p, width=width, height=height, units="in", dpi=png_dpi, bg="white")
}
write_data <- function(x, name) fwrite(x, file.path(data_dir, name), sep="\t", na="")
bool <- function(x) { z <- toupper(trimws(as.character(x))); !is.na(z) & z %in% c("TRUE","T","1","YES") }

dt <- fread(master_path, sep="\t", quote="", na.strings=c("", "-"), showProgress=TRUE)
required <- c("sample_id","protein_id","protein_length","eggnog_description","GO_terms","KEGG_ko",
              "pfam_domains","cazy_families","signalp_positive","secretome_candidate","effector_candidate",
              "effectorp_primary_class","effectorp_dual_localized","effectorp_apoplastic_probability",
              "effectorp_cytoplasmic_probability")
missing <- setdiff(required, names(dt)); if (length(missing)) stop("Missing columns: ", paste(missing, collapse=", "))
if (anyDuplicated(dt[,.(sample_id,protein_id)])) stop("Duplicated protein keys")
dt[, sample_id := factor(sample_id, levels=samples)]
for (v in c("signalp_positive","secretome_candidate","effector_candidate","effectorp_dual_localized")) set(dt,j=v,value=bool(dt[[v]]))
dt[, protein_length := as.numeric(protein_length)]
if (any(dt$secretome_candidate & !dt$signalp_positive)) stop("Non-nested secretome detected")
if (any(dt$effector_candidate & !dt$secretome_candidate)) stop("Non-secretome effector detected")

split_terms <- function(column, pattern="[;,]") {
  x <- dt[!is.na(get(column)) & nzchar(get(column)), .(sample_id,protein_id,value=get(column))]
  x[, term := strsplit(value, pattern, perl=TRUE)]
  unique(x[,.(term=trimws(unlist(term))),by=.(sample_id,protein_id)])[nzchar(term) & term!="-"]
}

# 01 — UpSet of annotation-layer intersections.
up <- as.data.frame(dt[, .(
  eggNOG=!is.na(eggnog_description) & nzchar(eggnog_description),
  GO=!is.na(GO_terms) & nzchar(GO_terms), KEGG=!is.na(KEGG_ko) & nzchar(KEGG_ko),
  Pfam=!is.na(pfam_domains) & nzchar(pfam_domains), CAZy=!is.na(cazy_families) & nzchar(cazy_families),
  Secretome=secretome_candidate, Effector=effector_candidate
)])
p_upset <- ComplexUpset::upset(
  up, intersect=c("eggNOG","GO","KEGG","Pfam","CAZy","Secretome","Effector"),
  min_size=10, width_ratio=0.18, sort_intersections_by=c("cardinality","degree"),
  base_annotations=list("Intersection size"=ComplexUpset::intersection_size(text=list(size=3))),
  set_sizes=ComplexUpset::upset_set_size() + scale_y_continuous(labels=comma)
) + labs(title="Intersections among functional annotation layers",
         subtitle="Each protein contributes to exactly one displayed intersection pattern",
         caption="Small intersections (<10 proteins) are omitted from display but retained in the master table.")
save_plot(p_upset, "figure_01_annotation_upset", 13, 7.5)

# Prepare term profiles used by dot plots, clustering and PCA.
pfam <- split_terms("pfam_domains")
cazy <- split_terms("cazy_families", "[;,|]"); cazy[,term:=sub("\\(.*$","",term)]; cazy<-cazy[grepl("^(GH|GT|PL|CE|AA|CBM)[0-9]+",term)]
go <- split_terms("GO_terms", "[,;]")[grepl("^GO:[0-9]{7}$",term)]
ko <- split_terms("KEGG_ko", "[,;]"); ko[,term:=sub("^ko:","",term,ignore.case=TRUE)]; ko<-ko[grepl("^K[0-9]{5}$",term)]

go_map <- data.table(term=character(), name=character())
if (nzchar(go_obo)) {
  if (!file.exists(go_obo)) stop("GO OBO not found: ",go_obo)
  lines<-readLines(go_obo,warn=FALSE); ids<-grep("^id: GO:[0-9]{7}$",lines); ni<-grep("^name: ",lines)
  pos<-findInterval(ids,ni)+1L; nextid<-c(ids[-1L],length(lines)+1L); valid<-pos<=length(ni) & ni[pmin(pos,length(ni))]<nextid
  go_map<-data.table(term=sub("^id: ","",lines[ids[valid]]),name=sub("^name: ","",lines[ni[pos[valid]]]))
}

term_counts <- function(long, layer) {
  z <- unique(long[,.(sample_id,protein_id,term)])[, .N, by=.(sample_id,term)]
  totals <- z[,.(total=sum(N)),by=term][order(-total,term)]
  selected <- head(totals,top_n)$term
  z <- z[term %in% selected]
  z[, layer:=layer]
  z[, total_proteins:=dt[,.N,by=sample_id]$N[match(sample_id,dt[,.N,by=sample_id]$sample_id)]]
  z[, percent:=100*N/total_proteins]
  z
}
dot_data <- rbindlist(list(term_counts(pfam,"Pfam"),term_counts(cazy,"CAZy"),term_counts(go,"GO"),term_counts(ko,"KEGG KO")))
dot_data[, term_label := term]
if (nrow(go_map)) dot_data[layer=="GO", term_label:=ifelse(is.na(go_map$name[match(term,go_map$term)]),term,paste0(go_map$name[match(term,go_map$term)]," (",term,")"))]
dot_data[is.na(term_label),term_label:=term]
dot_data[,term_label:=factor(term_label,levels=rev(unique(term_label[order(layer,-N)])))]
write_data(dot_data,"figure_02_functional_dotplot_data.tsv")
dot_one <- function(layer_name) ggplot(dot_data[layer==layer_name],aes(sample_id,term_label,size=N,colour=percent))+
  geom_point(alpha=.88)+scale_size_area(max_size=9,labels=comma)+
  scale_colour_viridis_c(option="C",labels=label_percent(scale=1))+
  labs(title=layer_name,x=NULL,y=NULL,size="Proteins",colour="Proteome coverage")+
  theme(panel.grid.major=element_line(colour="grey90"))
p_dot <- (dot_one("Pfam")|dot_one("CAZy"))/(dot_one("GO")|dot_one("KEGG KO"))+
  plot_annotation(title="Functional profiles across five proteomes",subtitle="Point area shows protein count; colour shows percentage of the complete proteome",tag_levels="A")
save_plot(p_dot,"figure_02_functional_dotplots",16,13)

# 03 — Clustered heatmap of combined functional profiles.
heat <- copy(dot_data)
heat[, feature:=paste(layer,term_label,sep=" · ")]
heat[, zscore:=if (.N>1 && sd(percent)>0) as.numeric(scale(percent)) else 0,by=feature]
wide <- dcast(heat,feature~sample_id,value.var="zscore",fill=0)
row_order <- if(nrow(wide)>1) hclust(dist(as.matrix(wide[,-1])))$order else 1L
sample_matrix <- t(as.matrix(wide[,-1])); colnames(sample_matrix)<-wide$feature
col_order <- if(nrow(sample_matrix)>1) hclust(dist(sample_matrix))$order else seq_len(nrow(sample_matrix))
heat[,feature:=factor(feature,levels=rev(wide$feature[row_order]))]
heat[,sample_id:=factor(sample_id,levels=rownames(sample_matrix)[col_order])]
write_data(heat,"figure_03_clustered_heatmap_data.tsv")
p_heat <- ggplot(heat,aes(sample_id,feature,fill=zscore))+geom_tile(colour="white",linewidth=.3)+
  scale_fill_gradient2(low=pal["blue"],mid="white",high=pal["vermillion"],midpoint=0,oob=squish)+
  labs(title="Clustered functional-profile heatmap",subtitle="Row-standardized proteome coverage; clustering uses Euclidean distance",x=NULL,y=NULL,fill="Row z-score")+
  theme(axis.text.y=element_text(size=6),panel.grid=element_blank())
save_plot(p_heat,"figure_03_clustered_functional_heatmap",9,14)

# 04 — PCA of Hellinger-transformed profiles using all observed terms.
all_profiles <- rbindlist(list(
  unique(pfam)[,feature:=paste0("Pfam:",term)], unique(cazy)[,feature:=paste0("CAZy:",term)],
  unique(go)[,feature:=paste0("GO:",term)], unique(ko)[,feature:=paste0("KO:",term)]
),use.names=TRUE,fill=TRUE)[,.(N=.N),by=.(sample_id,feature)]
profile_wide <- dcast(all_profiles,sample_id~feature,value.var="N",fill=0)
mat <- as.matrix(profile_wide[,-1]); rownames(mat)<-profile_wide$sample_id
mat <- sqrt(mat/rowSums(mat)); keep <- apply(mat,2,sd)>0; mat<-mat[,keep,drop=FALSE]
pca <- prcomp(mat,center=TRUE,scale.=FALSE)
pca_dt <- data.table(sample_id=factor(rownames(pca$x),levels=samples),PC1=pca$x[,1],PC2=pca$x[,2])
variance <- 100*pca$sdev^2/sum(pca$sdev^2)
write_data(pca_dt,"figure_04_functional_pca_scores.tsv")
p_pca <- ggplot(pca_dt,aes(PC1,PC2,colour=sample_id,label=sample_id))+geom_hline(yintercept=0,colour="grey85")+
  geom_vline(xintercept=0,colour="grey85")+geom_point(size=4)+ggrepel::geom_text_repel(show.legend=FALSE,size=3.5)+
  scale_colour_manual(values=sample_pal,drop=FALSE)+
  labs(title="PCA of combined functional profiles",subtitle="Pfam, CAZy, GO and KEGG KO counts; Hellinger transformation",
       x=paste0("PC1 (",number(variance[1],accuracy=.1),"%)"),y=paste0("PC2 (",number(variance[2],accuracy=.1),"%)"),colour="Proteome",
       caption="Descriptive ordination of five proteomes; it is not an inferential test.")+
  coord_equal()
save_plot(p_pca,"figure_04_functional_profile_PCA",8,6.5)

# 05 — Alluvial SignalP → DeepTMHMM → EffectorP flow.
flow <- dt[signalp_positive==TRUE, .N, by=.(
  sample_id,
  DeepTMHMM=ifelse(secretome_candidate,"Final secretome","Excluded by DeepTMHMM"),
  EffectorP=ifelse(!secretome_candidate,"Not tested",ifelse(!effector_candidate,"Non-effector",as.character(effectorp_primary_class)))
)]
write_data(flow,"figure_05_secretome_alluvial_data.tsv")
p_alluvial <- ggplot(flow,aes(axis1=sample_id,axis2=DeepTMHMM,axis3=EffectorP,y=N))+
  ggalluvial::geom_alluvium(aes(fill=EffectorP),width=1/12,alpha=.78)+
  ggalluvial::geom_stratum(width=1/9,fill="white",colour="grey35")+
  geom_text(stat=ggalluvial::StatStratum,aes(label=after_stat(stratum)),size=2.8)+
  scale_x_discrete(limits=c("Proteome","DeepTMHMM","EffectorP"),expand=c(.08,.08))+
  scale_fill_manual(values=c("Not tested"=unname(pal["grey"]),"Non-effector"=unname(pal["sky"]),"Apoplastic"=unname(pal["orange"]),"Cytoplasmic"=unname(pal["purple"])))+
  scale_y_continuous(labels=comma)+labs(title="Secretome-to-effector alluvial flow",subtitle="Flow begins with SignalP-positive proteins",x=NULL,y="Proteins",fill="Final outcome")+
  theme(panel.grid=element_blank())
save_plot(p_alluvial,"figure_05_secretome_effector_alluvial",13,7)

# 06 — Raincloud distributions for biologically relevant subsets.
length_groups <- rbindlist(list(
  dt[,.(sample_id,protein_length,group="Complete proteome")],
  dt[(signalp_positive),.(sample_id,protein_length,group="SignalP positive")],
  dt[(secretome_candidate),.(sample_id,protein_length,group="Final secretome")],
  dt[(effector_candidate),.(sample_id,protein_length,group="Effector candidates")]
))
length_groups[,group:=factor(group,levels=c("Complete proteome","SignalP positive","Final secretome","Effector candidates"))]
write_data(length_groups[,.(sample_id,protein_length,group)],"figure_06_raincloud_length_data.tsv")
p_rain <- ggplot(length_groups,aes(group,protein_length,fill=group,colour=group))+
  ggdist::stat_halfeye(adjust=.7,width=.7,.width=0,justification=-.20,point_colour=NA,alpha=.65)+
  geom_boxplot(width=.13,outlier.shape=NA,alpha=.8,colour="grey20")+
  coord_cartesian(ylim=c(0,quantile(length_groups$protein_length,.99)))+
  facet_wrap(~sample_id,nrow=1)+scale_fill_manual(values=setNames(unname(pal[c("blue","sky","green","purple")]),levels(length_groups$group)),drop=FALSE)+
  scale_colour_manual(values=setNames(unname(pal[c("blue","sky","green","purple")]),levels(length_groups$group)),drop=FALSE)+
  labs(title="Protein-length raincloud distributions",subtitle="Violin density and boxplot; axis truncated at the combined 99th percentile",x=NULL,y="Protein length (aa)")+
  theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="none")
save_plot(p_rain,"figure_06_protein_length_raincloud",15,6)

# 07 — EffectorP probability plane with censored secondary probabilities.
prob <- dt[effector_candidate==TRUE,.(sample_id,protein_id,effectorp_primary_class,effectorp_dual_localized,
  apo=as.numeric(effectorp_apoplastic_probability),cyto=as.numeric(effectorp_cytoplasmic_probability))]
prob[,apo_censored:=is.na(apo)]; prob[,cyto_censored:=is.na(cyto)]
prob[is.na(apo),apo:=.5]; prob[is.na(cyto),cyto:=.5]
prob[,reporting:=ifelse(apo_censored|cyto_censored,"One probability censored at ≤0.5","Both probabilities reported")]
write_data(prob,"figure_07_effectorp_probability_data.tsv")
p_prob <- ggplot(prob,aes(apo,cyto,colour=effectorp_primary_class,shape=reporting))+
  geom_hline(yintercept=.5,linetype=2,colour="grey60")+geom_vline(xintercept=.5,linetype=2,colour="grey60")+
  geom_point(alpha=.65,size=1.8)+facet_wrap(~sample_id)+
  scale_colour_manual(values=c("Apoplastic"=unname(pal["orange"]),"Cytoplasmic"=unname(pal["purple"])))+
  scale_shape_manual(values=c("Both probabilities reported"=16,"One probability censored at ≤0.5"=17))+
  coord_equal(xlim=c(.49,1),ylim=c(.49,1))+labs(title="EffectorP 3 probability landscape",
    subtitle="Unreported secondary probabilities are displayed at the 0.5 censoring boundary, not interpreted as exact values",
    x="Apoplastic probability",y="Cytoplasmic probability",colour="Primary class",shape="Portal reporting")
save_plot(p_prob,"figure_07_effectorp_probability_scatter",12,7)

# 08 — Hierarchical CAZy sunburst (class → family), combined proteomes.
sun <- unique(cazy[,.(protein_id,sample_id,term)])[,.(count=.N),by=term]
sun[,class:=sub("^([A-Z]+).*","\\1",term)]
sun<-sun[class %in% c("AA","CBM","CE","GH","GT","PL")]
sun[,rank:=frank(-count,ties.method="first"),by=class]
sun[rank>10,term:="Other"]; sun<-sun[,.(count=sum(count)),by=.(class,term)]
setorder(sun,class,-count,term); sun[,xmax:=cumsum(count)]; sun[,xmin:=shift(xmax,fill=0)]
classes<-sun[,.(xmin=min(xmin),xmax=max(xmax),count=sum(count)),by=class]
sun[,ymin:=1]; sun[,ymax:=2]; classes[,`:=`(ymin=0,ymax=1)]
write_data(sun,"figure_08_cazy_sunburst_families.tsv")
p_sun <- ggplot()+geom_rect(data=classes,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=class),colour="white")+
  geom_rect(data=sun,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=class),colour="white",linewidth=.2,alpha=.78)+
  geom_text(data=classes,aes(x=(xmin+xmax)/2,y=.5,label=class),size=3,fontface="bold")+
  geom_text(data=sun[term!="Other" & count>=quantile(count,.55)],aes(x=(xmin+xmax)/2,y=1.52,label=term),size=2.2,angle=0,check_overlap=TRUE)+
  coord_polar(theta="x")+xlim(0,sum(classes$count))+ylim(-.25,2.05)+
  scale_fill_manual(values=c("AA"=unname(pal["orange"]),"CBM"=unname(pal["yellow"]),"CE"=unname(pal["purple"]),"GH"=unname(pal["blue"]),"GT"=unname(pal["green"]),"PL"=unname(pal["vermillion"])))+
  labs(title="Consensus CAZyme hierarchy",subtitle="Inner ring: CAZy class; outer ring: families (top 10 per class plus Other)",fill="CAZy class")+
  theme_void()+theme(plot.title=element_text(face="bold"),plot.subtitle=element_text(colour="grey35"),legend.position="bottom")
save_plot(p_sun,"figure_08_cazy_sunburst",9,9)

manifest <- data.table(figure=sprintf("%02d",1:8),description=c(
  "Annotation-layer UpSet plot","Functional dot plots","Clustered functional heatmap","Functional-profile PCA",
  "Secretome-effector alluvial flow","Protein-length raincloud","EffectorP probability plane","CAZy sunburst"))
write_data(manifest,"advanced_viz_manifest.tsv")
write_data(data.table(input=normalizePath(master_path),md5=unname(tools::md5sum(master_path)),generated=format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"),top_n=top_n,png_dpi=png_dpi),"run_metadata.tsv")
capture.output(sessionInfo(),file=file.path(output_dir,"R_sessionInfo.txt"))
cat("ADVANCED VIZ COMPLETED\nOutput:",output_dir,"\n")
