library(tidyverse)
library(ggh4x)
library(ggpubr)
library(ggpmisc)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggsci)
library(CellPlot)
library(circlize)
library(RColorBrewer)
library(fgsea)
library(cowplot)
library(magick)
library(GeneOverlap)

results <- read_rds("./data/epigenetic_age_acceleration.rds") %>% filter(baseline==1&timepoint==2)
protein <- read_rds("./data/proteomics_difference.rds")
metabo <- read_rds("./data/metabolomics_difference.rds")
lipid <- read_rds("./data/lipidomics_difference.rds")
glycan <- read_rds("./data/glycomics_difference.rds")
cell <- read_rds("./data/cytomics_difference.rds")

combined <- rbind(results %>% left_join(reshape2::melt(protein) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Proteomics"),
                  results %>% left_join(reshape2::melt(metabo) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Metabolomics"),
                  results %>% left_join(reshape2::melt(lipid) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Lipidomics"),
                  results %>% left_join(reshape2::melt(glycan) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Glycomics"),
                  results %>% left_join(reshape2::melt(cell) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Cytomics"))

combined_r <- combined %>% group_by(omics,group,clock,feature) %>% summarise(n = nrow(tibble(age_acc_difference,value)%>% na.omit), r = cor(age_acc_difference,value)) %>% filter(n>3)

compare_sham <- combined_r %>% filter(group!="Sham") %>% dplyr::select(omics,group,clock,feature,r) %>% 
  left_join(combined_r %>% filter(group=="Sham") %>% ungroup %>% dplyr::select(omics,clock,feature,r) %>% set_names("omics","clock","feature","sham_r"))
compare_sham <- compare_sham %>% group_by(omics,group,feature) %>% na.omit %>% summarise(mean_r = mean(r), mean_sham = mean(sham_r), p = t.test(psych::fisherz(r),mu = 0)$p.value, p_sham = t.test(psych::fisherz(r),psych::fisherz(sham_r))$p.value)
compare_sham <- compare_sham %>% group_by(omics) %>% summarise(mean_r = mean_r, group = group, feature = feature, fdr = p.adjust(p), fdr_sham = p.adjust(p_sham))
writexl::write_xlsx(compare_sham %>% set_names("omics","r","group","feature","fdr_zero","fdr_sham"), path = "./output/tables/Sup_Table2.xlsx")
keep <- compare_sham %>% group_by(omics,group) %>% summarise(sig = feature[which(fdr<0.01&fdr_sham<0.01)]) %>% pull(sig) 

#figure 2a
mat <- reshape2::acast(combined_r %>% filter(feature%in%keep) %>% dplyr::group_by(omics,group,feature) %>% summarise(mean_r = mean(r)), feature~group, value.var = "mean_r")
mat[is.na(mat)] <- 0
col_fun1 <- colorRamp2(seq(-1,1,0.2),rev(brewer.pal(11, "RdBu")))

mat_sig <- combined_r %>% filter(feature%in%keep) %>% dplyr::group_by(omics,group,feature) %>% summarise(mean_r = mean(r))
mat_sig <- mat_sig %>% left_join(compare_sham %>% group_by(omics,group) %>% summarise(feature = feature[which(fdr<0.01&fdr_sham<0.01)]) %>% mutate(significant = 1))
mat_sig <- reshape2::acast(mat_sig, feature~group, value.var = "significant")
mat_sig[is.na(mat_sig)] <- 0

row_dend <- hclust(dist(mat))  
row_order <- row_dend$order    
mat <- mat[row_order, ]
mat_sig <- mat_sig[row_order, ]

mat_sig_proteomics <- mat_sig[which(sapply(rownames(mat), function(x) unique(compare_sham$omics[compare_sham$feature==x]) %>% na.omit)%in%"Proteomics"),]
mat_sig_metabolomics <- mat_sig[which(sapply(rownames(mat), function(x) unique(compare_sham$omics[compare_sham$feature==x]) %>% na.omit)%in%"Metabolomics"),]
mat_sig_cytomics <- mat_sig[which(sapply(rownames(mat), function(x) unique(compare_sham$omics[compare_sham$feature==x]) %>% na.omit)%in%"Cytomics"),]
mat_sig_glycomics <- mat_sig[which(sapply(rownames(mat), function(x) unique(compare_sham$omics[compare_sham$feature==x]) %>% na.omit)%in%"Glycomics"),]
mat_sig_lipidomics <- mat_sig[which(sapply(rownames(mat), function(x) unique(compare_sham$omics[compare_sham$feature==x]) %>% na.omit)%in%"Lipidomics"),]

#links
net_matrix <- combined %>% filter(group=="TPE + IVIG\n(B)") %>% group_by(patient,feature) %>% summarise(dif = mean(age_acc_difference), value = mean(value))
net_matrix <- reshape2::dcast(net_matrix, patient~feature, value.var = "value")
net_matrix$patient <- NULL
df <- net_matrix[,rownames(mat)]

df_link <- expand_grid(from_index = 1:ncol(df), to_index = 1:ncol(df))
df_link <- df_link %>% left_join(tibble(from_index = 1:ncol(df), from_name = colnames(df)))
df_link <- df_link %>% left_join(tibble(to_index = 1:ncol(df), to_name = colnames(df)))
df_link <- df_link %>% left_join(compare_sham[,c(4,1)] %>% set_names("from_name","from_omics") %>% unique)
df_link <- df_link %>% left_join(compare_sham[,c(4,1)] %>% set_names("to_name","to_omics") %>% unique)
df_link <- df_link %>% filter(from_index!=to_index)
df_link <- df_link %>% filter(from_omics!=to_omics)
df_link <- df_link %>% filter(from_omics=="Cytomics")

custom_corr <- function(x, y) {
  points <- tibble(x = df[,x], y = df[,y]) %>% na.omit()
  if (nrow(points) >= 3) {
    return(cor(points$x,points$y))
  } else {
    return(NA)
  }
}

df_link <- df_link %>% mutate(r = unlist(map2(from_index,to_index,custom_corr)))
df_link <- df_link %>% filter(abs(r)>0.90)
df_link <- df_link %>% filter(from_name%in%c(mat_sig[,1] %>% enframe %>% filter(value==1) %>% pull(name)))

pdf(file = "./output/figures/2a.pdf", height = 7, width = 7)
circos.par(start.degree = 80, gap.degree = c(5,5,5,5,20))
circos.heatmap(mat, 
               rownames.side = "outside",
               dend.side = "none",
               track.height = 0.2,
               cell.border = "white",
               cell.lty = 1,
               cluster = FALSE,
               cell.lwd = 0.5,
               rownames.cex = 0.5,
               split = sapply(rownames(mat), function(x) unique(compare_sham$omics[compare_sham$feature==x]) %>% na.omit), 
               show.sector.labels = TRUE,
               col = col_fun1)

for (r in seq_len(nrow(mat_sig_proteomics))) {
  for (c in seq_len(ncol(mat_sig_proteomics))) {
    if (mat_sig_proteomics[r, c] == 1) {
      circos.rect(xleft = r-1, xright = r,
                  ybottom = 4-c, ytop = 4-c+1,
                  sector.index = "Proteomics", 
                  track.index = 2,
                  border = "black", lwd = 1)
    }
  }
}

for (r in seq_len(nrow(mat_sig_metabolomics))) {
  for (c in seq_len(ncol(mat_sig_metabolomics))) {
    if (mat_sig_metabolomics[r, c] == 1) {
      circos.rect(xleft = r-1, xright = r,
                  ybottom = 4-c, ytop = 4-c+1,
                  sector.index = "Metabolomics", 
                  track.index = 2,
                  border = "black", lwd = 1)
    }
  }
}

for (r in seq_len(nrow(mat_sig_cytomics))) {
  for (c in seq_len(ncol(mat_sig_cytomics))) {
    if (mat_sig_cytomics[r, c] == 1) {
      circos.rect(xleft = r-1, xright = r,
                  ybottom = 4-c, ytop = 4-c+1,
                  sector.index = "Cytomics", 
                  track.index = 2,
                  border = "black", lwd = 1)
    }
  }
}

for (r in seq_len(nrow(mat_sig_glycomics))) {
  for (c in seq_len(ncol(mat_sig_glycomics))) {
    if (mat_sig_glycomics[r, c] == 1) {
      circos.rect(xleft = r-1, xright = r,
                  ybottom = 4-c, ytop = 4-c+1,
                  sector.index = "Glycomics", 
                  track.index = 2,
                  border = "black", lwd = 1)
    }
  }
}

for (r in seq_len(nrow(mat_sig_lipidomics))) {
  for (c in seq_len(ncol(mat_sig_lipidomics))) {
    if (mat_sig_lipidomics[r, c] == 1) {
      circos.rect(xleft = r-1, xright = r,
                  ybottom = 4-c, ytop = 4-c+1,
                  sector.index = "Lipidomics", 
                  track.index = 2,
                  border = "black", lwd = 1)
    }
  }
}

for(i in seq_len(nrow(df_link))) {
  circos.heatmap.link(df_link$from_index[i],
                      df_link$to_index[i],
                      h.ratio = 0.7,
                      col = ifelse(df_link$r[i]>0.9, "#F4A582","#92C5DE"),
                      rou = 0.5)
}

circos.clear()
dev.off()

#figure 2d
# pval_permutation <- function(x, y){
#   if(nrow(tibble(x,y)%>%na.omit)<3){
#     return(NA)
#   } else {
#     observed_correlation <- cor(x, y, method = "pearson")
#     permute_test <- function(x, y, num_permutations = 10000) {
#       n <- length(x)
#       set.seed(1234)
#       permuted_correlations <- replicate(num_permutations, {
#         permuted_y <- sample(y)  
#         cor(x, permuted_y, method = "pearson")  
#       })
#       p_value <- mean(abs(permuted_correlations) >= abs(observed_correlation))
#       return(p_value)
#     }
#     p_value <- permute_test(x, y)
#     return(p_value)
#   }
# }

top <- compare_sham %>% filter(group=="TPE + IVIG\n(B)"&fdr_sham<0.05) %>% mutate(sign = sign(mean_r)) %>% group_by(omics,sign) %>% summarise(feature = feature[which.min(fdr)])
correlations_top <- combined %>% filter(group=="TPE + IVIG\n(B)"&feature%in%top$feature) %>% group_by(omics,patient,feature) %>% summarise(dif = mean(age_acc_difference), value = mean(value))
correlations_top <- correlations_top %>% left_join(top[,c("feature","sign")])
correlations_top$direction <- ifelse(correlations_top$sign==1,"Positive","Negative")
correlations_top <- correlations_top %>% left_join(correlations_top %>% group_by(direction,omics) %>% summarise(cor = cor(dif,value)))

p1 <- ggplot(correlations_top, aes(dif,value))+
  geom_point(shape = 21, fill = "#009988")+
  facet_grid2(direction ~ omics, scales = "free", independent = "all")+
  geom_text(aes(label = feature, x = -Inf, y = Inf), hjust = -0.1, vjust = 1.5, color = "#009988", fontface="bold", stat = "unique")+
  stat_correlation(data = correlations_top %>% filter(direction=="Negative"), label.x = 0.05, label.y = 0.05, method = "pearson", aes(label = paste0("R == ", round(cor,2))))+
  stat_correlation(data = correlations_top %>% filter(direction=="Positive"), label.x = 0.95, label.y = 0.05, method = "pearson", aes(label = paste0("R == ", round(cor,2))))+
  geom_smooth(method = "lm", color = "#009988", alpha = 0.2)+
  geom_ribbon(stat = "smooth",method = "lm",se = TRUE, alpha = 0, colour = "black", linetype = "dotted")+
  theme_pubr(border = TRUE)+
  labs(x = "Age acceleration difference", y = "Feature level")
p1
write_rds(p1, file = "./output/figures/2d.rds")

#figure 2b
proteomics <- compare_sham %>% filter(omics=="Proteomics"&group=="TPE + IVIG\n(B)"&fdr<0.05&fdr_sham<0.05)
enrich_result <- enrichGO(
    gene = proteomics$feature,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    maxGSSize = 500,
    minGSSize = 50,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1)@result
enrich_result <- enrich_result[nchar(enrich_result$Description)<50,]
enrich_result$geneID <- lapply(enrich_result$geneID, function(x) strsplit(x, split = "\\/")[[1]])
enrich_result$r <- lapply(enrich_result$geneID, function(x) proteomics$mean_r[proteomics$feature%in%x])
enrich_result <- enrich_result[1:30,]
npg_colors <- pal_npg("nrc")(10)

pdf(file = "./output/figures/2b.pdf",width=7.31,height=7.42)
cell.plot(x = setNames(-log10(enrich_result$p.adjust), enrich_result$Description),
          cells = enrich_result$r,
          main ="",
          xlab = "-log10(FDR)",
          x.mar = c(0.7, 0),
          cell.col	= c(npg_colors[4],"white",npg_colors[8]),
          key.n = 5,
          xlab.ticks	= 2,
          y.mar = c(.1, 0),
          cex = 1.6,
          key = TRUE,
          sym = TRUE,
          key.lab	= "Correlation",
          cell.limit = 1,
          cell.outer = 1,
          cell.bounds = c(-1,1),
          bar.scale = .1,
          space = 0.1)
dev.off()

#figure 2c
geneset <- read_rds("../hoa/data/geneset_chatgpt_4o_mini.rds")
hyper <- tibble(pathway = names(geneset), pval = sapply(1:12, function(x) testGeneOverlap(newGeneOverlap(listA = proteomics$feature, listB = geneset[[x]]))@pval))
hyper$score <- -log10(hyper$pval)
hyper$pathway <- str_to_sentence(hyper$pathway)
hyper$pathway <- gsub(" ","\n",hyper$pathway)
hyper <- hyper %>% mutate(pathway = factor(pathway, levels = unique(pathway)))
hyper$fdr <- p.adjust(hyper$pval)

p2 <- ggplot(hyper, aes(x = pathway, y = score, fill = score)) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  coord_polar(clip = "off") +  # Ensure labels are not clipped
  geom_hline(yintercept = max(hyper$score) * 1.05, color = "black", linetype = "solid") +  # Add a circular line
  geom_segment(aes(x = as.numeric(pathway), xend = as.numeric(pathway), y = 0, yend = max(score) * 1.1), linetype = "dotted", color = "black") +  # Add dotted lines
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, face = "bold"),  # Keep text size
    axis.text.y = element_blank(),  # Optional: remove y-axis labels
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",  # Keep the legend on the right side
    legend.justification = "center",  # Center the legend
    legend.title = element_text(hjust = 0.5, size = 12),  
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Increase plot margins
  ) +
  scale_y_continuous(limits = c(0, max(hyper$score) * 1.3)) +  # Adjust y-axis limits to reduce plot size
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd")) +  # Use YlOrRd palette
  labs(x = "", y = "", fill = "-log10(p-value)")
p2
write_rds(p2, file = "./output/figures/2c.rds")

#combine panels
a <- ggdraw() + draw_image(image_read_pdf("./output/figures/2a.pdf", pages = 1) %>% image_trim()) 
b <- ggdraw() + draw_image(image_read_pdf("./output/figures/2b.pdf", pages = 1) %>% image_trim()) 
c <- read_rds("./output/figures/2c.rds")
d <- read_rds("./output/figures/2d.rds")

pdf(file = "./output/figures/figure2.pdf", width=11.54, height=10.71)
plot_grid(plot_grid(plot_grid(a,labels = c("a")),
          plot_grid(b,c, ncol = 1, labels = c("b","c"), rel_heights = c(1,0.7)),
          rel_widths = c(1,0.7)),plot_grid(d,labels = c("d")),nrow = 2, rel_heights = c(1,0.5))
dev.off()
