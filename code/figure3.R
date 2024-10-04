library(tidyverse)
library(impute)
library(glmnet)

results <- read_rds("./data/epigenetic_age_acceleration.rds") %>% filter(baseline==1&timepoint==2)
biomarker <- read_rds("./data/clinical_baseline.rds")
biomarker$feature[biomarker$feature=="EGFR"] <- "eGFR"
protein <- read_rds("./data/proteomics_baseline.rds")
metabo <- read_rds("./data/metabolomics_baseline.rds")
metabo$feature[metabo$feature=="Creatinine"] <- "Creatinine (M)"
lipid <- read_rds("./data/lipidomics_baseline.rds")
glycan <- read_rds("./data/glycomics_baseline.rds")
cell <- read_rds("./data/cytomics_baseline.rds")

combined <- rbind(results %>% left_join(reshape2::melt(biomarker) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Clinical"),
                  results %>% left_join(reshape2::melt(protein) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Proteomics"),
                  results %>% left_join(reshape2::melt(metabo) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Metabolomics"),
                  results %>% left_join(reshape2::melt(lipid) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Lipidomics"),
                  results %>% left_join(reshape2::melt(glycan) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Glycomics"),
                  results %>% left_join(reshape2::melt(cell) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Cytomics"))

combined_r <- combined %>% group_by(omics,group,clock,feature) %>% summarise(n = nrow(tibble(age_acc_difference,value)%>% na.omit), r = cor(age_acc_difference,value)) %>% filter(n>3)

compare_sham <- combined_r %>% filter(group!="Sham") %>% dplyr::select(omics,group,clock,feature,r) %>% 
  left_join(combined_r %>% filter(group=="Sham") %>% ungroup %>% dplyr::select(omics,clock,feature,r) %>% set_names("omics","clock","feature","sham_r"))
compare_sham <- compare_sham %>% group_by(omics,group,feature) %>% na.omit %>% summarise(mean_r = mean(r), mean_sham = mean(sham_r), p = t.test(psych::fisherz(r),mu = 0)$p.value, p_sham = t.test(psych::fisherz(r),psych::fisherz(sham_r))$p.value)
compare_sham <- compare_sham %>% group_by(omics) %>% summarise(mean_r = mean_r, group = group, feature = feature, fdr = p.adjust(p), fdr_sham = p.adjust(p_sham))
writexl::write_xlsx(compare_sham %>% set_names("omics","r","group","feature","fdr_zero","fdr_sham"), path = "./output/tables/Sup_Table3.xlsx")
keep <- compare_sham %>% group_by(omics,group) %>% summarise(sig = feature[which(fdr<0.01&fdr_sham<0.01)]) %>% pull(sig) 

#figure 3a
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

mat_sig_clinical <- mat_sig[which(sapply(rownames(mat), function(x) unique(compare_sham$omics[compare_sham$feature==x]) %>% na.omit)%in%"Clinical"),]
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
df_link <- df_link %>% filter(from_omics=="Clinical")

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

pdf(file = "./output/figures/3a.pdf", height = 7, width = 7)
circos.par(start.degree = 80, gap.degree = c(5,5,5,5,5,20))
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

for (r in seq_len(nrow(mat_sig_clinical))) {
  for (c in seq_len(ncol(mat_sig_clinical))) {
    if (mat_sig_clinical[r, c] == 1) {
      circos.rect(xleft = r-1, xright = r,
                  ybottom = 4-c, ytop = 4-c+1,
                  sector.index = "Clinical", 
                  track.index = 2,
                  border = "black", lwd = 1)
    }
  }
}

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

#figure 3b
top <- compare_sham %>% filter(group!="Sham"&omics=="Clinical"&fdr_sham<0.05) %>% mutate(sign = sign(mean_r)) %>% group_by(group,sign) %>% summarise(feature = feature[which.min(fdr)])
correlations_top <- top %>% left_join(combined) %>% filter(group!="Sham"&omics=="Clinical") %>% group_by(group,patient,feature,sign) %>% summarise(dif = mean(age_acc_difference), value = mean(value))
correlations_top$direction <- ifelse(correlations_top$sign==1,"Positive","Negative")
correlations_top <- correlations_top %>% left_join(correlations_top %>% group_by(direction,group) %>% summarise(cor = cor(dif,value)))

p1 <- ggplot(correlations_top %>% na.omit, aes(dif,value))+
  geom_point(shape = 21, fill = "#009988")+
  facet_grid2(direction ~ group, scales = "free", independent = "all")+
  geom_text(aes(label = feature), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, color = "#009988", fontface="bold", stat = "unique")+
  stat_correlation(data = correlations_top %>% filter(direction=="Negative"), label.x = 0.05, label.y = 0.05, method = "pearson", aes(label = paste0("R == ", round(cor,2))))+
  stat_correlation(data = correlations_top %>% filter(direction=="Positive"), label.x = 0.95, label.y = 0.05, method = "pearson", aes(label = paste0("R == ", round(cor,2))))+
  geom_smooth(method = "lm", color = "#009988", alpha = 0.2)+
  geom_ribbon(stat = "smooth",method = "lm",se = TRUE, alpha = 0, colour = "black", linetype = "dotted")+
  theme_pubr(border = TRUE)+
  labs(x = "Age acceleration difference", y = "Feature level")
p1
write_rds(p1, file = "./output/figures/3b.rds")

clinical <- biomarker
clinical <- clinical[,colSums(is.na(clinical))<134]
clinical <- clinical[rowSums(is.na(clinical[,2:ncol(clinical)]))<(38*0.5),]
clinical <- clinical[,colSums(is.na(clinical))<(nrow(clinical)*0.5)]
clinical <- clinical %>% column_to_rownames(var = "feature")
imputed <- impute.knn(as.matrix(clinical))$data %>% t()
lambdas <- 10^seq(-5, 5, by = 0.1)
preds <- tibble()
coefs <- tibble()

#predictors
age_acc <- results %>% filter(group=="TPE + IVIG\n(B)") %>% filter(patient%in%rownames(imputed))
imputed_expanded <- imputed[match(age_acc$patient,rownames(imputed)),]
set.seed(1234)
age_model <- cv.glmnet(x = imputed_expanded,
                       y = ifelse(age_acc$age_acc_difference<0,1,0),
                       type.measure="auc",
                       family="binomial",
                       alpha = 0.5,
                       keep=TRUE,
                       lambda = lambdas,
                       nfolds = 10)
preds <- rbind(preds, tibble(real = ifelse(age_acc$age_acc_difference<0,1,0), predicted = age_model$fit.preval[,age_model$lambda==age_model$lambda.min] %>% enframe %>% pull(2)) %>% mutate(group = "TPE + IVIG\n(B)"))
coefs <- rbind(coefs, coef(age_model, s = age_model$lambda.min) %>% {tibble(gene = rownames(.), coefficient = as.matrix(.)[,1])} %>% filter(coefficient != 0) %>% mutate(group = "TPE + IVIG\n(B)"))

age_acc <- results %>% filter(group=="TPE\n(B)") %>% filter(patient%in%rownames(imputed))
imputed_expanded <- imputed[match(age_acc$patient,rownames(imputed)),]
set.seed(1234)
age_model <- cv.glmnet(x = imputed_expanded,
                       y = ifelse(age_acc$age_acc_difference<0,1,0),
                       type.measure="auc",
                       family="binomial",
                       alpha = 0.5,
                       keep=TRUE,
                       lambda = lambdas,
                       nfolds = 10)
preds <- rbind(preds, tibble(real = ifelse(age_acc$age_acc_difference<0,1,0), predicted = age_model$fit.preval[,age_model$lambda==age_model$lambda.min] %>% enframe %>% pull(2)) %>% mutate(group = "TPE\n(B)"))
coefs <- rbind(coefs, coef(age_model, s = age_model$lambda.min) %>% {tibble(gene = rownames(.), coefficient = as.matrix(.)[,1])} %>% filter(coefficient != 0) %>% mutate(group = "TPE\n(B)"))

age_acc <- results %>% filter(group=="TPE\n(M)") %>% filter(patient%in%rownames(imputed))
imputed_expanded <- imputed[match(age_acc$patient,rownames(imputed)),]
set.seed(1234)
age_model <- cv.glmnet(x = imputed_expanded,
                       y = ifelse(age_acc$age_acc_difference<0,1,0),
                       type.measure="auc",
                       family="binomial",
                       alpha = 0.5,
                       keep=TRUE,
                       lambda = lambdas,
                       nfolds = 10)
preds <- rbind(preds, tibble(real = ifelse(age_acc$age_acc_difference<0,1,0), predicted = age_model$fit.preval[,age_model$lambda==age_model$lambda.min] %>% enframe %>% pull(2)) %>% mutate(group = "TPE\n(M)"))
coefs <- rbind(coefs, coef(age_model, s = age_model$lambda.min) %>% {tibble(gene = rownames(.), coefficient = as.matrix(.)[,1])} %>% filter(coefficient != 0) %>% mutate(group = "TPE\n(M)"))
preds$predicted_outcome <- ifelse(preds$predicted > 0.5, 1, 0)

all_curves <- tibble()
tpe_ivig <- preds %>% filter(group=="TPE + IVIG\n(B)")
auc <- auc_roc(tpe_ivig$predicted, tpe_ivig$real)
curve <- rbind(tibble(FPR = 0, TPR = 0), auc_roc(tpe_ivig$predicted, tpe_ivig$real, returnDT=TRUE)[,c("CumulativeFPR","CumulativeTPR")] %>% set_names("FPR","TPR")) %>% mutate(group = "TPE + IVIG (B)", auc)
all_curves <- rbind(all_curves,curve)

tpe_b <- preds %>% filter(group=="TPE\n(B)")
auc <- auc_roc(tpe_b$predicted, tpe_b$real)
curve <- rbind(tibble(FPR = 0, TPR = 0), auc_roc(tpe_b$predicted, tpe_b$real, returnDT=TRUE)[,c("CumulativeFPR","CumulativeTPR")] %>% set_names("FPR","TPR")) %>% mutate(group = "TPE (B)", auc)
all_curves <- rbind(all_curves,curve)

tpe_m <- preds %>% filter(group=="TPE\n(M)")
auc <- auc_roc(tpe_m$predicted, tpe_m$real)
curve <- rbind(tibble(FPR = 0, TPR = 0), auc_roc(tpe_m$predicted, tpe_m$real, returnDT=TRUE)[,c("CumulativeFPR","CumulativeTPR")] %>% set_names("FPR","TPR")) %>% mutate(group = "TPE (M)", auc)
all_curves <- rbind(all_curves,curve)

all_curves$group <- factor(all_curves$group, levels = unique(all_curves$group))
p2 <- ggplot(all_curves, aes(FPR, TPR, color = group, group = group))+
  geom_point()+
  geom_line()+
  geom_abline(slope = 1, lty = 2)+
  stat_correlation(data = all_curves %>% filter(group=="TPE + IVIG (B)"), label.x = 0.8, label.y = 0.3, method = "pearson", aes(label = paste0("AUC == ", round(auc,2))))+
  stat_correlation(data = all_curves %>% filter(group=="TPE (B)"), label.x = 0.8, label.y = 0.2, method = "pearson", aes(label = paste0("AUC == ", round(auc,2))))+
  stat_correlation(data = all_curves %>% filter(group=="TPE (M)"), label.x = 0.8, label.y = 0.1, method = "pearson", aes(label = paste0("AUC == ", round(auc,2))))+
  theme_pubr()+
  scale_color_manual(values = c("#ddaa33","#bb5566","#004488"))+
  labs(x = "False Positive Rate", y = "True Positive Rate", color = "")
write_rds(p2, file = "./output/figures/3c.rds")
writexl::write_xlsx(coefs %>% set_names("biomarker","coefficient","group"), path = "./output/tables/Sup_Table4.xlsx")

#coefficients network
library(GGally)
library(network)

pairs <- coefs %>% filter(gene != "(Intercept)") %>% mutate(dir = sign(coefficient)) %>% dplyr::select(gene, group, dir)
pairs$gene <- gsub(" ","\n", pairs$gene)
edge_list <- data.frame(from = pairs$gene, to = pairs$group)
net <- network(edge_list, directed = FALSE)
net %e% "dir" <- pairs$dir
all_nodes <- unique(c(pairs$gene, pairs$group))
node_types <- ifelse(all_nodes %in% pairs$gene, "Gene", "Group")
net %v% "type" <- node_types
#node_colors <- ifelse(net %v% "type" == "Gene", "grey", "gold")
node_colors <- c(rep("grey",25),"#ddaa33","#bb5566","#004488")
edge_colors <- ifelse(net %e% "dir" == 1, "tomato", "steelblue")
set.seed(1234)
p3 <- ggnet2(net, 
       node.color = node_colors,
       #edge.color = edge_colors,
       edge.size = 1, 
       node.size = 6,
       label = TRUE,
       label.size = 3) +
  theme_void() + 
  labs(x = "", y = "")+
  coord_cartesian(clip = "off") 
p3
write_rds(p3, file = "./output/figures/3d.rds")


a <- ggdraw() + draw_image(image_read_pdf("./output/figures/3a.pdf", pages = 1) %>% image_trim()) 
b <- read_rds("./output/figures/3b.rds")
c <- read_rds("./output/figures/3c.rds")
d <- read_rds("./output/figures/3d.rds")

#plot_grid(plot_grid(a,b, nrow = 1, labels = c("a","b"), rel_widths = c(1,0.8)),
#          plot_grid(c,d, nrow = 1, labels = c("c","d"), rel_widths = c(1,1)), nrow = 2, rel_heights = c(1,0.7))
          
pdf(file = "./output/figures/figure3.pdf", width=14.52, height=7.64)
plot_grid(plot_grid(a, nrow = 1, labels = c("a")),
          plot_grid(plot_grid(b, nrow = 1, labels = c("b")),
                    plot_grid(c,d, nrow = 1, labels = c("c","d"), rel_widths = c(1,1)), nrow = 2))
dev.off()
