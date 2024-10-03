library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(khroma)
library(ggpmisc)
library(cowplot)

#load age acceleration 
results <- read_rds("./output/epigenetic_age_acceleration.rds")

#figure 1b
avg_groups <- results %>% filter(baseline==1) %>% group_by(group, timepoint, category, clock) %>% summarise(mean = mean(age_acc_difference), p = wilcox.test(age_acc_timepoint,age_acc_baseline, paired = TRUE)$p.val)
avg_groups$label <- paste0(avg_groups$group,"_",avg_groups$timepoint)

mat <- reshape2::acast(avg_groups, clock~label, value.var = "mean")
mat_p <- reshape2::acast(avg_groups, clock~label, value.var = "p")
mat_p <- ifelse(mat_p<0.05,"*","")

col_fun <- colorRamp2(seq(-5,5,1),rev(brewer.pal(11, "RdBu")))
lgd = Legend(col_fun = col_fun, title = "foo", at = c(-5, 0, 5), labels = c("< -5", "0", "> 5"))
p0 <- Heatmap(mat,
        column_gap = unit(4, "mm"),
        column_split = factor(sapply(colnames(mat), function(x) strsplit(x, split = "\\_")[[1]][1]), levels = c("TPE + IVIG\n(B)","TPE\n(B)","TPE\n(M)","Sham")),
        column_title_gp = gpar(fontsize = 10),
        rect_gp = gpar(col = "black", lwd = 0.2, lty = 1),
        row_split = sapply(rownames(mat), function(x) avg_groups$category[avg_groups$clock==x] %>% unique),
        row_title_rot = 0,
        border = TRUE,
        use_raster = FALSE,
        row_dend_side = "right",
        row_names_side = "left",
        column_names_side = "top",
        column_names_rot = 0,
        column_labels = sapply(colnames(mat), function(x) strsplit(x, split = "\\_")[[1]][2]),
        cluster_row_slices = TRUE,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        heatmap_legend_param = list(direction = "vertical", title = "Age\nacceleration\ndifference", title_position = "topleft", at = c(-5, 0, 5), labels = c("< -5", "0", "> +5")),
        col = col_fun)
p0
write_rds(p0, file = "./output/figures/1b.rds")

#figure 1c
avg_groups_p <- avg_groups
avg_groups_p$timepoint <- ifelse(avg_groups_p$timepoint=="2","2 vs 1","3 vs 1")
my_comparisons <- list(c("Sham", "TPE\n(M)"), c("Sham", "TPE\n(B)"), c("Sham", "TPE + IVIG\n(B)"))
avg_groups_p %>% group_by(group,timepoint) %>% summarise(mean = mean(mean))

p1 <- ggplot(avg_groups_p, aes(group, mean, group = interaction(group, timepoint), fill = timepoint)) + 
  stat_summary(aes(group = timepoint), fun = mean, geom = 'bar', linewidth = 3, width = 0.75, position = position_dodge(width = 0.8)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', width = 0.2, size = 0.5, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0)+
  stat_compare_means(data = avg_groups_p %>% filter(timepoint=="2 vs 1"), tip.length = 0, comparisons = my_comparisons, label.y = c(-4.9, -5.45, -6), size = 3)+
  stat_compare_means(data = avg_groups_p %>% filter(timepoint=="3 vs 1"), tip.length = 0, comparisons = my_comparisons, label.y = c(1.1, 1.55, 2.0),  size = 3)+
  theme_pubr(border = TRUE)+
  theme(legend.position = c(0.2,0.75))+
  scale_fill_manual(values = brewer.pal(10, "Paired")[1:2])+
  labs(x = "Treatment", y = "Age acceleration difference", fill = "Time point")
p1
write_rds(p1, file = "./output/figures/1c.rds")

#figure 1d
avg_groups_p <- avg_groups_p %>% group_by(category,group,label) %>% summarise(p = t.test(mean, mu = 0)$p.value, mean = mean(mean)) %>% group_by(category) %>% mutate(fdr = p.adjust(p)) %>% dplyr::select(category, group, label, mean, p, fdr)
avg_groups_p$significant <- ifelse(avg_groups_p$fdr<0.05, "*","")
avg_groups_p$timepoint <- sapply(avg_groups_p$label, function(x) strsplit(x, split = "\\_")[[1]][2])
avg_groups_p$category <- gsub("\n"," ",avg_groups_p$category)
avg_groups_p$category <- factor(avg_groups_p$category, levels = c("Adaptation, Causal and Damage clocks","Epigenetic clocks","Fitness Age","PC clocks","Stochastic clocks","Systems Ages"))
avg_groups_p$timepoint <- ifelse(avg_groups_p$timepoint=="2","Time point 2 vs 1","Time point 3 vs 1")

discrete_colors <- color("light")
avg_groups_p$just <- ifelse(avg_groups_p$mean>0,-0.2,1.2)

p2 <- ggplot(avg_groups_p, aes(group, mean, fill = category))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black", lwd = 0.3) +
  geom_text(aes(label = significant, vjust = just), position = position_dodge(width = 0.8), size = 6) +
  scale_fill_manual(values = (discrete_colors(6) %>% as.character())[1:6])+
  theme_pubr(border = TRUE)+
  facet_wrap(.~timepoint)+
  geom_hline(yintercept = 0)+
  theme(legend.position = c(0.75,0.25), legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 8))+
  scale_x_discrete(labels = c("A" = "TPE + IVIG (B)", "B" = "TPE (B)", "D" = "TPE (M)"))+
  labs(x = "Treatment", y = "Age acceleration difference", fill = "")
p2
write_rds(p2, file = "./output/figures/1d.rds")

#figure 1e
pval_permutation <- function(x, y){
  if(nrow(tibble(x,y)%>%na.omit)<3){
    return(NA)
  } else {
    observed_correlation <- cor(x, y, method = "pearson")
    permute_test <- function(x, y, num_permutations = 10000) {
      n <- length(x)
      set.seed(1234)
      permuted_correlations <- replicate(num_permutations, {
        permuted_y <- sample(y)  
        cor(x, permuted_y, method = "pearson")  
      })
      p_value <- mean(abs(permuted_correlations) >= abs(observed_correlation))
      return(p_value)
    }
    p_value <- permute_test(x, y)
    return(p_value)
  }
}

pairwise <- reshape2::dcast(results %>% group_by(patient,group,comparison) %>% summarise(dif = mean(age_acc_difference)), group+patient~comparison, value.var = "dif") %>% na.omit
colnames(pairwise) <- c("group","patient","Time point 2 vs 1","Time point 3 vs 1","Time point 3 vs 2")
pairwise <- pairwise %>% dplyr::group_by(group) %>% mutate(p = unlist(map2(list(`Time point 2 vs 1`),list(`Time point 3 vs 2`),pval_permutation)))

p3 <- ggplot(pairwise, aes(`Time point 2 vs 1`,`Time point 3 vs 2`))+
  geom_point(shape = 21, fill = "#009988")+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  facet_wrap(.~group, nrow = 1, scales = "free_x")+
  stat_correlation(label.x = 0.05, label.y = 0.93) +
  stat_correlation(label.x = 0.05, label.y = 0.83, method = "pearson", aes(label = paste0("p == ", p)))+
  geom_smooth(method = "lm", color = "#009988", alpha = 0.2)+
  geom_ribbon(stat = "smooth",method = "lm",se = TRUE, alpha = 0, colour = "black", linetype = "dotted")+
  theme_pubr(border = TRUE)
write_rds(p3, file = "./output/figures/1e.rds")

#generate final figure
a <- ggdraw() + draw_image(magick::image_read_pdf("./output/figures/1a.pdf", pages = 1)) 
b <- grid.grabExpr(draw(read_rds("./output/figures/1b.rds")))
c <- read_rds("./output/figures/1c.rds")
d <- read_rds("./output/figures/1d.rds")
e <- read_rds("./output/figures/1e.rds")

pdf(file = "./output/figures/figure1.pdf", width=16, height=11)
plot_grid(plot_grid(a,nrow = 1, labels = c("a")),
          plot_grid(
            plot_grid(b,nrow = 1, labels = c("b")),
            plot_grid(plot_grid(c,d,nrow = 1, labels = c("c","d"), rel_widths = c(0.6,1)), plot_grid(e,nrow = 1, labels = c("e")), nrow = 2, rel_heights = c(1,0.9)), ncol = 2, rel_widths = c(0.7,1)), nrow = 2, rel_heights = c(0.6,1))
dev.off()




