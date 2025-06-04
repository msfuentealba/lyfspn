library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(khroma)
library(ggpmisc)
library(cowplot)
library(doRNG)
library(doParallel)
registerDoParallel(40)
library(exactRankTests)
library(reshape2)

# load age acceleration
results <- read_rds("./data/epigenetic_age_acceleration.rds")
writexl::write_xlsx(results, path = "./output/tables/Sup_Table1.xlsx")


tempa <- foreach(i = 1:length(unique(results$clock)), .combine = rbind) %dorng% {
  temp <- results[results$clock %in% unique(results$clock)[i], ]

  # figure 1b
  avg_groups <- temp %>%
    filter(baseline == 1) %>%
    group_by(group, timepoint, category, clock) %>%
    summarise(mean = mean(age_acc_difference), p.wilcox = wilcox.exact(age_acc_timepoint, age_acc_baseline, paired = T)$p.value, p = cor.test(age_acc_timepoint, age_acc_baseline, method = "spearman")$p.value, r = cor(age_acc_timepoint, age_acc_baseline, method = "spearman"))
  avg_groups
}
avg_groupsa <- tempa
avg_groupsa$label <- paste0(avg_groupsa$group, "_", avg_groupsa$timepoint)

tp2 <- results %>%
  group_by(group, timepoint, category, clock) %>%
  filter(timepoint == 2)

results <- read_rds("./data/epigenetic_age_acceleration.rds") %>% filter(baseline == 1 & timepoint == 2)
biomarker <- read_rds("./data/clinical_baseline.rds")
biomarker$feature[biomarker$feature == "EGFR"] <- "eGFR"
protein <- read_rds("./data/proteomics_baseline.rds")
metabo <- read_rds("./data/metabolomics_baseline.rds")
metabo$feature[metabo$feature == "Creatinine"] <- "Creatinine (M)"
lipid <- read_rds("./data/lipidomics_baseline.rds")
glycan <- read_rds("./data/glycomics_baseline.rds")
cell <- read_rds("./data/cytomics_baseline.rds")
iage <- read_rds("./data/iage_base.rds")

combined <- rbind(
  results %>% left_join(reshape2::melt(biomarker) %>% na.omit() %>% set_names("feature", "patient", "value")) %>% mutate(omics = "Clinical"),
  results %>% left_join(reshape2::melt(protein) %>% na.omit() %>% set_names("feature", "patient", "value")) %>% mutate(omics = "Proteomics"),
  results %>% left_join(reshape2::melt(metabo) %>% na.omit() %>% set_names("feature", "patient", "value")) %>% mutate(omics = "Metabolomics"),
  results %>% left_join(reshape2::melt(lipid) %>% na.omit() %>% set_names("feature", "patient", "value")) %>% mutate(omics = "Lipidomics"),
  results %>% left_join(reshape2::melt(glycan) %>% na.omit() %>% set_names("feature", "patient", "value")) %>% mutate(omics = "Glycomics"),
  results %>% left_join(reshape2::melt(cell) %>% na.omit() %>% set_names("feature", "patient", "value")) %>% mutate(omics = "Cytomics"),
  results %>% left_join(reshape2::melt(iage) %>% na.omit() %>% set_names("feature", "patient", "value")) %>% mutate(omics = "iAge")
)


# combined <- combined %>%
# filter(timepoint == 2) %>%
# left_join(tp2, by = c( "group","timepoint","clock", "category"))


avg_groupsa$p.wilcox.adj <- p.adjust(avg_groupsa$p.wilcox)
mat <- reshape2::acast(avg_groupsa, clock ~ label, value.var = "mean")
mat_p <- reshape2::acast(avg_groupsa, clock ~ label, value.var = "p.wilcox")
mat_p.pergroup <- sapply(as.data.frame(mat_p), p.adjust)
mat_p <- ifelse(mat_p < 0.05, "*", "")

col_fun <- colorRamp2(seq(-5, 5, 1), rev(brewer.pal(11, "RdBu")))
lgd <- Legend(col_fun = col_fun, title = "foo", at = c(-5, 0, 5), labels = c("< -5", "0", "> 5"))
p0 <- Heatmap(mat,
  column_gap = unit(4, "mm"),
  column_split = factor(sapply(colnames(mat), function(x) strsplit(x, split = "\\_")[[1]][1]), levels = c("TPE + IVIG\n(B)", "TPE\n(B)", "TPE\n(M)", "Sham")),
  column_title_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "black", lwd = 0.2, lty = 1),
  row_split = sapply(rownames(mat), function(x) avg_groupsa$category[avg_groupsa$clock == x] %>% unique()),
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
  col = col_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%s", mat_p[i, j]), x, y - height * .25, gp = gpar(fontsize = 20))
  }
)
# p0
write_rds(p0, file = "./output/figures/1b.rds")

# figure 1c

avg_groups_p <- avg_groupsa
avg_groups_p$timepoint <- ifelse(avg_groups_p$timepoint == "2", "2 vs 1", "3 vs 1")
my_comparisons <- list(c("Sham", "TPE\n(M)"), c("Sham", "TPE\n(B)"), c("Sham", "TPE + IVIG\n(B)"))
stats <- rbind(
  data.frame(timepoint = "2 vs 1", group1 = "TPE + IVIG\n(B)", group2 = "Sham", p.val = wilcox.exact(dcast(avg_groups_p, clock ~ group + timepoint, value.var = "mean")[, 2], dcast(avg_groups_p, clock ~ group + timepoint, value.var = "mean")[, 8])$p.value),
  data.frame(timepoint = "2 vs 1", group1 = "TPE\n(B)", group2 = "Sham", p.val = wilcox.exact(dcast(avg_groups_p, clock ~ group + timepoint, value.var = "mean")[, 4], dcast(avg_groups_p, clock ~ group + timepoint, value.var = "mean")[, 8])$p.value),
  data.frame(timepoint = "2 vs 1", group1 = "TPE\n(M)", group2 = "Sham", p.val = wilcox.exact(dcast(avg_groups_p, clock ~ group + timepoint, value.var = "mean")[, 6], dcast(avg_groups_p, clock ~ group + timepoint, value.var = "mean")[, 8])$p.value),
  data.frame(timepoint = "3 vs 1", group1 = "TPE + IVIG\n(B)", group2 = "Sham", p.val = wilcox.exact(dcast(avg_groups_p, clock ~ group + timepoint, value.var = "mean")[, 3], dcast(avg_groups_p, clock ~ group + timepoint, value.var = "mean")[, 9])$p.value),
  data.frame(timepoint = "3 vs 1", group1 = "TPE\n(B)", group2 = "Sham", p.val = wilcox.exact(dcast(avg_groups_p, clock ~ group + timepoint, value.var = "mean")[, 5], dcast(avg_groups_p, clock ~ group + timepoint, value.var = "mean")[, 9])$p.value),
  data.frame(timepoint = "3 vs 1", group1 = "TPE\n(M)", group2 = "Sham", p.val = wilcox.exact(dcast(avg_groups_p, clock ~ group + timepoint, value.var = "mean")[, 7], dcast(avg_groups_p, clock ~ group + timepoint, value.var = "mean")[, 9])$p.value)
)




stats$y.position <- c(-4, -3.5, -3, 4, 3.5, 3)
stats$xmin <- c(.75, 1.75, 2.75, 1.25, 2.25, 3.25)
stats$xmax <- c(3.75, 3.75, 3.75, 4.25, 4.25, 4.25)


stats$p.adj <- round(p.adjust(stats$p.val), 3)


p1 <- ggplot(avg_groups_p, aes(group, mean, group = interaction(group, timepoint), fill = timepoint)) +
  stat_summary(aes(group = timepoint), fun = mean, geom = "bar", linewidth = 3, width = 0.75, position = position_dodge(width = 0.8)) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, size = 0.5, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0) +
  stat_pvalue_manual(stats, label = "p.adj", tip.length = 0, size = 3) +

  #  stat_compare_means(data = avg_groups_p %>% filter(timepoint=="2 vs 1"), tip.length = 0,
  # comparisons = my_comparisons, label.y = c(-4.9, -5.45, -6), size = 3)+
  #  stat_compare_means(data = avg_groups_p %>% filter(timepoint=="3 vs 1"), tip.length = 0, #comparisons = my_comparisons, label.y = c(1.1, 1.55, 2.0),  size = 3)+
  theme_pubr(border = TRUE) +
  theme(legend.position = c(0.2, 0.75)) +
  scale_fill_manual(values = brewer.pal(10, "Paired")[1:2]) +
  labs(x = "Treatment", y = "Age acceleration difference", fill = "Time point")

write_rds(p1, file = "./output/figures/1c.rds")

# figure 1d
avg_groups_p <- avg_groupsa %>%
  group_by(category, group, label) %>%
  summarise(
    p.wilcox = wilcox.exact(mean, mu = 0)$p.value,
    mean = mean(mean)
  ) %>%
  group_by(category)
avg_groups_p$fdr.wilcox <- p.adjust(avg_groups_p$p.wilcox)

avg_groups_p$significant <- ifelse(avg_groups_p$fdr.wilcox < 0.05, "*", "")
avg_groups_p$timepoint <- sapply(avg_groups_p$label, function(x) strsplit(x, split = "\\_")[[1]][2])
avg_groups_p$category <- gsub("\n", " ", avg_groups_p$category)
avg_groups_p$category <- factor(avg_groups_p$category, levels = c("Adaptation, Causal and Damage clocks", "Epigenetic clocks", "Fitness Age", "PC clocks", "Stochastic clocks", "Systems Ages"))
avg_groups_p$timepoint <- ifelse(avg_groups_p$timepoint == "2", "Time point 2 vs 1", "Time point 3 vs 1")
library(Seurat.utils)
discrete_colors <- colors()[grep("light", colors())][c(1, 26, 12, 21, 32, 31)] # color("light")
avg_groups_p$just <- ifelse(avg_groups_p$mean > 0, -0.2, 1.2)

p2 <- ggplot(avg_groups_p, aes(group, mean, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black", lwd = 0.3) +
  geom_text(aes(label = significant, vjust = just), position = position_dodge(width = 0.8), size = 6) +
  scale_fill_manual(values = (discrete_colors %>% as.character())[1:6]) +
  theme_pubr(border = TRUE) +
  facet_wrap(. ~ timepoint, ncol = 1) +
  geom_hline(yintercept = 0) +
  theme(legend.position = c(0.75, 0.1), legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 8)) +
  scale_x_discrete(labels = c("A" = "TPE + IVIG (B)", "B" = "TPE (B)", "D" = "TPE (M)")) +
  labs(x = "Treatment", y = "Age acceleration difference", fill = "")

write_rds(p2, file = "./output/figures/1d.rds")
# figure 1e

results <- read_rds("./data/epigenetic_age_acceleration.rds")
library(jmuOutlier)
pairwise <- reshape2::dcast(results %>% group_by(patient, group, comparison) %>% summarise(dif = mean(age_acc_difference)), group + patient ~ comparison, value.var = "dif") %>% na.omit()
colnames(pairwise) <- c("group", "patient", "Time point 2 vs 1", "Time point 3 vs 1", "Time point 3 vs 2")
pairwise <- pairwise %>% left_join(pairwise %>% dplyr::group_by(group) %>% summarise(cor = cor(`Time point 2 vs 1`, `Time point 3 vs 2`, method = "spearman"), p = perm.cor.test(`Time point 2 vs 1`, `Time point 3 vs 2`, method = "spearman")$p.value))
pairwise$p <- format(pairwise$p, scientific = TRUE, digits = 3)

p3 <- ggplot(pairwise, aes(`Time point 2 vs 1`, `Time point 3 vs 2`)) +
  geom_point(shape = 21, fill = "#009988") +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  facet_wrap(. ~ group, nrow = 1, scales = "free_x") +
  stat_correlation(label.x = 0.05, label.y = 0.93, method = "spearman", aes(label = paste0("R == ", round(cor, 2)))) +
  stat_correlation(label.x = 0.05, label.y = 0.83, method = "spearman", aes(label = paste0("p = ", p)), parse = FALSE) +
  geom_smooth(method = "lm", color = "#009988", alpha = 0.2) +
  geom_ribbon(stat = "smooth", method = "lm", se = TRUE, alpha = 0, colour = "black", linetype = "dotted") +
  theme_pubr(border = TRUE)

b <- grid.grabExpr(draw(read_rds("./output/figures/1b.rds")))
c <- read_rds("./output/figures/1c.rds")
d <- read_rds("./output/figures/1d.rds")


pdf(file = "./output/figures/figure1.pdf", width = 16, height = 11)

plot_grid(
  plot_grid(b, nrow = 1, labels = c("a")),
  plot_grid(plot_grid(c, d, nrow = 2, labels = c("b", "c"), rel_heights = c(0.6, 1))),
  ncol = 2, rel_widths = c(1, .9)
)
dev.off()
