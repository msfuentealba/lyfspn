plot_circos_baseline <- function(compare_sham, mat, df_link){
  compare_sham <- compare_sham %>% filter(feature%in%c(compare_sham %>% filter(group=="A"&fdr<0.05&fdr_sham<0.05) %>% pull(feature)))
  
  col_fun1 <- colorRamp2(seq(-1,1,0.2),rev(brewer.pal(11, "RdBu")))
  pdf(file = "./figures/circos_baseline.pdf", height = 7, width = 7)
  circos.par(start.degree = 80, gap.degree = c(5,5,5,5,5,20))
  circos.heatmap(mat, 
                 rownames.side = "outside",
                 dend.side = "none",
                 track.height = 0.2,
                 cell.border = "black",
                 cell.lty = 1,
                 cell.lwd = 0.5,
                 rownames.cex = 0.5,
                 split = sapply(rownames(mat), function(x) unique(compare_sham$omics[compare_sham$feature==x]) %>% na.omit), 
                 show.sector.labels = TRUE,
                 col = col_fun1)
  for(i in seq_len(nrow(df_link))) {
    circos.heatmap.link(df_link$from_index[i],
                        df_link$to_index[i],
                        h.ratio = 0.7,
                        #col = col_fun1(df_link$r[i]),
                        col = ifelse(df_link$r[i]>0.9, "#F4A582","#92C5DE"),
                        #lwd	= 2,
                        rou = 0.5)
  }
  circos.clear()
  dev.off()
}

plot_top_baseline <- function(epigenomics,clinical,proteomics,metabolomics,lipidomics,glycomics,cytomics,compare_sham){
  combined <- rbind(epigenomics %>% left_join(reshape2::melt(clinical) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Clinical"),
                    epigenomics %>% left_join(reshape2::melt(proteomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Proteomics"),
                    epigenomics %>% left_join(reshape2::melt(metabolomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Metabolomics"),
                    epigenomics %>% left_join(reshape2::melt(lipidomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Lipidomics"),
                    epigenomics %>% left_join(reshape2::melt(glycomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Glycomics"),
                    epigenomics %>% left_join(reshape2::melt(cytomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Cytomics"))
  
  compare_sham <- compare_sham %>% filter(feature%in%c(compare_sham %>% filter(group=="A"&fdr<0.05&fdr_sham<0.05) %>% pull(feature)))
  
  top <- compare_sham %>% filter(group=="A"&omics%in%"Clinical") 
  correlations_top <- combined %>% filter(group=="A"&feature%in%top$feature) %>% group_by(patient,feature) %>% summarise(dif = mean(dif), value = mean(value))
  correlations_top$feature <- factor(correlations_top$feature, levels = correlations_top %>% group_by(feature) %>% summarise(r = cor(dif,value)) %>% arrange(r) %>% pull(feature))

  ggplot(correlations_top, aes(dif,value))+
    geom_smooth(method = "lm", color = "#0077bb")+
    geom_point(shape = 21, fill = "#0077bb")+
    #geom_hline(yintercept = 0, lty = 2)+
    #geom_vline(xintercept = 0, lty = 2)+
    #facet_grid2(. ~ feature, scales = "free", independent = "all")+
    facet_wrap(. ~ feature, nrow = 4, scales = "free")+
    #geom_text(aes(label = feature, x = -Inf, y = Inf), hjust = -0.1, vjust = 1.5, color = "black", fontface="bold")+
    stat_correlation(label.x = 0.5, label.y = 0.9) +
    theme_pubr(border = TRUE)+
    labs(x = "Age acceleration difference", y = "Biomarker level")
}

print_plot <- function(plot,name, width, height){
    ggsave(here::here("figures", name), plot, width = 7, height = 7)
}




