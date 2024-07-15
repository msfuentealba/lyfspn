
calculate_correlation_omics <- function(epigenomics,clinical,proteomics,metabolomics,lipidomics,glycomics,cytomics){
  combined <- rbind(epigenomics %>% left_join(reshape2::melt(clinical) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Clinical"),
                    epigenomics %>% left_join(reshape2::melt(proteomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Proteomics"),
                    epigenomics %>% left_join(reshape2::melt(metabolomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Metabolomics"),
                    epigenomics %>% left_join(reshape2::melt(lipidomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Lipidomics"),
                    epigenomics %>% left_join(reshape2::melt(glycomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Glycomics"),
                    epigenomics %>% left_join(reshape2::melt(cytomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Cytomics"))
  
  combined_r <- combined %>% group_by(omics,group,clock,feature) %>% summarise(n = nrow(tibble(dif,value)%>% na.omit), r = cor(dif,value)) %>% filter(n>=3)  
  combined_r
}

stats_correlations_baseline <- function(correlations){
  compare_sham <- correlations %>% dplyr::select(omics,group,clock,feature,r) %>% 
    left_join(correlations %>% filter(group=="C") %>% ungroup %>% dplyr::select(omics,clock,feature,r) %>% set_names("omics","clock","feature","sham_r"))
  compare_sham <- compare_sham %>% group_by(omics,group,feature) %>% na.omit %>% summarise(mean_r = mean(r), mean_sham = mean(sham_r), p = t.test(psych::fisherz(r),mu = 0)$p.value, p_sham = t.test(psych::fisherz(r),psych::fisherz(sham_r))$p.value)
  compare_sham <- compare_sham %>% group_by(omics) %>% summarise(mean_r = mean_r, group = group, feature = feature, fdr = p.adjust(p), fdr_sham = p.adjust(p_sham))
  compare_sham
}
  
extract_significant_baseline <- function(compare_sham){
  compare_sham <- compare_sham %>% filter(feature%in%c(compare_sham %>% filter(group=="A"&fdr<0.05&fdr_sham<0.05) %>% pull(feature)))
  mat <- reshape2::acast(compare_sham %>% filter(feature!="EGFR"), feature~group, value.var = "mean_r")
  mat[is.na(mat)] <- 0
  #mat <- mat[nchar(rownames(mat))<50,]
  mat <- mat[,c("A","B","D","C")]
  colnames(mat) <- c("TPE + IVIG (B)","TPE (B)","TPE (M)","Sham")
  mat
}

find_baseline_links <- function(epigenomics,clinical,proteomics,metabolomics,lipidomics,glycomics,cytomics,compare_sham, mat){
  
  combined <- rbind(epigenomics %>% left_join(reshape2::melt(clinical) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Clinical"),
                    epigenomics %>% left_join(reshape2::melt(proteomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Proteomics"),
                    epigenomics %>% left_join(reshape2::melt(metabolomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Metabolomics"),
                    epigenomics %>% left_join(reshape2::melt(lipidomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Lipidomics"),
                    epigenomics %>% left_join(reshape2::melt(glycomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Glycomics"),
                    epigenomics %>% left_join(reshape2::melt(cytomics) %>% na.omit %>% set_names("feature","patient","value")) %>% mutate(omics = "Cytomics"))
  
  net_matrix <- combined %>% filter(group=="A") %>% group_by(patient,feature) %>% summarise(dif = mean(dif), value = mean(value))
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
  df_link
}

