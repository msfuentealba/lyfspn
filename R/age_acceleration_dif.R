extract_age_acceleration_dif <- function(file){
  results <- read_csv(file)
  results$Group <- paste0(results$Group,"_",results$Test.Number)
  results$Patient.ID <- results$Label.1
  results <- results[,grepl("EAA|Group|Patient.ID",colnames(results))] 
  results <- reshape2::melt(results) %>% separate(Group, into = c("group","timepoint"), sep = "_") %>% set_names("patient", "group","timepoint","clock","value")
  results$clock <- gsub("\\.EAA","",results$clock)
  results <- results %>% filter(timepoint%in%c(2)) %>% left_join(results %>% filter(timepoint==1) %>% dplyr::select(patient,group,clock,value) %>% set_names("patient","group","clock","baseline")) %>% na.omit
  results$dif <- results$value-results$baseline
  results <- results %>% dplyr::select(patient,group,clock,dif)
}