extract_baseline_clinical <- function(file) {
  biomarker <- readxl::read_xlsx(file, skip = 1)[,1:43]
  biomarker <- biomarker[-118,]
  biomarker$ID[biomarker$ID=="%CD4"] <- "% CD4"
  biomarker <- biomarker %>% filter(!is.na(ID))
  biomarker <- biomarker %>% column_to_rownames(var = "ID")
  biomarker <- biomarker %>% mutate_all(as.numeric)
  biomarker <- biomarker %>% t() %>% t %>% data.frame(check.names = FALSE) %>% rownames_to_column(var = "feature")
  biomarker[biomarker==0] <- NA
  biomarker
}

extract_baseline_proteomics <- function(file){
  protein <- readxl::read_xlsx(file, sheet = 1, skip = 2) 
  protein <- protein[,c(1,271:401)]
  protein <- protein[,!grepl("BLIND|POST",protein[1,] %>% as.vector())]
  protein <- protein %>% as.data.frame() %>% filter(!is.na(...1))
  protein <- protein %>% column_to_rownames(var = "...1")
  colnames(protein) <- sapply(protein[1,] %>% as.vector, function(x) strsplit(x, split = "\\_")[[1]][7])
  colnames(protein) <- paste(sapply(protein[1,] %>% as.character(), function(x) strsplit(x, split = "_")[[1]][6]),
                             ifelse(grepl("^A",colnames(protein)),"A",ifelse(grepl("^B",colnames(protein)),"B",ifelse(grepl("^C",colnames(protein)),"C",ifelse(grepl("^D",colnames(protein)),"D",NA)))),
                             ifelse(grepl("p1$",colnames(protein)),1,ifelse(grepl("p2$",colnames(protein)),2,ifelse(grepl("p3$",colnames(protein)),3,ifelse(grepl("p4$",colnames(protein)),4,NA)))), sep = "_")
  protein <- protein[-1,] 
  protein <- protein %>% mutate_all(as.numeric)
  protein <- log(protein+1)
  protein <- reshape2::melt(protein %>% as.matrix) %>% separate(col = Var2, into = c("patient","group","timepoint"), sep = "_")
  protein$Var1 <- sapply(protein$Var1 %>% as.character(), function(x) strsplit(x, split = ";")[[1]][1])
  protein <- protein %>% filter(timepoint==1) %>% dplyr::select(Var1,patient,group,value) %>% set_names("feature","patient","group","tp1") %>% left_join(protein %>% filter(timepoint==2) %>% dplyr::select(Var1,patient,group,value) %>% set_names("feature","patient","group","tp2"))
  head(protein)
  protein <- protein %>% group_by(patient,feature,group) %>% summarise(tp1 = mean(tp1), tp2 = mean(tp2))
  protein <- reshape2::dcast(protein, feature~patient, value.var = "tp1")
  protein <- protein[!grepl("^KRT", protein$feature),]
  protein
}

extract_baseline_metabolomics <- function(file){
  metabo <- readxl::read_xlsx(file, sheet = 3, skip = 2) 
  metabo <- metabo[-1,]
  metabo$Species <- NULL
  metabo <- metabo[!grepl("^p400 HR",metabo$...1),]
  metabo <- metabo[!grepl("POST",metabo$...1),]
  metabo <- metabo[!grepl("BLIND",metabo$...1),]
  metabo <- metabo[!grepl("1a$",metabo$...1),]
  metabo$...1 <- gsub("HAAS ","HAAS_",metabo$...1)
  metabo$...1 <- gsub(" ","",metabo$...1)
  metabo$...1 <- gsub("__","_",metabo$...1)
  metabo <- metabo[-40,]
  metabo <- metabo[-20,]
  metabo <- metabo[-15,]
  metabo$...1 <- sapply(metabo$...1, function(x) paste(strsplit(x, split = "_")[[1]][c(1,2,4)], collapse = "_"))
  metabo <- metabo %>% column_to_rownames(var = "...1")
  metabo[metabo=="NA"] <- NA
  metabo <- metabo %>% mutate_all(as.numeric)
  metabo[metabo==0] <- NA
  #metabo <- log(metabo+1)
  metabo <- reshape2::melt(metabo %>% as.matrix %>% t()) %>% separate(col = Var2, into = c("patient","group","timepoint"), sep = "_")
  metabo <- metabo %>% filter(timepoint==1) %>% dplyr::select(Var1,patient,group,value) %>% set_names("feature","patient","group","tp1") %>% left_join(metabo %>% filter(timepoint==2) %>% dplyr::select(Var1,patient,group,value) %>% set_names("feature","patient","group","tp2"))
  head(metabo)
  metabo <- metabo %>% group_by(patient,feature,group) %>% summarise(tp1 = mean(tp1), tp2 = mean(tp2))
  metabo <- reshape2::dcast(metabo, feature~patient, value.var = "tp1")
  metabo
}

extract_baseline_lipidomics <- function(file){
  lipid <- readxl::read_xlsx(file, sheet = 1) 
  lipid <- lipid[17:nrow(lipid),]
  lipid$...2 <- NULL
  lipid <- lipid[!grepl("^p400 HR",lipid$...1),]
  lipid <- lipid[!grepl("POST",lipid$...1),]
  lipid <- lipid[!grepl("BLIND",lipid$...1),]
  lipid <- lipid[!grepl("1a$",lipid$...1),]
  lipid$...1 <- gsub("HAAS ","HAAS_",lipid$...1)
  lipid$...1 <- gsub(" ","",lipid$...1)
  lipid$...1 <- gsub("__","_",lipid$...1)
  lipid <- lipid[-38,]
  lipid$...1 <- sapply(lipid$...1, function(x) paste(strsplit(x, split = "_")[[1]][c(1,2,4)], collapse = "_"))
  lipid <- lipid %>% column_to_rownames(var = "...1")
  lipid[lipid=="NA"] <- NA
  lipid <- lipid %>% mutate_all(as.numeric)
  lipid[lipid==0] <- NA
  #lipid <- log(lipid+1)
  lipid <- reshape2::melt(lipid %>% as.matrix %>% t()) %>% separate(col = Var2, into = c("patient","group","timepoint"), sep = "_")
  lipid <- lipid %>% filter(timepoint==1) %>% dplyr::select(Var1,patient,group,value) %>% set_names("feature","patient","group","tp1") %>% left_join(lipid %>% filter(timepoint==2) %>% dplyr::select(Var1,patient,group,value) %>% set_names("feature","patient","group","tp2"))
  head(lipid)
  lipid <- lipid %>% group_by(patient,feature,group) %>% summarise(tp1 = mean(tp1), tp2 = mean(tp2))
  lipid <- reshape2::dcast(lipid, feature~patient, value.var = "tp1")
  lipid
}

extract_baseline_glycomics <- function(file){
  glycan <- readxl::read_xlsx(file, sheet = 1) 
  glycan <- glycan[grepl(".1$|.2$|.3$",glycan$Name),]
  glycan$Name <- paste(glycan$ID,
                       glycan$Group <- gsub("0|1|2|3|4|5|6|7|8|9|10","",glycan$Group),
                       sapply(glycan$Name, function(x) strsplit(x, split = "\\.")[[1]][2]), sep = "_")
  glycan <- glycan[,c(4,79:ncol(glycan))]
  glycan <- glycan %>% column_to_rownames(var = "Name")
  glycan <- glycan %>% mutate_all(as.numeric)
  #glycan <- log(glycan+1)
  glycomics_annotation <- readxl::read_xlsx("./data/glycomics_annotation.xlsx", col_names = FALSE)[,1:2] %>% set_names("feature","name")
  colnames(glycan) <- sapply(colnames(glycan), function(x) glycomics_annotation$name[glycomics_annotation$feature==x])
  glycan <- reshape2::melt(glycan %>% as.matrix %>% t()) %>% separate(col = Var2, into = c("patient","group","timepoint"), sep = "_")
  glycan <- glycan %>% filter(timepoint==1) %>% dplyr::select(Var1,patient,group,value) %>% set_names("feature","patient","group","tp1") %>% left_join(glycan %>% filter(timepoint==2) %>% dplyr::select(Var1,patient,group,value) %>% set_names("feature","patient","group","tp2"))
  head(glycan)
  glycan <- glycan %>% group_by(patient,feature,group) %>% summarise(tp1 = mean(tp1), tp2 = mean(tp2))
  glycan <- reshape2::dcast(glycan, feature~patient, value.var = "tp1")
  glycan
}

extract_baseline_cytomics <- function(file){
  cell <- readxl::read_xlsx(file)
  cell <- cell[c(6,8,11:nrow(cell)),c(9,13:ncol(cell))]
  cell <- cell[c(1:nrow(cell)),cell[2,]%in%c(NA,1,2)]
  cell <- (cell %>% t()) %>% data.frame
  cell <- cell[,c(1,2,4:ncol(cell))]
  cols <- cell[1,3:ncol(cell)] %>% unlist
  cols <- gsub("_"," ",cols)
  cell <- cell[2:nrow(cell),]
  colnames(cell) <- c("patient","timepoint",cols)
  cell$patient <- gsub("_","",cell$patient)
  cell <- cell[,-28] #check
  long_cell <- pivot_longer(cell, 
                            cols = 3:ncol(cell), 
                            names_to = "cell_type", 
                            values_to = "percentage_live")
  long_cell$percentage_live <- as.numeric(long_cell$percentage_live)
  long_cell <- long_cell %>% filter(timepoint==1) %>% dplyr::select(patient,cell_type,percentage_live) %>% set_names("patient","feature","tp1") %>% 
    left_join(long_cell %>% filter(timepoint==2) %>% dplyr::select(patient,cell_type,percentage_live) %>% set_names("patient","feature","tp2"))
  #long_cell$dif <- long_cell$tp2-long_cell$tp1
  long_cell <- reshape2::dcast(long_cell, feature~patient, value.var = "tp1")
  long_cell
}





