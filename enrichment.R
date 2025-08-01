cal.age.enrichment <- function(age, cluster){
  # Perform Kruskal-Wallis test for age differences across clusters
  age.test.result = kruskal.test(as.numeric(unlist(age)), as.numeric(unlist(cluster)))
  return(age.test.result)
}

cal.discrete.enrichment <- function(file_dir){
  # Perform chi-square test on discrete variables
  para_tb = table(as.data.frame(file_dir[,2:3]))
  enrichment.result <- chisq.test(para_tb)
  return(enrichment.result)
}

enrich <- function(folder, sign=F, label){
  labels <- label
  age_na <- paste0("data/",folder,"/age.txt")
  
  para_num <- 0
  age <- read.table(age_na, row.names = 1)
  age_test_res = cal.age.enrichment(age, labels)
  p_val = as.matrix(age_test_res)
  p_value = as.data.frame(p_val[3,])
  if(p_value < 0.05){
    para_num <- para_num + 1
  }
  
  para_gen_na <- paste0("data/",folder,"/gender.txt")
  para_gen <- read.table(para_gen_na, sep = "\t", header = T)
  para_gen.pval = cal.discrete.enrichment(cbind(para_gen, labels))
  p_val_gen = as.matrix(para_gen.pval)
  p_value_gen = as.data.frame(p_val_gen[3,])
  if(p_value_gen < 0.05){
    para_num <- para_num + 1
  }

  if(sign == T){
    para_M_na <- paste0("data/",folder,"/pathologic_M.txt")
    para_M <- read.table(para_M_na, sep = "\t", header = T)
    para_M.pval = cal.discrete.enrichment(cbind(para_M, labels))
    p_val_M = as.matrix(para_M.pval)
    p_value_M = as.data.frame(p_val_M[3,])
    if(p_value_M < 0.05){
      para_num <- para_num + 1
    }
    
    para_N_na <- paste0("data/",folder,"/pathologic_N.txt")
    para_N <- read.table(para_N_na, sep = "\t", header = T)
    para_N.pval = cal.discrete.enrichment(cbind(para_N, labels))
    p_val_N = as.matrix(para_N.pval)
    p_value_N = as.data.frame(p_val_N[3,])
    if(p_value_N < 0.05){
      para_num <- para_num + 1
    }
    
    para_T_na <- paste0("data/",folder,"/pathologic_T.txt")
    para_T <- read.table(para_T_na, sep = "\t", header = T)
    para_T.pval = cal.discrete.enrichment(cbind(para_T, labels))
    p_val_T = as.matrix(para_T.pval)
    p_value_T = as.data.frame(p_val_T[3,])
    if(p_value_T < 0.05){
      para_num <- para_num + 1
    }
    
    para_stage_na <- paste0("data/",folder,"/pathologic_stage.txt")
    para_stage <- read.table(para_stage_na, sep = "\t", header = T)
    para_stage.pval = cal.discrete.enrichment(cbind(para_stage, labels))
    p_val_stage = as.matrix(para_stage.pval)
    p_value_stage = as.data.frame(p_val_stage[3,])
    if(p_value_stage < 0.05){
      para_num <- para_num + 1
    }
  }
  return(para_num)
}
