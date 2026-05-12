p_value_assign <- function(p_value_in, custom_lims = c(0.05, 0.01, 0.001, 0.0001)) {
  if (p_value_in < custom_lims[4]) {
    p_out = as.character("****")
  } else if ((p_value_in > custom_lims[4]) &  (p_value_in < custom_lims[3])) {
    p_out = as.character("***")
  } else if ((p_value_in > custom_lims[3]) &  (p_value_in < custom_lims[2])) {
    p_out = as.character("**")
  } else if ((p_value_in > custom_lims[2]) &  (p_value_in < custom_lims[1])) {
    p_out = as.character("*")
  } else if (p_value_in > custom_lims[1]) {
    p_out = as.character("ns")
  }
  return(p_out)
}

clean_filename <- function(x) {
  x = gsub("'", "", x)
  x = gsub("[^A-Za-z0-9]+", "_", x)
  x = gsub("^_|_$", "", x)
  return(x)
}

U_Mann_Whitney_peptide_test <- function(my_data, grouping_data) {
  out_test = data.frame(matrix(data = NA, ncol = 3, nrow = nrow(my_data)))
  colnames(out_test) = c("Peptide", "p_val", "signif")
  for (i_test in 1:nrow(my_data)) {
    my_data_tmp = data.frame(matrix(data = NA, ncol = 2, nrow = length(grouping_data$Condition)))
    colnames(my_data_tmp) = c("Group", "Value")
    my_data_tmp$Group = grouping_data$Condition
    my_data_tmp$Value = as.numeric(my_data[i_test, 2:ncol(my_data)])
    
    p_val_tmp = wilcox.test(Value ~ Group, data=my_data_tmp, exact = FALSE)[["p.value"]]
    
    if (is.na(p_val_tmp)) {
      p_val_tmp = 1
    }
    
    out_test$Peptide[i_test] = my_data$Peptide[i_test]
    out_test$p_val[i_test] = p_val_tmp
    out_test$signif[i_test] = p_value_assign(p_val_tmp)
  }
  return(out_test)
}

my_data_input_prepare <- function(my_data, grouping_data, n_first, n_last, 
                                  rm_precursors, rm_zero_na, rm_blanks) {
  
  if (colnames(my_data)[1] == "Protein") {
    colnames(my_data)[1] = "Protein Name"
  }

  # Remove precursors
  if (rm_precursors == TRUE) {
    my_data = my_data[-(which(my_data$`Fragment Ion` == "precursor")),]
    my_data = my_data[-(which(my_data$`Fragment Ion` == "precursor [M+1]")),]
    my_data = my_data[-(which(my_data$`Fragment Ion` == "precursor [M+2]")),]
  }
  
  # Remove zeros and NAs
  if (rm_zero_na == TRUE) {
    my_data[my_data == 0] = NA
    my_data[my_data == "#N/A"] = NA
  }
  
  # Remove blank samples
  if (rm_blanks == TRUE) {
    tb_rm_blanks = which(colnames(my_data) %in% grouping_data$Replicate[which(grouping_data$Condition == 0)])
    if (sum(which(colnames(my_data) %in% grouping_data$Replicate[which(grouping_data$Condition == 0)])) > 0) {
      my_data = my_data[,-tb_rm_blanks]
      grouping_data = grouping_data[-which(grouping_data$Condition == 0),]
    }
  }
  
  
  my_data_additional = my_data[, -which(colnames(my_data) %in% grouping_data$Replicate)]
  my_data = my_data[, which(colnames(my_data) %in% grouping_data$Replicate)]
  
  my_data = cbind(my_data_additional$`Protein Name`,  my_data)
  my_data = as.data.frame(t(my_data))
  colnames(my_data) = as.character(my_data[1,])
  my_data = my_data[-1,]
  my_data = my_data %>% janitor::clean_names()
  my_data_additional$`Cleaned Protein Name` = colnames(my_data)
  metabolite_names = colnames(my_data)
  my_data$Group = NaN
  tb_removed = c()
  for (i_sample in 1:nrow(my_data)) {
    if (is.na(grouping_data$Condition[which(grouping_data$Replicate == rownames(my_data)[i_sample])])) {
      tb_removed = c(tb_removed, i_sample)
    } else {
      my_data$Group[i_sample] = grouping_data$Condition[which(grouping_data$Replicate == rownames(my_data)[i_sample])]
    }
  }
  
  if (!is.null(tb_removed)) {
    my_data = my_data[-tb_removed,]
  }
  labels = as.character(my_data$Group)
  
  
  return(list(my_data, my_data_additional, grouping_data, labels, metabolite_names))
}

get_grouping_data <- function(my_data) {
  groups_tbd = colnames(my_data)[10:length(colnames(my_data))]
  split_list = strsplit(groups_tbd, "_")
  Condition = sapply(split_list, function(x) {
    nums = x[grepl("^[0-9]+$", x)]
    if (length(nums) == 0) return(NA)
    as.numeric(paste(nums, collapse = "."))
  })
  
  grouping_data = data.frame(
    Replicate = groups_tbd,
    Condition = Condition
  )
  
  return(grouping_data)
}

MSstats_prepare_own <- function(my_data, grouping_data, level = "Protein") {
  out_data = data.frame(matrix(ncol = 10, nrow = 0))
  colnames(out_data) = c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge", "IsotopeLabelType", 
                         "Condition", "BioReplicate", "Run", "Intensity")
  for (i_row in 1:nrow(my_data)) {
    for (i_column in 11:ncol(my_data)) {
      out_data[nrow(out_data)+1,] = NA
      i_current_out_data = nrow(out_data)
      
      if (level == "Protein") {
        out_data$ProteinName[i_current_out_data] = my_data$`Protein Name`[i_row]
      } else if (level == "Peptide") {
        out_data$ProteinName[i_current_out_data] = my_data$Peptide[i_row]
      } else {
        print("Error at MSstats data preparation")
        return(NULL)
      }
      
      out_data$PeptideSequence[i_current_out_data] = my_data$Peptide[i_row]
      out_data$PrecursorCharge[i_current_out_data] = my_data$`Precursor Charge`[i_row]
      out_data$FragmentIon[i_current_out_data] = my_data$`Fragment Ion`[i_row]
      out_data$ProductCharge[i_current_out_data] = my_data$`Product Charge`[i_row]
      out_data$IsotopeLabelType[i_current_out_data] = "L"
      out_data$Condition[i_current_out_data] = grouping_data$Condition[which(grouping_data$Replicate == colnames(my_data)[i_column])]
      out_data$BioReplicate[i_current_out_data] = substr(colnames(my_data)[i_column], nchar(colnames(my_data)[i_column]), nchar(colnames(my_data)[i_column]))
      out_data$Run[i_current_out_data] = colnames(my_data)[i_column]
      out_data$Intensity[i_current_out_data] = my_data[i_row, i_column]
    }
  }
  return(out_data)
}

MSstats_return_protein_own <- function(my_data_new_MSstats_return) {
  out_data = data.frame(matrix(ncol = (1+length(levels(my_data_new_MSstats_return$ProteinLevelData$originalRUN))),
                               nrow = (length(levels(my_data_new_MSstats_return$ProteinLevelData$Protein)))))
  colnames(out_data) = c("Protein", as.character(unique(my_data_new_MSstats_return$ProteinLevelData$originalRUN)))
  out_data$Protein = as.character(unique(my_data_new_MSstats_return$ProteinLevelData$Protein))
  for (i_protein in 1:nrow(out_data)) {
    tmp_data = my_data_new_MSstats_return$ProteinLevelData[which(my_data_new_MSstats_return$ProteinLevelData$Protein == out_data$Protein[i_protein]), ]
    out_data[i_protein, 2:ncol(out_data)] = 2^tmp_data$LogIntensities
  }
  return(out_data)
}

MSstats_return_feature_own <- function(my_data_new_MSstats_return, my_data) {
  out_data = my_data
  out_data[, 11:ncol(out_data)] = NA
  tmp_data = my_data_new_MSstats_return$FeatureLevelData
  for (i_row in 1:nrow(tmp_data)) {
    tmp_data_row = tmp_data[i_row, ]
    protein_set = which(out_data$`Protein Name` == as.character(tmp_data_row$PROTEIN))
    peptide_set = which(out_data$Peptide == strsplit(as.character(tmp_data_row$PEPTIDE), "_")[[1]][1])
    precursor_charge_set = which(out_data$`Precursor Charge` == strsplit(as.character(tmp_data_row$PEPTIDE), "_")[[1]][2])
    fragment_set = which(out_data$'Fragment Ion' == strsplit(as.character(tmp_data_row$TRANSITION), "_")[[1]][1])
    fragment_charge_set = which(out_data$`Product Charge` == strsplit(as.character(tmp_data_row$TRANSITION), "_")[[1]][2])
    row_set = intersect(intersect(intersect(intersect(protein_set,peptide_set),precursor_charge_set), fragment_set), fragment_charge_set)
    column_set = which(colnames(out_data) == tmp_data_row$originalRUN)
    out_data[row_set, column_set] = 2^tmp_data_row$newABUNDANCE
  }
  return(out_data)
}

MDI_HCA_removal <- function(my_data, labels, metabolite_names, my_data_additional, grouping_data, 
                            print_MD_procentage = FALSE, print_frag_HCA = FALSE, print_sample_HCA = FALSE) {
  library(pheatmap, quietly=T, warn.conflicts=F)
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(dendextend))
  
  df = as.data.frame(!is.na(my_data[,1:(ncol(my_data)-1)])) * 1
  
  df_tmp = as.matrix(df)
  dend = as.dendrogram(hclust(dist(t(df_tmp))))
  dend_test = color_branches(dend, k = 2)
  dend_labels = dend_test %>% labels
  
  # HCA for ion fragments (i.e. attribiutes)
  no_labels = (dend_test %>% get_nodes_attr("members", 2))
  if (no_labels < (ncol(df_tmp) - no_labels)) {
    dend_labels = dend_labels[1:no_labels]
    dend = color_branches(dend, k = 2, col = c("brown1", "cornflowerblue"))
  } else {
    dend_labels = dend_labels[no_labels:ncol(df_tmp)]
    dend = color_branches(dend, k = 2, col = c("cornflowerblue", "brown1"))
  }
  tb_rm_names = dend_labels
  
  if (print_frag_HCA == TRUE) {
    col_fun = colorRamp2(c(0, 1), c("brown1", "cornflowerblue"))
    cp = ComplexHeatmap::Heatmap(df_tmp, cluster_rows = F,  show_column_names = F,
                                 show_heatmap_legend = F, row_names_gp = grid::gpar(fontsize = 8),
                                 column_dend_height = unit(4, "cm"), column_title=NULL, col = col_fun,
                                 column_split = 2, cluster_columns = dend)
    draw(cp, column_title = "Clustering for ion fragments (i.e. attributes)", column_title_gp=grid::gpar(fontsize=16))
  }
  
  
  # HCA for samples
  dend = as.dendrogram(hclust(dist(df_tmp)))
  dend_test = color_branches(dend, k = 2)
  dend_labels = dend_test %>% labels
  
  no_labels = (dend_test %>% get_nodes_attr("members", 2))
  if (no_labels < (ncol(df_tmp) - no_labels)) {
    dend_labels = dend_labels[1:no_labels]
    dend = color_branches(dend, k = 2, col = c("brown1", "cornflowerblue"))
  } else {
    dend_labels = dend_labels[no_labels:ncol(df_tmp)]
    dend = color_branches(dend, k = 2, col = c("cornflowerblue", "brown1"))
  }
  
  if (print_sample_HCA == TRUE) {
    col_fun = colorRamp2(c(0, 1), c("brown1", "cornflowerblue"))
    cp = ComplexHeatmap::Heatmap(df_tmp, cluster_rows = dend, cluster_columns = F, show_column_names = F,
                                 show_heatmap_legend = F, row_names_gp = grid::gpar(fontsize = 8),
                                 column_dend_height = unit(4, "cm"), row_split = 2, row_title=NULL, column_title=NULL,
                                 col = col_fun)
    draw(cp, column_title = "Clustering for samples", column_title_gp=grid::gpar(fontsize=16))
  }
  
  group_row = as.data.frame(t(my_data))[nrow(as.data.frame(t(my_data))),]
  metabolite_names = metabolite_names[-(which(metabolite_names %in% tb_rm_names))]
  my_data_core = as.data.frame(t(my_data))[-nrow(as.data.frame(t(my_data))),]
  my_data_core = mutate_all(my_data_core, function(x) as.numeric(as.character(x)))
  my_data = cbind(my_data_additional, my_data_core)
  my_data_additional = my_data_additional[-(which(my_data_additional$`Cleaned Protein Name` %in% tb_rm_names)),]
  my_data = my_data[-(which(my_data$`Cleaned Protein Name` %in% tb_rm_names)),]
  
  if (print_MD_procentage == TRUE) {
    miss_data_proc = (sum(is.na(my_data)) / (nrow(my_data) * length(which(colnames(my_data) %in% grouping_data$Replicate)))) * 100
    print(paste("Missing data (step 1 HCA): ",  round(miss_data_proc,3), "%", sep=""))
  }
  
  return(my_data)
}

MDI_kNN <- function(my_data, grouping_data, set_kNN_k_parameter) {
  library(multiUS, quietly=T, warn.conflicts=F)
  
  my_data_core_before = my_data[, which(colnames(my_data) %in% grouping_data$Replicate)]
  
  my_data_new_additional = my_data[, -(which(colnames(my_data) %in% grouping_data$Replicate))]
  my_data_new_core = my_data[, which(colnames(my_data) %in% grouping_data$Replicate)]
  my_data_new_core = KNNimp(my_data_new_core, k = set_kNN_k_parameter, scale = TRUE, meth = "weighAvg")
  
  my_data_new = cbind(my_data_new_additional, my_data_new_core)
  
  changed_values = is.na(my_data_core_before) * my_data_new_core
  changed_values[changed_values == 0] = NA
  
  return_list = list("my_data_new" = my_data_new, "changed_values" = changed_values)
  return(return_list)
}

MDI_data_driven <- function(my_data, grouping_data,
                            print_MD_procentage = FALSE, show_pb = TRUE) {
  
  library(progress, quietly=T, warn.conflicts=F)
  library(EnvStats, quietly=T, warn.conflicts=F)
  
  my_data_core_before = my_data[, which(colnames(my_data) %in% grouping_data$Replicate)]
  
  my_data_cleaned = my_data
  my_data_cleaned$Fragment = paste(my_data_cleaned$`Fragment Ion`, "+", my_data_cleaned$`Product Charge`, sep="")
  all_intensity_ratios = data.frame(matrix(NaN, nrow = 0, ncol = ncol(my_data_cleaned)+1))
  colnames(all_intensity_ratios)[1:4] = c("Protein_name","Peptide_name","Ratio_name","CV")
  colnames(all_intensity_ratios)[5:ncol(all_intensity_ratios)] = colnames(my_data_cleaned)[4:ncol(my_data_cleaned)]
  
  if (show_pb == TRUE) {
    pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                          total = (2*length(unique(my_data_cleaned$Peptide))),
                          complete = "=",   # Completion bar character
                          incomplete = "-", # Incomplete bar character
                          current = ">",    # Current bar character
                          clear = FALSE,    # If TRUE, clears the bar when finish
                          width = 100)      # Width of the progress bar
  }
  
  for (i_peptide in 1:length(unique(my_data_cleaned$Peptide))) {
    if (show_pb == TRUE) {
      pb$tick()
    }
    
    peptide_tmp = unique(my_data_cleaned$Peptide)[i_peptide]
    my_data_tmp = my_data_cleaned[my_data_cleaned$Peptide == peptide_tmp, ]
    fragments_unique = my_data_tmp$Fragment
    my_data_tmp_additional = my_data_tmp[, -(which(colnames(my_data_tmp) %in% grouping_data$Replicate))]
    protein_tmp = my_data_tmp_additional$`Protein Name`[1]
    my_data_tmp = my_data_tmp[, which(colnames(my_data_tmp) %in% grouping_data$Replicate)]
    my_data_tmp = cbind(as.factor(fragments_unique), my_data_tmp)
    colnames(my_data_tmp)[1] = "Fragment"
    
    if (nrow(my_data_tmp) > 1) {
      all_combinations = combn(as.character(my_data_tmp$Fragment),2)
      
      out_data = data.frame(matrix(NaN, nrow = ncol(all_combinations), ncol = ncol(my_data_tmp)+3))
      colnames(out_data) = c("Protein_name","Peptide","Ratio_name","CV", colnames(my_data_tmp)[2:ncol(my_data_tmp)])
      for (i2 in 1:ncol(all_combinations)) {
        frag_tmp_1 = all_combinations[1,i2]
        frag_tmp_2 = all_combinations[2,i2]
        
        data_tmp_1 = my_data_tmp[which(my_data_tmp == frag_tmp_1), -1]
        data_tmp_2 = my_data_tmp[which(my_data_tmp == frag_tmp_2), -1]
        
        for (i_ratio in 1:(ncol(my_data_tmp)-1)) {
          if (!is.na(data_tmp_1[i_ratio]) & !is.na(data_tmp_2[i_ratio])) {
            out_data[i2,i_ratio+4] = as.numeric(data_tmp_1[i_ratio]) / as.numeric(data_tmp_2[i_ratio])
          } else {
            out_data[i2,i_ratio+4] = NA
          }
        }
        out_data$Ratio_name[i2] = paste(frag_tmp_1 , " / ", frag_tmp_2, sep="")
        out_data$CV[i2] = cv(as.numeric(out_data[i2, (5:ncol(out_data))]), na.rm = TRUE)
        out_data$Protein_name[i2] = protein_tmp
        out_data$Peptide[i2] = peptide_tmp
      }
      all_intensity_ratios = rbind(all_intensity_ratios, out_data)
    }
  }
  rm(my_data_cleaned, data_tmp_1, data_tmp_2, i_peptide, out_data)
  
  # Imput missing data based on the ratios
  my_data_new_additional = my_data[, -(which(colnames(my_data) %in% grouping_data$Replicate))]
  my_data_new_core = my_data[, which(colnames(my_data) %in% grouping_data$Replicate)]
  my_data_new_core[,] = NA
  
  used_intensity_ratios = data.frame(matrix(NaN, nrow = 0, ncol = 3))
  colnames(used_intensity_ratios) = c("Peptide_name","frag_1", "frag_1")
  
  for (i_peptide in 1:length(unique(my_data$Peptide))) {
    if (show_pb == TRUE) {
      pb$tick()
    }
    
    peptide_tmp = unique(my_data$Peptide)[i_peptide]
    my_data_tmp = my_data[my_data$Peptide == peptide_tmp, ]
    data_imputation_rows = which((my_data_new_additional$Peptide == peptide_tmp) == TRUE)
    fragments_unique = my_data_tmp$`Fragment Ion`
    my_data_tmp_additional = my_data_tmp[, -(which(colnames(my_data_tmp) %in% grouping_data$Replicate))]
    my_data_tmp = my_data_tmp[, which(colnames(my_data_tmp) %in% grouping_data$Replicate)]
    my_data_tmp$Fragment = paste(as.character(fragments_unique), "+", my_data_tmp_additional$`Product Charge`, sep="")
    fragments_unique = my_data_tmp$Fragment
    my_data_tmp = cbind(as.factor(fragments_unique), my_data_tmp)
    colnames(my_data_tmp)[1] = "Fragment"
    my_data_tmp = my_data_tmp[,-ncol(my_data_tmp)]
    
    if (nrow(my_data_tmp) <= 1) {
      my_data_new_core = my_data_new_core[-(which((my_data_new_core$Peptide == peptide_tmp) == TRUE)), ]
      my_data_new_additional = my_data_new_additional[-(which((my_data_new_additional$Peptide == peptide_tmp) == TRUE)), ]
    } else {
      
      for (i_sample in 2:ncol(my_data_tmp)) {
        sample_tmp = my_data_tmp[,i_sample]
        fragments_tmp = as.character(my_data_tmp[,1])
        
        if (sum(is.na(my_data_tmp[,i_sample])) == 0) {
          my_data_new_core[data_imputation_rows, i_sample-1] = sample_tmp
        } else {
          fragments_combinations_tmp = combn(fragments_tmp,2)
          
          for (i_row in 1:length(sample_tmp)) {
            sample_NA_tmp = sample_tmp[i_row]
            fragments_NA_tmp = fragments_tmp[i_row]
            data_imputation_rows_tmp = data_imputation_rows[1] + i_row - 1
            
            if (!is.na(sample_NA_tmp)) {
              my_data_new_core[data_imputation_rows_tmp, i_sample-1] = sample_NA_tmp
            } else {
              
              comb_inx = c((which(fragments_NA_tmp == fragments_combinations_tmp[1,])),
                           (which(fragments_NA_tmp == fragments_combinations_tmp[2,])))
              comb_inx = fragments_combinations_tmp[,comb_inx]
              if (!is.null(ncol(comb_inx))) {
                comb_inx = paste(comb_inx[1,], " / ", comb_inx[2,], sep="")
              } else {
                comb_inx = paste(comb_inx[1], " / ", comb_inx[2], sep="")
              }
              
              tmp_intensity_ratios = all_intensity_ratios[which(peptide_tmp == all_intensity_ratios$Peptide),]
              
              if (nrow(tmp_intensity_ratios) == 0) {
                my_data_new_core[data_imputation_rows_tmp, i_sample-1] = sample_NA_tmp
              } else {
                
                comb_inx_matched = c()
                for (i_ratios_tmp in 1:nrow(tmp_intensity_ratios)) {
                  comb_inx_matched = c(comb_inx_matched, which(comb_inx[i_ratios_tmp] == tmp_intensity_ratios$Ratio_name))
                }
                if (length(comb_inx_matched) == 0) {
                  my_data_new_core[data_imputation_rows_tmp, i_sample-1] = sample_NA_tmp
                } else {
                  tmp_intensity_ratios = tmp_intensity_ratios[comb_inx_matched, ]
                  tmp_intensity_ratios = tmp_intensity_ratios[order(tmp_intensity_ratios$CV, decreasing = FALSE),]
                  fragmet_ratio_inx = 1
                  fragmet_ratio_flag = 0
                  while (fragmet_ratio_flag == 0) {
                    if (fragmet_ratio_inx > nrow(tmp_intensity_ratios)) {
                      my_data_new_core[data_imputation_rows_tmp, i_sample-1] = NA
                      break
                    }
                    
                    tmp_intensity_ratios_min_CV = tmp_intensity_ratios[fragmet_ratio_inx, ]
                    matched_fragments = strsplit(tmp_intensity_ratios_min_CV$Ratio_name, split = " / ")[[1]]
                    
                    used_intensity_ratios[nrow(used_intensity_ratios)+1,1] = peptide_tmp
                    used_intensity_ratios[nrow(used_intensity_ratios),2] = matched_fragments[1]
                    used_intensity_ratios[nrow(used_intensity_ratios),3] = matched_fragments[2]
                    
                    if (which(fragments_NA_tmp == matched_fragments) == 1) {
                      if (!is.na(my_data_tmp[which(matched_fragments[2] == my_data_tmp$Fragment), i_sample])) {
                        median_ratio = median(as.numeric(tmp_intensity_ratios_min_CV[,5:ncol(tmp_intensity_ratios_min_CV)]), na.rm = TRUE)
                        my_data_new_core[data_imputation_rows_tmp, i_sample-1] = my_data_tmp[which(matched_fragments[2] == my_data_tmp$Fragment), i_sample] * median_ratio
                        fragmet_ratio_flag = 1
                      } else {
                        fragmet_ratio_inx = fragmet_ratio_inx + 1
                      }
                    } else if (which(fragments_NA_tmp == matched_fragments) == 2) {
                      if (!is.na(my_data_tmp[which(matched_fragments[1] == my_data_tmp$Fragment), i_sample])) {
                        median_ratio = median(as.numeric(tmp_intensity_ratios_min_CV[,5:ncol(tmp_intensity_ratios_min_CV)]), na.rm = TRUE)
                        my_data_new_core[data_imputation_rows_tmp, i_sample-1] = my_data_tmp[which(matched_fragments[1] == my_data_tmp$Fragment), i_sample] / median_ratio
                        fragmet_ratio_flag = 1
                      } else {
                        fragmet_ratio_inx = fragmet_ratio_inx + 1
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  my_data_new = cbind(my_data_new_additional, my_data_new_core)
  
  changed_values = is.na(my_data_core_before) * my_data_new_core
  changed_values[changed_values == 0] = NA
  
  if (print_MD_procentage == TRUE) {
    miss_data_proc = (sum(is.na(my_data_new)) / (nrow(my_data_new) * length(which(colnames(my_data_new) %in% grouping_data$Replicate)))) * 100
    print(paste("Missing data (step 2 Fragment ratios): ",  round(miss_data_proc,3), "%", sep=""))
  }
  
  return_list = list("my_data_new" = my_data_new, "all_intensity_ratios" = all_intensity_ratios, 
                     "changed_values" = changed_values, "used_intensity_ratios" = used_intensity_ratios)
  return(return_list)
  
}

MDI_AI <- function(my_data, grouping_data, choose_ref = "CVs",
                   set_Koina_model = "ms2pip_timsTOF2024",
                   print_MD_procentage = FALSE) {
  
  # Remove doubly charged fragments
  if (length(which(my_data$`Product Charge` == 2)) != 0) {
    my_data = my_data[-which(my_data$`Product Charge` == 2),]
    my_data_additional = my_data_additional[-which(my_data_additional$`Product Charge` == 2),]
  }
  
  # Get CVs from data (and ratios)
  my_data_cleaned = my_data
  my_data_cleaned$Fragment = paste(my_data_cleaned$`Fragment Ion`, "+", my_data_cleaned$`Product Charge`, sep="")
  all_intensity_ratios_DD = data.frame(matrix(NaN, nrow = 0, ncol = ncol(my_data_cleaned)+1))
  colnames(all_intensity_ratios_DD)[1:4] = c("Protein_name","Peptide_name","Ratio_name","CV")
  colnames(all_intensity_ratios_DD)[5:ncol(all_intensity_ratios_DD)] = colnames(my_data_cleaned)[4:ncol(my_data_cleaned)]
  
  for (i_peptide in 1:length(unique(my_data_cleaned$Peptide))) {
    peptide_tmp = unique(my_data_cleaned$Peptide)[i_peptide]
    my_data_tmp = my_data_cleaned[my_data_cleaned$Peptide == peptide_tmp, ]
    fragments_unique = my_data_tmp$Fragment
    my_data_tmp_additional = my_data_tmp[, -(which(colnames(my_data_tmp) %in% grouping_data$Replicate))]
    protein_tmp = my_data_tmp_additional$`Protein Name`[1]
    my_data_tmp = my_data_tmp[, which(colnames(my_data_tmp) %in% grouping_data$Replicate)]
    my_data_tmp = cbind(as.factor(fragments_unique), my_data_tmp)
    colnames(my_data_tmp)[1] = "Fragment"
    
    if (nrow(my_data_tmp) > 1) {
      all_combinations = combn(as.character(my_data_tmp$Fragment),2)
      
      out_data = data.frame(matrix(NaN, nrow = ncol(all_combinations), ncol = ncol(my_data_tmp)+3))
      colnames(out_data) = c("Protein_name","Peptide","Ratio_name","CV", colnames(my_data_tmp)[2:ncol(my_data_tmp)])
      for (i2 in 1:ncol(all_combinations)) {
        frag_tmp_1 = all_combinations[1,i2]
        frag_tmp_2 = all_combinations[2,i2]
        
        data_tmp_1 = my_data_tmp[which(my_data_tmp == frag_tmp_1), -1]
        data_tmp_2 = my_data_tmp[which(my_data_tmp == frag_tmp_2), -1]
        
        for (i_ratio in 1:(ncol(my_data_tmp)-1)) {
          if (!is.na(data_tmp_1[i_ratio]) & !is.na(data_tmp_2[i_ratio])) {
            out_data[i2,i_ratio+4] = as.numeric(data_tmp_1[i_ratio]) / as.numeric(data_tmp_2[i_ratio])
          } else {
            out_data[i2,i_ratio+4] = NA
          }
        }
        out_data$Ratio_name[i2] = paste(frag_tmp_1 , " / ", frag_tmp_2, sep="")
        out_data$CV[i2] = cv(as.numeric(out_data[i2, (5:ncol(out_data))]), na.rm = TRUE)
        out_data$Protein_name[i2] = protein_tmp
        out_data$Peptide[i2] = peptide_tmp
      }
      all_intensity_ratios_DD = rbind(all_intensity_ratios_DD, out_data)
    }
  }
  rm(my_data_cleaned, data_tmp_1, data_tmp_2, i_peptide, out_data)
  
  
  # Assign CVs to AI ratios
  all_intensity_ratios = get_Koina_all_ratios(my_data, grouping_data, set_Koina_model)
  all_intensity_ratios = cbind(all_intensity_ratios, rep(NA, nrow(all_intensity_ratios)))
  colnames(all_intensity_ratios)[length(colnames(all_intensity_ratios))] = "CV"
  all_intensity_ratios_final = data.frame(matrix(NA, nrow=0, ncol=5))
  colnames(all_intensity_ratios_final) = colnames(all_intensity_ratios)
  for (i_ratio_peptide in 1:length(unique(all_intensity_ratios$Peptide_name))) {
    peptide_tmp = unique(all_intensity_ratios$Peptide_name)[i_ratio_peptide]
    all_intensity_ratios_tmp = all_intensity_ratios[all_intensity_ratios$Peptide_name == peptide_tmp, ]
    
    all_intensity_ratios_DD_tmp = all_intensity_ratios_DD[all_intensity_ratios_DD$Peptide == peptide_tmp, 1:4]
    
    for (i_ratio in 1:nrow(all_intensity_ratios_tmp)) {
      ratio_name_tmp = all_intensity_ratios_tmp$Ratio_name[i_ratio]
      
      ratio_name_tmp_part_1 = gsub(" ", "", strsplit(strsplit(ratio_name_tmp, "\\/")[[1]][1], "\\|")[[1]][2])
      ratio_name_tmp_part_2 = gsub(" ", "", strsplit(strsplit(ratio_name_tmp, "\\/")[[1]][2], "\\|")[[1]][2])
      ratio_name_tmp = paste(ratio_name_tmp_part_1, " / ", ratio_name_tmp_part_2, sep="")
      
      if (length(which(all_intensity_ratios_DD_tmp$Ratio_name == ratio_name_tmp)) != 0) {
        CV_tmp = all_intensity_ratios_DD_tmp$CV[which(all_intensity_ratios_DD_tmp$Ratio_name == ratio_name_tmp)]
        all_intensity_ratios_tmp$CV[i_ratio] = CV_tmp
      }
    }
    all_intensity_ratios_final = rbind(all_intensity_ratios_final, all_intensity_ratios_tmp)
  }
  all_intensity_ratios = all_intensity_ratios_final
  
  
  # Run MDI
  my_data_core_before = my_data[, which(colnames(my_data) %in% grouping_data$Replicate)]
  my_data_new_additional = my_data[, -(which(colnames(my_data) %in% grouping_data$Replicate))]
  my_data_new_core = my_data[, which(colnames(my_data) %in% grouping_data$Replicate)]
  my_data_new_core[,] = NA
  
  for (i_peptide in 1:length(unique(my_data$Peptide))) {
    peptide_tmp = unique(my_data$Peptide)[i_peptide]
    my_data_tmp = my_data[my_data$Peptide == peptide_tmp, ]
    my_data_additional_tmp = my_data_additional[my_data_additional$Peptide == peptide_tmp, ]
    data_imputation_rows = which((my_data_new_additional$Peptide == peptide_tmp) == TRUE)
    fragments_unique = paste("prec+", my_data_additional_tmp$`Precursor Charge`, " | ", my_data_additional_tmp$`Fragment Ion`, 
                             "+", my_data_additional_tmp$`Product Charge`, sep="")
    my_data_tmp_additional = my_data_tmp[, -(which(colnames(my_data_tmp) %in% grouping_data$Replicate))]
    my_data_tmp = my_data_tmp[, which(colnames(my_data_tmp) %in% grouping_data$Replicate)]
    my_data_tmp$Fragment = paste("prec+", my_data_additional_tmp$`Precursor Charge`, " | ", my_data_additional_tmp$`Fragment Ion`, 
                                 "+", my_data_additional_tmp$`Product Charge`, sep="")
    fragments_unique = my_data_tmp$Fragment
    my_data_tmp = cbind(as.factor(fragments_unique), my_data_tmp)
    colnames(my_data_tmp)[1] = "Fragment"
    my_data_tmp = my_data_tmp[,-ncol(my_data_tmp)]
    
    if (nrow(my_data_tmp) <= 1) {
      my_data_new_core = my_data_new_core[-(which((my_data_new_core$Peptide == peptide_tmp) == TRUE)), ]
      my_data_new_additional = my_data_new_additional[-(which((my_data_new_additional$Peptide == peptide_tmp) == TRUE)), ]
    } else {
      
      for (i_sample in 2:ncol(my_data_tmp)) {
        sample_tmp = my_data_tmp[,i_sample]
        fragments_tmp = as.character(my_data_tmp[,1])
        
        if (sum(is.na(my_data_tmp[,i_sample])) == 0) {
          my_data_new_core[data_imputation_rows, i_sample-1] = sample_tmp
        } else {
          fragments_combinations_tmp = combn(fragments_tmp,2)
          
          for (i_row in 1:length(sample_tmp)) {
            sample_NA_tmp = sample_tmp[i_row]
            fragments_NA_tmp = fragments_tmp[i_row]
            data_imputation_rows_tmp = data_imputation_rows[1] + i_row - 1
            
            if (!is.na(sample_NA_tmp)) {
              my_data_new_core[data_imputation_rows_tmp, i_sample-1] = sample_NA_tmp
            } else {
              
              comb_inx = c((which(fragments_NA_tmp == fragments_combinations_tmp[1,])),
                           (which(fragments_NA_tmp == fragments_combinations_tmp[2,])))
              comb_inx = fragments_combinations_tmp[,comb_inx]
              if (!is.null(ncol(comb_inx))) {
                comb_inx = paste(comb_inx[1,], " / ", comb_inx[2,], sep="")
              } else {
                comb_inx = paste(comb_inx[1], " / ", comb_inx[2], sep="")
              }
              
              tmp_intensity_ratios = all_intensity_ratios[which(peptide_tmp == all_intensity_ratios$Peptide),]
              
              if (nrow(tmp_intensity_ratios) == 0) {
                my_data_new_core[data_imputation_rows_tmp, i_sample-1] = sample_NA_tmp
              } else {
                
                comb_inx_matched = c()
                for (i_ratios_tmp in 1:nrow(tmp_intensity_ratios)) {
                  comb_inx_matched = c(comb_inx_matched, which(comb_inx[i_ratios_tmp] == tmp_intensity_ratios$Ratio_name))
                }
                
                if (length(comb_inx_matched) == 0) {
                  my_data_new_core[data_imputation_rows_tmp, i_sample-1] = sample_NA_tmp
                } else {
                  tmp_intensity_ratios = tmp_intensity_ratios[comb_inx_matched, ]
                  
                  if (choose_ref == "heaviest") {
                    ratios_weight_1 = sapply(strsplit(sapply(strsplit(sapply(strsplit(tmp_intensity_ratios$Ratio_name, " | "), "[", 3), "\\+"), "[", 1), "[a-z]+"), "[", 2)
                    ratios_weight_2 = sapply(strsplit(sapply(strsplit(sapply(strsplit(tmp_intensity_ratios$Ratio_name, " | "), "[", 7), "\\+"), "[", 1), "[a-z]+"), "[", 2)
                    
                    # Take ratios of heaviest fragments
                    ratios_weight = as.numeric(ratios_weight_1) * as.numeric(ratios_weight_2)
                    frag_ref_order = order(ratios_weight, decreasing = TRUE)
                    
                  } else if (choose_ref == "CVs") {
                  # Order ratios by CVs
                    frag_ref_order = order(tmp_intensity_ratios$CV, decreasing = FALSE)
                  }
                  
                  frag_ref_flag = 0
                  i_frag_ref = 0
                  while (frag_ref_flag == 0) {
                    i_frag_ref = i_frag_ref + 1
                    
                    if (is.na(frag_ref_order[i_frag_ref])) {
                      my_data_new_core[data_imputation_rows_tmp, i_sample-1] = NA
                      writeLines("Some fragments could not be matched using AI")
                      break
                    }
                    
                    frag_ref = tmp_intensity_ratios$Ratio_name[frag_ref_order[i_frag_ref]]
                    matched_fragments = strsplit(frag_ref, split = " / ")[[1]]
                    if (which(fragments_NA_tmp == matched_fragments) == 1) {
                      input_value = my_data_tmp[which(matched_fragments[2] == my_data_tmp$Fragment), i_sample] * tmp_intensity_ratios$Ratio[which(tmp_intensity_ratios$Ratio_name == frag_ref)]
                      
                      if (!is.na(input_value)) {
                        my_data_new_core[data_imputation_rows_tmp, i_sample-1] = input_value
                        frag_ref_flag = 1
                      }
                    } else if (which(fragments_NA_tmp == matched_fragments) == 2) {
                      input_value = my_data_tmp[which(matched_fragments[1] == my_data_tmp$Fragment), i_sample] / tmp_intensity_ratios$Ratio[which(tmp_intensity_ratios$Ratio_name == frag_ref)]
                      
                      if (!is.na(input_value)) {
                        my_data_new_core[data_imputation_rows_tmp, i_sample-1] = input_value
                        frag_ref_flag = 1
                      }
                    } else {
                      print("Something went wrong with MDI AI")
                    }
                  }
                }
              }
            }
            
            if (is.na(my_data_new_core[data_imputation_rows_tmp, i_sample-1])) {
              print(paste(i_peptide, i_sample, i_row))
            }
            
          }
        }
      }
    }
  }
  
  my_data_new = cbind(my_data_new_additional, my_data_new_core)
  
  changed_values = is.na(my_data_core_before) * my_data_new_core
  changed_values[changed_values == 0] = NA
  
  if (print_MD_procentage == TRUE) {
    miss_data_proc = (sum(is.na(my_data_new)) / (nrow(my_data_new) * length(which(colnames(my_data_new) %in% grouping_data$Replicate)))) * 100
    print(paste("Missing data (step 2 Fragment ratios): ",  round(miss_data_proc,3), "%", sep=""))
  }
  
  return_list = list("my_data_new" = my_data_new, "changed_values" = changed_values)
  return(return_list)
  
}

get_Koina_all_ratios <- function(my_data, grouping_data, 
                                 set_Koina_model = "ms2pip_timsTOF2024") {
  
  library(koinar, quietly=T, warn.conflicts=F)
  
  my_data_additional = my_data[, -which(colnames(my_data) %in% grouping_data$Replicate)]
  
  # Get all possible fragment combinations
  if (set_Koina_model == "ms2pip_timsTOF2024") {
    my_data_export = my_data[,c(1,2,5,8,9)]
    my_data_export_wo_fragments = my_data_export[,-4]
    my_data_export_wo_fragments = unique(my_data_export_wo_fragments)
    
    model = koinar::Koina$new(model_name = set_Koina_model,
                              server_url = "koina.wilhelmlab.org:443",
                              ssl = TRUE)
    
    input_data = data.frame(list(list("peptide_sequences"= my_data_export_wo_fragments$Peptide),
                                 list("precursor_charges"= my_data_export_wo_fragments$`Precursor Charge`),
                                 list("collision_energies"= rep(31, length(my_data_export_wo_fragments$`Precursor Charge`)))))
    writeLines("Contacting Koina... Running Predictions...")
    predictions = model$predict(input_data)
    writeLines("All done!")
    writeLines("")
    
  } else if (set_Koina_model == "Prosit_2020_intensity_HCD") {
    my_data_export = my_data[,c(1,2,5,8,9)]
    my_data_export_wo_fragments = my_data_export[,-4]
    my_data_export_wo_fragments$`Collision Energy` = round(my_data_export_wo_fragments$`Collision Energy`, digits = 0)
    my_data_export_wo_fragments = unique(my_data_export_wo_fragments)

    model = koinar::Koina$new(model_name = "Prosit_2020_intensity_HCD",
                              server_url = "koina.wilhelmlab.org:443",
                              ssl = TRUE)

    input_data = data.frame(list(list("peptide_sequences"= my_data_export_wo_fragments$Peptide),
                                 list("precursor_charges"= my_data_export_wo_fragments$`Precursor Charge`),
                                 list("collision_energies"= my_data_export_wo_fragments$`Collision Energy`)))
    writeLines("Contacting Koina... Running Predictions...")
    predictions = model$predict(input_data)
    writeLines("All done!")
    writeLines("")

  } else {
    print("No AI model has been selected! Aborting...")
    break
  }
  
  my_data_core_before = my_data[, which(colnames(my_data) %in% grouping_data$Replicate)]
  
  my_data_cleaned = my_data
  my_data_cleaned$Fragment = paste("prec+", my_data_cleaned$`Precursor Charge`, " | ", my_data_cleaned$`Fragment Ion`, 
                                   "+", my_data_cleaned$`Product Charge`, sep="")
  all_intensity_ratios = data.frame(matrix(NaN, nrow = 0, ncol = 4))
  colnames(all_intensity_ratios) = c("Protein_name","Peptide_name","Ratio_name","Ratio")
  
  for (i_peptide in 1:length(unique(my_data_cleaned$Peptide))) {
    peptide_tmp = unique(my_data_cleaned$Peptide)[i_peptide]
    my_data_tmp = my_data_cleaned[my_data_cleaned$Peptide == peptide_tmp, ]
    fragments_unique = my_data_tmp$Fragment
    my_data_tmp_additional = my_data_tmp[, -(which(colnames(my_data_tmp) %in% grouping_data$Replicate))]
    protein_tmp = my_data_tmp_additional$`Protein Name`[1]
    my_data_tmp = my_data_tmp[, which(colnames(my_data_tmp) %in% grouping_data$Replicate)]
    my_data_tmp = cbind(as.factor(fragments_unique), my_data_tmp)
    colnames(my_data_tmp)[1] = "Fragment"
    
    my_data_predicted_tmp = predictions[which(predictions$peptide_sequences == peptide_tmp),]
    my_data_predicted_tmp$annotation = paste("prec+", my_data_predicted_tmp$precursor_charges, " | ", 
                                             my_data_predicted_tmp$annotation, sep="")
    my_data_predicted_tmp = my_data_predicted_tmp[,-2]
    my_data_predicted_tmp = my_data_predicted_tmp[which(my_data_predicted_tmp$annotation %in% my_data_tmp[,1]),]
    
    if (nrow(my_data_tmp) > 1) {
      all_combinations = combn(as.character(my_data_tmp$Fragment),2)
      
      out_data = data.frame(matrix(NaN, nrow = ncol(all_combinations), ncol = 4))
      colnames(out_data) = c("Protein_name","Peptide_name","Ratio_name","Ratio")
      for (i2 in 1:ncol(all_combinations)) {
        frag_tmp_1 = all_combinations[1,i2]
        frag_tmp_2 = all_combinations[2,i2]
        
        data_tmp_1 = my_data_predicted_tmp$intensities[which(my_data_predicted_tmp$annotation == frag_tmp_1)]
        data_tmp_2 = my_data_predicted_tmp$intensities[which(my_data_predicted_tmp$annotation == frag_tmp_2)]
        
        if (!(sum(nchar(as.character(data_tmp_1))) == 0) & !(sum(nchar(as.character(data_tmp_2))) == 0)) {
          out_data$Ratio[i2] = (as.numeric(data_tmp_1) / as.numeric(data_tmp_2))
        } else {
          out_data$Ratio[i2] = NA
        }
        out_data$Ratio_name[i2] = paste(frag_tmp_1 , " / ", frag_tmp_2, sep="")
        out_data$Protein_name[i2] = protein_tmp
        out_data$Peptide_name[i2] = peptide_tmp
      }
      all_intensity_ratios = rbind(all_intensity_ratios, out_data)
    }
  }
  
  all_intensity_ratios = na.omit(all_intensity_ratios)
  
  return(all_intensity_ratios)
}

data_combination_summing <- function(my_data, grouping_data, sum_level = "Protein") {
  my_data_additional = my_data[, -which(colnames(my_data) %in% grouping_data$Replicate)]
  
  out_data_combination = data.frame(matrix(NaN, nrow = 1, ncol = (ncol(my_data)-ncol(my_data_additional)+1)))
  colnames(out_data_combination) = c("Protein", colnames(my_data)[(ncol(my_data_additional)+1):ncol(my_data)])
  
  if (sum_level == "Protein") {
    for (i_protein in 1:length(unique(my_data$`Protein`))) {
      protein_tmp = (unique(my_data$`Protein Name`))[i_protein]
      my_data_tmp = my_data[my_data$`Protein Name` == protein_tmp, ]
      my_data_tmp = my_data_tmp[,(ncol(my_data_additional)+1):ncol(my_data_tmp)]
      out_data_combination[i_protein, 1] = protein_tmp
      out_data_combination[i_protein, 2:ncol(out_data_combination)] = colSums(as.data.frame(sapply(my_data_tmp, as.numeric)))
    }
  } else if (sum_level == "Peptide") {
    colnames(out_data_combination)[1] = "Peptide"
    for (i_peptide in 1:length(unique(my_data$`Peptide`))) {
      peptide_tmp = (unique(my_data$`Peptide`))[i_peptide]
      my_data_tmp = my_data[my_data$`Peptide` == peptide_tmp, ]
      my_data_tmp = my_data_tmp[,(ncol(my_data_additional)+1):ncol(my_data_tmp)]
      out_data_combination[i_peptide, 1] = peptide_tmp
      out_data_combination[i_peptide, 2:ncol(out_data_combination)] = colSums(as.data.frame(sapply(my_data_tmp, as.numeric)))
    }
  }
  return(out_data_combination)
}

data_combination_choose_best <- function(my_data, grouping_data, 
                                         sum_level = "Protein", max_type = "global") {
  
  library(DescTools, quietly=T, warn.conflicts=F)
  
  my_data_additional = my_data[, -which(colnames(my_data) %in% grouping_data$Replicate)]
  
  out_data_combination = data.frame(matrix(NaN, nrow = 1, ncol = (ncol(my_data)-ncol(my_data_additional)+1)))
  colnames(out_data_combination) = c("Protein", colnames(my_data)[(ncol(my_data_additional)+1):ncol(my_data)])
  
  if (sum_level == "Protein") {
    for (i_protein in 1:length(unique(my_data$`Protein`))) {
      for (i_protein in 1:length(unique(my_data$`Protein Name`))) {
        protein_tmp = (unique(my_data$`Protein Name`))[i_protein]
        my_data_tmp = my_data[my_data$`Protein Name` == protein_tmp, ]
        my_data_tmp = my_data_tmp[,(ncol(my_data_additional)+1):ncol(my_data_tmp)]
        out_data_combination[i_protein, 1] = protein_tmp
        
        if (max_type == "global") {
          max_frag_int = which.max(as.vector(rowMedians(as.matrix(my_data_tmp))))
          out_data_combination[i_protein, 2:ncol(out_data_combination)] = as.vector(sapply(my_data_tmp[max_frag_int,], as.numeric))
        } else if (max_type == "local") {
          
          which_max_vector = c()
          for (i_col in 1:ncol(my_data_tmp)) {
            which_max_vector = c(which_max_vector, which.max(my_data_tmp[,i_col]))
          }
          
          if (length(unique(which_max_vector)) == 1) {
            out_data_combination[i_protein, 2:ncol(out_data_combination)] = as.vector(sapply(my_data_tmp[which_max_vector[1],], as.numeric))
          } else {
            
            which_most_often = as.numeric(sort(table(which_max_vector),decreasing=TRUE)[1])
            
            for (i_vec in 1:length(which_max_vector)) {
              if (which_max_vector[i_vec] == which_most_often) {
                out_data_combination[i_protein, i_vec+1] = as.vector(sapply(my_data_tmp[which_max_vector[i_vec], i_vec], as.numeric))
              } else {
                ratio_tmp = mean(as.numeric(my_data_tmp[which_max_vector[i_vec],] / my_data_tmp[which_most_often,]))
                out_data_combination[i_protein, i_vec+1] = as.numeric(my_data_tmp[which_max_vector[i_vec],i_vec]) / ratio_tmp
              }
            }
          }
        }
      }
    }
  } else if (sum_level == "Peptide") {
    colnames(out_data_combination)[1] = "Peptide"
    for (i_peptide in 1:length(unique(my_data$`Peptide`))) {
      peptide_tmp = (unique(my_data$`Peptide`))[i_peptide]
      my_data_tmp = my_data[my_data$`Peptide` == peptide_tmp, ]
      my_data_tmp = my_data_tmp[,(ncol(my_data_additional)+1):ncol(my_data_tmp)]
      out_data_combination[i_peptide, 1] = peptide_tmp
      
      if (max_type == "global") {
        max_frag_int = which.max(as.vector(rowMedians(as.matrix(my_data_tmp))))
        out_data_combination[i_peptide, 2:ncol(out_data_combination)] = as.vector(sapply(my_data_tmp[max_frag_int,], as.numeric))
      } else if (max_type == "local") {
        
        which_max_vector = c()
        for (i_col in 1:ncol(my_data_tmp)) {
          which_max_vector = c(which_max_vector, which.max(my_data_tmp[,i_col]))
        }
        
        if (length(unique(which_max_vector)) == 1) {
          out_data_combination[i_peptide, 2:ncol(out_data_combination)] = as.vector(sapply(my_data_tmp[which_max_vector[1],], as.numeric))
        } else {
          
          which_most_often = as.numeric(DescTools::Mode(which_max_vector))
          
          if (length(which_most_often) > 1) {
            which_most_often = which_most_often[1]
          }
          
          
          for (i_vec in 1:length(which_max_vector)) {
            if (which_max_vector[i_vec] == which_most_often) {
              out_data_combination[i_peptide, i_vec+1] = as.vector(sapply(my_data_tmp[which_max_vector[i_vec], i_vec], as.numeric))
            } else {
              ratio_tmp = mean(as.numeric(my_data_tmp[which_max_vector[i_vec],] / my_data_tmp[which_most_often,]))
              out_data_combination[i_peptide, i_vec+1] = as.numeric(my_data_tmp[which_max_vector[i_vec],i_vec]) / ratio_tmp
            }
          }
        }
      }
    }
  }
  return(out_data_combination)
}

data_rescaling_AI_without_change <- function(my_data, grouping_data, set_Koina_model = "ms2pip_timsTOF2024",
                                             choose_ref = "highest", consolidation_method = "summing", 
                                             print_stacked_peptides = FALSE, print_fragment_variation = FALSE,
                                             print_fragment_ratio_variation = FALSE, show_double_charge_warn = TRUE) {
  all_intensity_ratios = get_Koina_all_ratios(my_data, grouping_data, set_Koina_model)
  
  # Change fragment ratios
  my_data_out_tmp = my_data
  my_data_out_tmp[,(ncol(my_data_additional)+1):ncol(my_data_out_tmp)] = NA
  
  
  if (length(which(my_data_out_tmp$`Product Charge` == 2)) != 0) {
    my_data_out_tmp = my_data_out_tmp[-which(my_data_out_tmp$`Product Charge` == 2),]
  }

  for (i_peptide in 1:length(unique(my_data$Peptide))) {
    peptide_tmp = unique(my_data$Peptide)[i_peptide]
    my_data_tmp = my_data[my_data$Peptide == peptide_tmp, ]
    protein_tmp = my_data_tmp[1,1]
    my_data_additional_tmp = my_data_additional[my_data_additional$Peptide == peptide_tmp, ]
    my_data_tmp = my_data_tmp[,(ncol(my_data_additional)+1):ncol(my_data_tmp)]
    my_data_tmp_fragments = paste("prec+", my_data_additional_tmp$`Precursor Charge`, " | ", my_data_additional_tmp$`Fragment Ion`, 
                                  "+", my_data_additional_tmp$`Product Charge`, sep="")
    my_data_tmp = cbind(my_data_additional_tmp$`Protein Name`, my_data_additional_tmp$Peptide, my_data_tmp_fragments, my_data_tmp)
    colnames(my_data_tmp)[1:3] = c("Protein", "Peptide", "Fragment")
    all_intensity_ratios_tmp = all_intensity_ratios[all_intensity_ratios$Peptide == peptide_tmp, ]
    
    # Print fragment variation
    if (print_fragment_variation == TRUE) {
      print_fragment_variation_plot(my_data_tmp)
    }
    
    # Doubly charged fragments cannot be predicted - removing
    doubly_removed = c()
    for (i_row in 1:nrow(my_data_tmp)) {
      doubly_removed = c(doubly_removed, strsplit(strsplit(my_data_tmp$Fragment[i_row], " | ")[[1]][3], "+", fixed=TRUE)[[1]][2] != 1)
    }
    if (sum(doubly_removed) > 0) {
      my_data_tmp = my_data_tmp[-which(doubly_removed == TRUE),]
      
      
      if (show_double_charge_warn == TRUE) {
        writeLines(paste("Warning: Peptide ", peptide_tmp, 
                         " contains fragment that is doubly charged! Removing from analysis...", sep=""))
      }
    }
    
    # Print fragment ratio variation
    if (print_fragment_ratio_variation == TRUE) {
      print_fragment_ratio_variation_plot(my_data_tmp)
    }
    
    # Find a reference peptide - heaviest or higest peak area
    if (choose_ref == "heaviest") {
      frag_ref_number = sapply(strsplit(sapply(strsplit(my_data_tmp$Fragment, " | "), "[", 3), "\\+"), "[", 1)
      frag_ref_number = as.numeric(which.max(sapply(strsplit(frag_ref_number, "[a-z]+"), "[", 2)))
    } else if (choose_ref == "highest") {
      frag_ref_number = data.frame(sapply(my_data_tmp[,4:ncol(my_data_tmp)], function(x) as.numeric(as.character(x))))
      frag_ref_number = sort(rowMedians(as.matrix(frag_ref_number)), decreasing = TRUE)
      test_flag = 0
      test_inx = 1
      while (test_flag == 0) {
        frag_ref_name = my_data_tmp$Fragment[which(frag_ref_number == frag_ref_number[test_inx])]
        test_frag = which(paste(frag_ref_name, " / ", my_data_tmp$Fragment[-which(my_data_tmp$Fragment == frag_ref_name)], 
                                sep="") %in% all_intensity_ratios_tmp$Ratio_name)
        if (length(paste(frag_ref_name, " / ", my_data_tmp$Fragment[-which(my_data_tmp$Fragment == frag_ref_name)], sep="")) == length(test_frag)) {
          test_flag = 1
        } else {
          test_inx = test_inx + 1
        }
      }
      frag_ref_number = which(frag_ref_number == frag_ref_number[test_inx])
    }
    
    if (length(frag_ref_number) > 1) {
      frag_ref_number = min(frag_ref_number)
    }
    
    
    for (i_col in 4:ncol(my_data_tmp)) {
      for (i_row in (1:nrow(my_data_tmp))[! (1:nrow(my_data_tmp)) %in% frag_ref_number]) {
        if (!(strsplit(strsplit(my_data_tmp$Fragment[frag_ref_number], " | ")[[1]][3], "+", fixed=TRUE)[[1]][2] == 2)) {
          frag_ref_name = my_data_tmp$Fragment[frag_ref_number]
          
          if (!(strsplit(strsplit(my_data_tmp$Fragment[i_row], " | ")[[1]][3], "+", fixed=TRUE)[[1]][2] == 2)) {
            current_frag_name = my_data_tmp$Fragment[i_row]
            
            current_ratio_name = paste(frag_ref_name, " / ", current_frag_name, sep="")
            current_ratio = all_intensity_ratios_tmp$Ratio[all_intensity_ratios_tmp$Ratio_name == current_ratio_name]
            
            if (length(current_ratio) != 0) {
              my_data_tmp[i_row, i_col] = as.numeric(my_data_tmp[my_data_tmp$Fragment == frag_ref_name, ][i_col]) / current_ratio
            } else {
              current_ratio_name = paste(current_frag_name, " / ", frag_ref_name, sep="")
              current_ratio = all_intensity_ratios_tmp$Ratio[all_intensity_ratios_tmp$Ratio_name == current_ratio_name]
              
              my_data_tmp[i_row, i_col] = as.numeric(my_data_tmp[my_data_tmp$Fragment == frag_ref_name, ][i_col]) * current_ratio
            }
          }
        } else {
          writeLines("Warning: Reference fragment is doubly charged!")
        }
      }
    }
    
    my_data_out_tmp[rownames((my_data_tmp)),(ncol(my_data_additional)+1):ncol(my_data_out_tmp)] = my_data_tmp[,-c(1:3)]
    
    # Print stacked plots
    if (print_stacked_peptides == TRUE) {
      print_stacked_plots_peptides(my_data_tmp)
    }
  }
  
  # Consolidate data either by simple sum or using highest (by median) peaks
  my_data_out_tmp = my_data_out_tmp[complete.cases(my_data_out_tmp[,(ncol(my_data_additional)+1):ncol(my_data_out_tmp)]), ]
  if (consolidation_method == "median") {
    out_data_combination = data.frame(matrix(NaN, nrow = 0, ncol = (ncol(my_data_out_tmp)-ncol(my_data_additional)+1)))
    colnames(out_data_combination)[1] = c("Protein Name")
    colnames(out_data_combination)[2:(ncol(out_data_combination))] = colnames(my_data_out_tmp)[(max(
      which(colnames(my_data_out_tmp) %in% colnames(my_data_additional))) + 1) : ncol(my_data_out_tmp)]
    for (i_protein in 1:length(unique(my_data_out_tmp$`Protein Name`))) {
      protein_tmp = unique(my_data_out_tmp$`Protein Name`)[i_protein]
      my_data_tmp = my_data_out_tmp[my_data_out_tmp$`Protein Name` == protein_tmp, ]
      my_data_additional_tmp = my_data_additional[my_data_additional$`Protein Name` == protein_tmp, ]
      my_data_tmp = my_data_tmp[,(ncol(my_data_additional)+1):ncol(my_data_tmp)]
      my_data_tmp = data.frame(sapply(my_data_tmp, function(x) as.numeric(as.character(x))))
      median_values = rowMedians(as.matrix(my_data_tmp))
      out_data_combination[nrow(out_data_combination)+1, 1] = protein_tmp
      out_data_combination[nrow(out_data_combination), 2:ncol(out_data_combination)] = my_data_tmp[which.max(median_values), ]
    }
  } else if (consolidation_method == "summing") {
    out_data_combination = data_combination_summing(my_data_out_tmp, grouping_data)
  } else {
    out_data_combination = my_data_out_tmp
  }
  
  return(out_data_combination)
  
}

data_rescaling_AI_to_normalized <- function(my_data, grouping_data, set_Koina_model, 
                                            choose_ref, consolidation_method, 
                                            print_stacked_peptides, print_fragment_variation,
                                            print_fragment_ratio_variation) {
  
  all_intensity_ratios = get_Koina_all_ratios(my_data, grouping_data, set_Koina_model)
  
  # Change fragment ratios
  my_data_out_tmp = my_data
  my_data_out_tmp[,(ncol(my_data_additional)+1):ncol(my_data_out_tmp)] = NA
  for (i_peptide in 1:length(unique(my_data$Peptide))) {
    peptide_tmp = unique(my_data$Peptide)[i_peptide]
    my_data_tmp = my_data[my_data$Peptide == peptide_tmp, ]
    protein_tmp = my_data_tmp[1,1]
    my_data_additional_tmp = my_data_additional[my_data_additional$Peptide == peptide_tmp, ]
    my_data_tmp = my_data_tmp[,(ncol(my_data_additional)+1):ncol(my_data_tmp)]
    my_data_tmp_fragments = paste("prec+", my_data_additional_tmp$`Precursor Charge`, " | ", my_data_additional_tmp$`Fragment Ion`, 
                                  "+", my_data_additional_tmp$`Product Charge`, sep="")
    my_data_tmp = cbind(my_data_additional_tmp$`Protein Name`, my_data_additional_tmp$Peptide, my_data_tmp_fragments, my_data_tmp)
    colnames(my_data_tmp)[1:3] = c("Protein", "Peptide", "Fragment")
    all_intensity_ratios_tmp = all_intensity_ratios[all_intensity_ratios$Peptide == peptide_tmp, ]
    
    # Print fragment variation
    if (print_fragment_variation == TRUE) {
      print_fragment_variation_plot(my_data_tmp)
    }
    
    # Doubly charged fragments cannot be predicted - removing
    doubly_removed = c()
    for (i_row in 1:nrow(my_data_tmp)) {
      doubly_removed = c(doubly_removed, strsplit(strsplit(my_data_tmp$Fragment[i_row], " | ")[[1]][3], "+", fixed=TRUE)[[1]][2] != 1)
    }
    if (sum(doubly_removed) > 0) {
      my_data_tmp = my_data_tmp[-which(doubly_removed == TRUE),]
      writeLines(paste("Warning: Peptide ", peptide_tmp, 
                       " contains fragment that is doubly charged! Removing from analysis...", sep=""))
    }
    
    # Print fragment ratio variation
    if (print_fragment_ratio_variation == TRUE) {
      print_fragment_ratio_variation_plot(my_data_tmp)
    }
    
    # Find a reference peptide - heaviest or highest peak area
    if (choose_ref == "heaviest") {
      frag_ref_number = sapply(strsplit(sapply(strsplit(my_data_tmp$Fragment, " | "), "[", 3), "\\+"), "[", 1)
      frag_ref_number = as.numeric(which.max(sapply(strsplit(frag_ref_number, "[a-z]+"), "[", 2)))
    } else if (choose_ref == "highest") {
      frag_ref_number = data.frame(sapply(my_data_tmp[,4:ncol(my_data_tmp)], function(x) as.numeric(as.character(x))))
      frag_ref_number = which.max(rowMedians(as.matrix(frag_ref_number)))
    }
    frag_ref_name = my_data_tmp$Fragment[frag_ref_number]
    
    # Re-scale data to the reference fragment
    for (i_col in 4:ncol(my_data_tmp)) {
      for (i_row in 1:nrow(my_data_tmp)) {
        current_frag_name = my_data_tmp$Fragment[i_row]
        if (current_frag_name != frag_ref_name) {
          current_ratio_name = paste(frag_ref_name, " / ", current_frag_name, sep="")
          current_ratio = all_intensity_ratios_tmp$Ratio[all_intensity_ratios_tmp$Ratio_name == current_ratio_name]
          if (length(current_ratio != 0)) {
            my_data_tmp[i_row, i_col] = as.numeric(my_data_tmp[i_row, i_col]) * current_ratio
          } else {
            current_ratio_name = paste(current_frag_name, " / ", frag_ref_name, sep="")
            current_ratio = all_intensity_ratios_tmp$Ratio[all_intensity_ratios_tmp$Ratio_name == current_ratio_name]
            if (length(current_ratio != 0)) {
              my_data_tmp[i_row, i_col] = as.numeric(my_data_tmp[i_row, i_col]) / current_ratio
            } else {
              writeLines("Something went wrong... Can't set ratios.")
              writeLines("")
            }
          }
        }
      }
    }
    
    my_data_out_tmp[rownames((my_data_tmp)),(ncol(my_data_additional)+1):ncol(my_data_out_tmp)] = my_data_tmp[,-c(1:3)]
    
    # Print stacked plots
    if (print_stacked_peptides == TRUE) {
      if (peptide_tmp == "GNPTVEVELTTEK" | peptide_tmp == "TGHIAADGSVYNR") {
        print_stacked_plots_peptides(my_data_tmp)
      }
    }
  }
  writeLines("")
  
  # Consolidate data either by simple sum or using highest (by median) peaks
  my_data_out_tmp = my_data_out_tmp[complete.cases(my_data_out_tmp[,(ncol(my_data_additional)+1):ncol(my_data_out_tmp)]), ]
  if (consolidation_method == "median") {
    out_data_combination = data.frame(matrix(NaN, nrow = 0, ncol = (ncol(my_data_out_tmp)-ncol(my_data_additional)+1)))
    colnames(out_data_combination)[1] = c("Protein Name")
    colnames(out_data_combination)[2:(ncol(out_data_combination))] = colnames(my_data_out_tmp)[(max(
      which(colnames(my_data_out_tmp) %in% colnames(my_data_additional))) + 1) : ncol(my_data_out_tmp)]
    for (i_protein in 1:length(unique(my_data_out_tmp$`Protein Name`))) {
      protein_tmp = unique(my_data_out_tmp$`Protein Name`)[i_protein]
      my_data_tmp = my_data_out_tmp[my_data_out_tmp$`Protein Name` == protein_tmp, ]
      my_data_additional_tmp = my_data_additional[my_data_additional$`Protein Name` == protein_tmp, ]
      my_data_tmp = my_data_tmp[,(ncol(my_data_additional)+1):ncol(my_data_tmp)]
      my_data_tmp = data.frame(sapply(my_data_tmp, function(x) as.numeric(as.character(x))))
      median_values = rowMedians(as.matrix(my_data_tmp))
      out_data_combination[nrow(out_data_combination)+1, 1] = protein_tmp
      out_data_combination[nrow(out_data_combination), 2:ncol(out_data_combination)] = my_data_tmp[which.max(median_values), ]
    }
  } else if (consolidation_method == "summing") {
    out_data_combination = data_combination_summing(my_data_out_tmp, grouping_data)
  } else {
    writeLines("No data consolidation method chosen - using simple sum.")
    writeLines("")
    out_data_combination = data_combination_summing(my_data_out_tmp, grouping_data)
  }
  
  return(out_data_combination)
  
}

data_combination_p_int_Browns <- function(my_data, grouping_data, int_level = "Peptide") {
  library(EmpiricalBrownsMethod)
  
  if (int_level == "Protein") {
    
    writeLines("P-value integration non-reliable! Fix this later")
    
    my_data_out = data.frame(matrix(NaN, nrow = 0, ncol = 3))
    colnames(my_data_out) = c("Protein Name","p_value", "Signif")
    for (i_protein in (1:length(unique(my_data$`Protein Name`)))) {
      protein_tmp = unique(my_data$`Protein Name`)[i_protein]
      my_data_tmp = my_data[which(my_data$`Protein Name` == protein_tmp),]
      my_data_tmp_additional = my_data[, -which(colnames(my_data) %in% grouping_data$Replicate)]
      
      my_data_tmp_IDs = paste(my_data_tmp$`Protein Name`, " | ", my_data_tmp$Peptide,
                              " | ","prec+", my_data_tmp$`Precursor Charge`,
                              " | ", my_data_tmp$`Fragment Ion`, "+", my_data_tmp$`Product Charge`, sep="")
      
      my_data_tmp = my_data_tmp[,-which(colnames(my_data_tmp) %in% colnames(my_data_tmp_additional))]
      rownames(my_data_tmp) = my_data_tmp_IDs
      my_data_tmp = rbind(grouping_data$Condition[26:35], my_data_tmp)
      rownames(my_data_tmp)[1] = c("Group")
      my_data_tmp = as.data.frame(t(my_data_tmp))
      my_data_tmp[, (2:ncol(my_data_tmp))] = apply(my_data_tmp[, (2:ncol(my_data_tmp))], 2, function(x) as.numeric(as.character(x)))
      my_data_tmp = my_data_tmp %>% janitor::clean_names()
      my_data_tmp = sort_by(my_data_tmp, as.numeric(my_data_tmp$group))
      my_data_tmp[,1] = factor(my_data_tmp[,1], levels = as.numeric(unique(my_data_tmp[,1])))
      colnames(my_data_tmp)[1] = c("Group")
      
      p_values = c()
      for (i_fragment in 1:ncol(my_data_tmp[,-1])) {
        data_test = my_data_tmp[, c(1, i_fragment+1)]
        
        shapiro_p_value = shapiro.test(data_test[which(data_test$Group == levels(data_test$Group)[1]), 2])$p.value
        shapiro_p_value = c(shapiro_p_value, shapiro.test(data_test[which(data_test$Group == levels(data_test$Group)[2]), 2])$p.value)
        
        if (sum(shapiro_p_value <= 0.05) >= 1) {
          "Data not normally distributed! P-value integration non-reliable."
        }
        
        ttest_out = t.test(as.formula(paste(colnames(data_test)[2], " ~ Group", sep="")), data = data_test)
        p_values[i_fragment] = ttest_out$p.value
      }
      
      int_p_value = kostsMethod(data_matrix=t(my_data_tmp[,-1]), p_values = p_values, extra_info=TRUE)$P_test
      
      my_data_out[nrow(my_data_out) + 1, 1] = protein_tmp
      my_data_out[nrow(my_data_out), 2] = int_p_value
      my_data_out[nrow(my_data_out), 3] = p_value_assign(int_p_value)
    }
  } else if (int_level == "Peptide") {
    my_data_out = data.frame(matrix(NaN, nrow = 0, ncol = 3))
    colnames(my_data_out) = c("Peptide","p_value", "Signif")
    for (i_peptide in (1:length(unique(my_data$`Peptide`)))) {
      peptide_tmp = unique(my_data$`Peptide`)[i_peptide]
      my_data_tmp = my_data[which(my_data$`Peptide` == peptide_tmp),]
      my_data_tmp_additional = my_data[, -which(colnames(my_data) %in% grouping_data$Replicate)]
      
      my_data_tmp_IDs = paste(my_data_tmp$`Protein Name`, " | ", my_data_tmp$Peptide,
                              " | ","prec+", my_data_tmp$`Precursor Charge`,
                              " | ", my_data_tmp$`Fragment Ion`, "+", my_data_tmp$`Product Charge`, sep="")
      
      my_data_tmp = my_data_tmp[,-which(colnames(my_data_tmp) %in% colnames(my_data_tmp_additional))]
      rownames(my_data_tmp) = my_data_tmp_IDs
      my_data_tmp = rbind(grouping_data$Condition, my_data_tmp)
      rownames(my_data_tmp)[1] = c("Group")
      my_data_tmp = as.data.frame(t(my_data_tmp))
      
      my_data_tmp[, (2:ncol(my_data_tmp))] = apply(my_data_tmp[, (2:ncol(my_data_tmp))], 2, function(x) as.numeric(as.character(x)))
      my_data_tmp = my_data_tmp %>% janitor::clean_names()
      my_data_tmp$group = as.factor(my_data_tmp$group)
      
      p_values = c()
      for (i_frag in 2:ncol(my_data_tmp)) {
        my_data_tmp_test = data.frame(matrix(data = NA, ncol = 2, nrow = length(my_data_tmp$group)))
        colnames(my_data_tmp_test) = c("Group", "Value")
        my_data_tmp_test$Group = my_data_tmp$group
        my_data_tmp_test$Value = as.numeric(my_data_tmp[, i_frag])
        
        p_values = c(p_values, wilcox.test(Value ~ Group, data=my_data_tmp_test, exact = FALSE)[["p.value"]])
      }
      
      int_p_value = kostsMethod(data_matrix=t(my_data_tmp[,-1]), p_values = p_values, extra_info=TRUE)$P_test
      
      p_lim_0_05 = kostsMethod(data_matrix=t(my_data_tmp[,-1]), p_values = rep(0.05, length(p_values)), extra_info=TRUE)$P_test
      p_lim_0_01 = kostsMethod(data_matrix=t(my_data_tmp[,-1]), p_values = rep(0.01, length(p_values)), extra_info=TRUE)$P_test
      p_lim_0_001 = kostsMethod(data_matrix=t(my_data_tmp[,-1]), p_values = rep(0.001, length(p_values)), extra_info=TRUE)$P_test
      p_lim_0_0001 = kostsMethod(data_matrix=t(my_data_tmp[,-1]), p_values = rep(0.0001, length(p_values)), extra_info=TRUE)$P_test
      
      my_data_out[nrow(my_data_out) + 1, 1] = peptide_tmp
      my_data_out[nrow(my_data_out), 2] = int_p_value
      my_data_out[nrow(my_data_out), 3] = p_value_assign(int_p_value, custom_lims = c(p_lim_0_05, p_lim_0_01, p_lim_0_001, p_lim_0_0001))
    }
  }
  return(my_data_out)
}

data_combination_PERMANOVA <- function(my_data, grouping_data, comb_level = "Protein", plot_PCoA = FALSE) {
  library(vegan)
  flag_failed_homogeneity = 0
  
  if (comb_level == "Protein") {
    my_data_out = data.frame(matrix(NaN, nrow = 0, ncol = 3))
    colnames(my_data_out) = c("Protein Name","p_value", "Signif")
    for (i_protein in (1:length(unique(my_data$`Protein Name`)))) {
      protein_tmp = unique(my_data$`Protein Name`)[i_protein]
      my_data_tmp = my_data[which(my_data$`Protein Name` == protein_tmp),]
      my_data_tmp_additional = my_data[, -which(colnames(my_data) %in% grouping_data$Replicate)]
      
      my_data_tmp_IDs = paste(my_data_tmp$`Protein Name`, " | ", my_data_tmp$Peptide,
                              " | ","prec+", my_data_tmp$`Precursor Charge`,
                              " | ", my_data_tmp$`Fragment Ion`, "+", my_data_tmp$`Product Charge`, sep="")
      
      my_data_tmp = my_data_tmp[,-which(colnames(my_data_tmp) %in% colnames(my_data_tmp_additional))]
      rownames(my_data_tmp) = my_data_tmp_IDs
      sample_ids = which(grouping_data$Replicate %in% colnames(my_data_tmp))
      my_data_tmp = rbind(grouping_data$Condition[sample_ids], my_data_tmp)
      rownames(my_data_tmp)[1] = c("grouping_variable")
      my_data_tmp = as.data.frame(t(my_data_tmp))
      my_data_tmp[, (2:ncol(my_data_tmp))] = apply(my_data_tmp[, (2:ncol(my_data_tmp))], 2, function(x) as.numeric(as.character(x)))
      my_data_tmp = my_data_tmp %>% janitor::clean_names()
      my_data_tmp[,1] = factor(my_data_tmp[,1])
      colnames(my_data_tmp)[1] = c("grouping_variable")
      
      dependent_vars = as.matrix(my_data_tmp[,2:ncol(my_data_tmp)])
      independent_var = my_data_tmp$grouping_variable
      
      distance_matrix = vegdist(dependent_vars, method = "euclidean")
      dispersion_test = betadisper(distance_matrix, group = independent_var)
      anova_dispersion = anova(dispersion_test)
      
      if (flag_failed_homogeneity == 0) {
        if (anova_dispersion$`Pr(>F)`[1] < 0.05) {
          # writeLines("FAILED: homogeneity of multivariate dispersion for some of the proteins!")
          # writeLines("")
          # print(anova_dispersion)
          # pairwise_dispersion = permutest(dispersion_test, pairwise = TRUE)
          # print(pairwise_dispersion)
          # plot(dispersion_test, main="")
          flag_failed_homogeneity = 1
        }
      }
      
      if (plot_PCoA == TRUE) {
        plot(dispersion_test, main="")
      }
      
      adonis_result = adonis2(dependent_vars ~ grouping_variable,
                              data = my_data_tmp, permutations = 15000)
      my_data_out[nrow(my_data_out) + 1, 1] = protein_tmp
      my_data_out[nrow(my_data_out), 2] = adonis_result$`Pr(>F)`[1]
      my_data_out[nrow(my_data_out), 3] = p_value_assign(adonis_result$`Pr(>F)`[1])
    }
  } else if (comb_level == "Peptide") {
    my_data_out = data.frame(matrix(NaN, nrow = 0, ncol = 3))
    colnames(my_data_out) = c("Peptide","p_value", "Signif")
    for (i_peptide in (1:length(unique(my_data$`Peptide`)))) {
      peptide_tmp = unique(my_data$`Peptide`)[i_peptide]
      my_data_tmp = my_data[which(my_data$`Peptide` == peptide_tmp),]
      my_data_tmp_additional = my_data[, -which(colnames(my_data) %in% grouping_data$Replicate)]
      
      my_data_tmp_IDs = paste(my_data_tmp$`Protein Name`, " | ", my_data_tmp$Peptide,
                              " | ","prec+", my_data_tmp$`Precursor Charge`,
                              " | ", my_data_tmp$`Fragment Ion`, "+", my_data_tmp$`Product Charge`, sep="")
      
      my_data_tmp = my_data_tmp[,-which(colnames(my_data_tmp) %in% colnames(my_data_tmp_additional))]
      rownames(my_data_tmp) = my_data_tmp_IDs
      sample_ids = which(grouping_data$Replicate %in% colnames(my_data_tmp))
      my_data_tmp = rbind(grouping_data$Condition[sample_ids], my_data_tmp)
      rownames(my_data_tmp)[1] = c("grouping_variable")
      my_data_tmp = as.data.frame(t(my_data_tmp))
      my_data_tmp[, (2:ncol(my_data_tmp))] = apply(my_data_tmp[, (2:ncol(my_data_tmp))], 2, function(x) as.numeric(as.character(x)))
      my_data_tmp = my_data_tmp %>% janitor::clean_names()
      my_data_tmp[,1] = factor(my_data_tmp[,1])
      colnames(my_data_tmp)[1] = c("grouping_variable")
      
      dependent_vars = as.matrix(my_data_tmp[,2:ncol(my_data_tmp)])
      independent_var = my_data_tmp$grouping_variable
      
      distance_matrix = vegdist(dependent_vars, method = "euclidean")
      dispersion_test = betadisper(distance_matrix, group = independent_var)
      anova_dispersion = anova(dispersion_test)
      
      if (flag_failed_homogeneity == 0) {
        if (anova_dispersion$`Pr(>F)`[1] < 0.05) {
          # writeLines("FAILED: homogeneity of multivariate dispersion for some of the proteins!")
          # writeLines("")
          # print(anova_dispersion)
          # pairwise_dispersion = permutest(dispersion_test, pairwise = TRUE)
          # print(pairwise_dispersion)
          # plot(dispersion_test, main="")
          flag_failed_homogeneity = 1
        }
      }
      
      if (plot_PCoA == TRUE) {
        plot(dispersion_test, main="")
      }
      
      adonis_result = adonis2(dependent_vars ~ grouping_variable,
                              data = my_data_tmp, permutations = 15000)
      my_data_out[nrow(my_data_out) + 1, 1] = peptide_tmp
      my_data_out[nrow(my_data_out), 2] = adonis_result$`Pr(>F)`[1]
      my_data_out[nrow(my_data_out), 3] = p_value_assign(adonis_result$`Pr(>F)`[1])
    }
  }
  return(my_data_out)
}

print_calibration_curves <- function(my_data_function, grouping_data, fit, 
                                     concentration_lower_limit, concentration_upper_limit,
                                     skip_print_plot = FALSE) {
  
  library(ggplot2, quietly=T, warn.conflicts=F)
  library(yardstick, quietly=T, warn.conflicts=F)
  library(grid, quietly=T, warn.conflicts=F)
  library(plotly, quietly=T, warn.conflicts=F)
  library(blandr, quietly=T, warn.conflicts=F)
  
  R2_vector = c()
  Pearsons_r_vector = c()
  Lins_CCC_vector = c()
  MAE_vector = c()
  RMSE_vector = c()
  MAPE_vector = c()
  
  for (i_protein in 1:nrow(my_data_function)) {
    protein_tmp = my_data_function$'Protein'[i_protein]
    
    concentrations = grouping_data$Condition
    my_data_function_tmp = data.frame(matrix(NaN, nrow = length(concentrations), ncol = 2))
    colnames(my_data_function_tmp) = c("Concentration", "Intensity")
    my_data_function_tmp$Concentration = (as.numeric(concentrations) / 31.5) * 5
    my_data_function_tmp$Intensity = as.numeric(my_data_function[i_protein, 2:ncol(my_data_function)])
    
    if (fit == "quadratic") {
      p = ggplot(data = my_data_function_tmp, aes(x = Concentration, y = Intensity)) +
        geom_point(alpha = 1) + ylab("Normalized Intensity") +
        xlab("Protein loading [μg per injection]") +
        ggtitle(paste("Calibration curve: ", protein_tmp, sep="")) +
        theme(plot.title = element_text(size=14, face="bold", hjust = 0.5), axis.text.x = element_text(color="black", size=12),
              axis.ticks = element_line(color = "black"), axis.text.y = element_text(color="black", size=12), 
              axis.title = element_text(size=12)) + theme(legend.position = "none")
      
      model = lm(Concentration ~ poly(Intensity,2), data = my_data_function_tmp)
      R2 = summary(model)$adj.r.squared
      Peasrons_r = cor(my_data_function_tmp$Concentration, my_data_function_tmp$Intensity, method = c("pearson"))
      
      p = p + geom_smooth(method = "lm", formula = y~poly(x,2), col="black", linewidth = 1, alpha = 0.3) +
        annotate("Text", x = -Inf, y = Inf, label = sprintf("R^2 == '%.4f'", round(R2,4)), 
                 vjust = 1.5, hjust = -0.1, size=5, parse = TRUE)
      
      if (skip_print_plot != TRUE) {
        print(p)
      }
      R2_vector = c(R2_vector, R2)
      Pearsons_r_vector = c(Pearsons_r_vector, Pearsons_r)
      
    } else if (fit == "sigmoid") {
      my_data_function_tmp$Concentration = log(my_data_function_tmp$Concentration)
      my_data_function_tmp$Intensity = my_data_function_tmp$Intensity
      
      fm1 = nls(Intensity ~ SSlogis(Concentration, a, b, c), data = my_data_function_tmp)
      
      p = ggplot(data = my_data_function_tmp, aes(x = Concentration, y = Intensity)) +
        geom_point(alpha = 1) + ylab("Normalized Intensity") +
        xlab("Protein loading [μg per injection]") +
        ggtitle(paste("Calibration curve: ", protein_tmp, sep="")) +
        theme(plot.title = element_text(size=14, face="bold", hjust = 0.5), axis.text.x = element_text(color="black", size=12),
              axis.ticks = element_line(color = "black"), axis.text.y = element_text(color="black", size=12), 
              axis.title = element_text(size=12)) + theme(legend.position = "none")
      p = p + geom_smooth(method = "nls", se = FALSE,
                          formula = y ~ a/(1+exp(-b*(x-c))),
                          method.args = list(start = coef(fm1), algorithm = 'port'),
                          col="black", linewidth = 1, alpha = 0.3) 
      if (skip_print_plot != TRUE) {
        print(p)
      }
      
    } else if (fit == "find_linear" | fit == "find_linear_log") {
      
      blanks_tmp = my_data_function_tmp[which(my_data_function_tmp$Concentration == 0),]
      LOD_int = 3*sd(blanks_tmp$Intensity)
      LOQ_int = 10*sd(blanks_tmp$Intensity)
      
      lowest_above_LOQ = ((my_data_function_tmp$Concentration[min(which(my_data_function_tmp$Intensity > LOQ_int))] /5) * 31.5)
      writeLines(paste(protein_tmp, " | First concentration above LOQ = ", as.character(lowest_above_LOQ), sep=""))
      
      which_linear = intersect(which(my_data_function_tmp$Concentration >= ((concentration_lower_limit / 31.5) * 5)),
                               which(my_data_function_tmp$Concentration <= ((concentration_upper_limit / 31.5) * 5)))
      which_not_linear = setdiff((1:length(my_data_function_tmp$Concentration)), which_linear)
      my_data_function_tmp_not_linear = my_data_function_tmp[which_not_linear,]
      my_data_function_tmp_linear = my_data_function_tmp[which_linear,]
      
      lm_model = lm(formula = Intensity ~ Concentration, data = my_data_function_tmp_linear)
      R2 = summary(lm_model)$adj.r.squared
      Peasrons_r = cor(my_data_function_tmp$Concentration, my_data_function_tmp$Intensity, method = c("pearson"))
      
      my_data_function_tmp_model = my_data_function_tmp_linear
      my_data_function_tmp_model$Intensity = lm_model$fitted.values
      
      p = ggplot(data = my_data_function_tmp_linear, aes(x = Concentration, y = Intensity)) +
        geom_point() + ylab("Intensity") + xlab("Protein loading [μg per injection]") +
        ggtitle(paste("Calibration curve: ", protein_tmp, sep="")) +
        theme(plot.title = element_text(size=14, face="bold", hjust = 0.5), axis.text.x = element_text(color="black", size=12),
              axis.ticks = element_line(color = "black"), axis.text.y = element_text(color="black", size=12), 
              axis.title = element_text(size=12))
      
      if (fit == "find_linear_log") {
        p = p + scale_y_continuous(trans = scales::log_trans(), breaks = scales::log_breaks()) +
          scale_x_continuous(trans = scales::log_trans(), breaks = scales::log_breaks()) +
          ylab("log(Intensity)") + xlab("log(Protein loading [μg per injection])")
      }
      
      p = p + geom_line(data = my_data_function_tmp_model, aes(x = Concentration, y = Intensity))
      p = p + annotate("Text", x = -Inf, y = Inf, label = sprintf("R^2 == '%.4f'", round(R2,4)), 
                       vjust = 1.5, hjust = -0.1, size=5, parse = TRUE)
      
      p = p + geom_point(data = my_data_function_tmp_not_linear, aes(x = Concentration, y = Intensity),
                         shape=18, color="red", size = 2)
      
      p = p + geom_hline(yintercept = LOD_int, linetype="dashed", color = "red") + 
        geom_hline(yintercept = LOQ_int, linetype="dashed", color = "darkorange")
      
      if (skip_print_plot != TRUE) {
        print(p)
      }
      R2_vector = c(R2_vector, R2)
      Pearsons_r_vector = c(Pearsons_r_vector, Pearsons_r)
      
    } else {
      p = ggplot(data = my_data_function_tmp, aes(x = Concentration, y = Intensity)) +
        geom_point(alpha = 1) + ylab("Normalized Intensity") +
        xlab("Protein loading [μg per injection]") +
        ggtitle(paste("Calibration curve: ", protein_tmp, sep="")) +
        theme(plot.title = element_text(size=14, face="bold", hjust = 0.5), axis.text.x = element_text(color="black", size=12),
              axis.ticks = element_line(color = "black"), axis.text.y = element_text(color="black", size=12), 
              axis.title = element_text(size=12)) + theme(legend.position = "none")
      
      # model = lm(Concentration ~ Intensity, data = my_data_function_tmp)
      model = lm(Intensity ~ Concentration, data = my_data_function_tmp)
      R2 = summary(model)$adj.r.squared
      Pearsons_r = cor(my_data_function_tmp$Concentration, my_data_function_tmp$Intensity, method = c("pearson"))
      
      p = p + geom_smooth(method = "lm", formula = y~x, col="black", linewidth = 1, alpha = 0.3) 
      
      # p = p + annotate("Text", x = -Inf, y = Inf, label = sprintf("R^2 == '%.4f'", round(R2,4)), 
      #            vjust = 1.5, hjust = -0.1, size=5, parse = TRUE)
      
      if (skip_print_plot != TRUE) {
        print(p)
      }
      R2_vector = c(R2_vector, R2)
      Pearsons_r_vector = c(Pearsons_r_vector, Pearsons_r)

      predicted_x = c()
      for (i_predict in 1:length(my_data_function_tmp$Concentration)) {
        predicted_x = c(predicted_x, (my_data_function_tmp$Concentration[i_predict] * as.numeric(model$coefficients[2]) + as.numeric(model$coefficients[1])))
      }
      
      MAE_vector = c(MAE_vector, Metrics::mae(my_data_function_tmp$Intensity, predicted_x) / mean(my_data_function_tmp$Intensity))
      RMSE_vector = c(RMSE_vector, Metrics::rmse(my_data_function_tmp$Intensity, predicted_x) / mean(my_data_function_tmp$Intensity))
      MAPE_vector = c(MAPE_vector, Metrics::mape(my_data_function_tmp$Intensity, predicted_x))
      
      Lins_CCC = CCC(my_data_function_tmp$Intensity, predicted_x)
      Lins_CCC_vector = c(Lins_CCC_vector, as.numeric(Lins_CCC$rho.c[1]))

      #BA_summary_data = blandr.draw(my_data_function_tmp$Intensity, predicted_x, sig.level = 0.95)

    }
  }
  
  return(list("R2" = R2_vector, "Pearsons_r" = Pearsons_r_vector, "Lins_CCC" = Lins_CCC_vector,
              "MAE" = MAE_vector, "RMSE" = RMSE_vector, "MAPE" = MAPE_vector))
}

print_plot_Fig2 <- function(results_all, plot_name, ylim_lower, ylim_upper, path, NAs_proc) {
  results_summary = data.frame(matrix(data = NA, nrow = 0, ncol = 4))
  colnames(results_summary) = c("Method", "Protein", "result_mean", "result_sd")
  for (i_method in 1:length(results_all)) {
    tmp_data_mean = as.data.frame(colMeans(results_all[[i_method]]))
    tmp_data_sd = as.data.frame(colSds(results_all[[i_method]]))
    
    if (grepl("no_MDI", names(results_all)[i_method]) == 1) {
      fig_name = "no MDI"
    } else if (grepl("kNN_", names(results_all)[i_method]) == 1) {
      fig_name = "Classical kNN"
    } else if (grepl("DD", names(results_all)[i_method]) == 1) {
      fig_name = "Data-driven MDI"
    } else if (grepl("AI", names(results_all)[i_method]) == 1) {
      fig_name = "AI MDI"
    } else if (grepl("MSstats", names(results_all)[i_method]) == 1) {
      fig_name = "MSstats MDI"
    }
    tmp_data = cbind(rep(fig_name, 3), c("ADH1", "HXKB", "ENO1"), tmp_data_mean, tmp_data_sd)
    colnames(tmp_data) = c("Method", "Protein", "result_mean", "result_sd")
    results_summary = rbind(results_summary, tmp_data)
  }
  results_summary$Method = factor(results_summary$Method, levels = c("no MDI", "Classical kNN", "Data-driven MDI", "MSstats MDI", "AI MDI"))
  
  p = ggplot(data = results_summary, aes(x = Method, y = result_mean, group = Protein, color=Protein)) +
    geom_line() + geom_point() + ylab(plot_name) + xlab("Method") +
    geom_errorbar(aes(ymin = result_mean - result_sd,
                      ymax = result_mean + result_sd), width = 0.05) +
    theme(legend.position = "none", plot.title = element_text(size=18, face="bold", hjust = 0.5), axis.text.x = element_text(color="black", size=14),
          axis.ticks = element_line(color = "black"), axis.text.y = element_text(color="black", size=14),
          axis.title = element_text(size=14))
  if (max(results_summary$result_mean + results_summary$result_sd, na.rm=TRUE) > 1) {
    p = p + coord_cartesian(ylim = c(NA, 1)) + scale_y_continuous(expand = expand_scale(mult=c(0.03, 0)))
  }
  
  if (!is.na(ylim_lower) && !is.na(ylim_upper)) {
    p = p + ylim(ylim_lower, ylim_upper)
  }
  
  if (NAs_proc == 0) {
    p = p + ggtitle(paste("Missing data proportion: 0.33%", sep=""))
    tiff_path = paste(path, "figures/", clean_filename(plot_name), "_MDI_0_33.tiff", sep="")
    emf_path = paste(path, "figures/", clean_filename(plot_name), "_MDI_0_33.emf", sep="")
  } else {
    p = p + ggtitle(paste("Missing data proportion: ", as.character(round(NAs_proc * 100, 2)), "%", sep=""))
    tiff_path = paste(path, "figures/", clean_filename(plot_name), "_MDI_", as.character(round(NAs_proc * 100, 2)), ".tiff", sep="")
    emf_path = paste(path, "figures/", clean_filename(plot_name), "_MDI_", as.character(round(NAs_proc * 100, 2)), ".emf", sep="")
  }
  p
  
  ggsave(plot = p, filename = tiff_path, width = 20.8, height = 13.6, units = "cm", dpi = "retina")
  ggsave(plot = p, filename = emf_path, width = 20.8, height = 13.6, units = "cm", dpi = "retina")
}

library(ggrepel, quietly=T, warn.conflicts=F)
library(readxl, quietly=T, warn.conflicts=F)
library(janitor, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
library(RColorBrewer, quietly=T, warn.conflicts=F)
library(scales, quietly=T, warn.conflicts=F)
library(reshape2, quietly=T, warn.conflicts=F)
library(plyr, quietly=T, warn.conflicts=F)
library(dplyr, quietly=T, warn.conflicts=F)
library(matrixStats, quietly=T, warn.conflicts=F)
library(progress, quietly=T, warn.conflicts=F)
library(SiZer, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
library(yardstick, quietly=T, warn.conflicts=F)
library(grid, quietly=T, warn.conflicts=F)
library(plotly, quietly=T, warn.conflicts=F)
library(collapse, quietly=T, warn.conflicts=F)
library(DescTools, quietly=T, warn.conflicts=F)
library(Rmisc, quietly=T, warn.conflicts=F)
library(fmsb, quietly=T, warn.conflicts=F)
library(MSstats, quietly=T, warn.conflicts=F)
library(writexl, quietly=T, warn.conflicts=F)
library(missForest, quietly=T, warn.conflicts=F)
library(devEMF, quietly=T, warn.conflicts=F)


######################################################################################
# Set path
######################################################################################

# This codes requires input_data.xlsx and "figures" folder in the path to save output
path = "C:/"
data_path = paste(path, "input_data.xlsx", sep="")
setwd(path)

# ######################################################################################
# # Figure 1 - Dataset characteristics - use data for yeast only
# ######################################################################################
set_kNN_k_parameter = 3
TP_proteins = c("sp|P00330|ADH1_YEAST", "sp|P04807|HXKB_YEAST", "sp|P00924|ENO1_YEAST")

my_data = as.data.frame(read_excel(data_path))
grouping_data = get_grouping_data(my_data)
my_data = my_data[my_data$Protein %in% TP_proteins, ]

my_data_prepared = my_data_input_prepare(my_data, grouping_data, n_first = 3, n_last = 3,
                                         rm_precursors = TRUE,
                                         rm_zero_na = FALSE,
                                         rm_blanks = FALSE)
my_data = my_data_prepared[[1]]
my_data_additional = my_data_prepared[[2]]
grouping_data = my_data_prepared[[3]]
labels = my_data_prepared[[4]]
metabolite_names = my_data_prepared[[5]]

my_data = my_data[,-which(colnames(my_data) == "Group")]
my_data = cbind(my_data_additional,  t(my_data))

out_data_combination = data_combination_summing(my_data, grouping_data)
my_data = out_data_combination

quan_limits = data.frame(matrix(NaN, nrow = length(my_data$'Protein'), ncol = 4))
colnames(quan_limits) = c("Protein", "LOD", "LOQ", "LOL")
quan_limits$Protein = my_data$'Protein'

quan_limit_conc = data.frame(matrix(NaN, nrow = length(my_data$'Protein'), ncol = 4))
colnames(quan_limit_conc) = c("Protein", "LOD", "LOQ", "LOL")
quan_limit_conc$Protein = my_data$'Protein'

picewise_linear_fit = list()
for (i_protein in 1:nrow(my_data)) {
  protein_tmp = my_data$'Protein'[i_protein]

  concentrations = grouping_data$Condition
  my_data_tmp = data.frame(matrix(NaN, nrow = length(concentrations), ncol = 2))
  colnames(my_data_tmp) = c("Concentration", "Intensity")
  my_data_tmp$Concentration = (as.numeric(concentrations) / 31.5) * 5
  my_data_tmp$Intensity = as.numeric(my_data[i_protein, 2:ncol(my_data)])

  blanks_tmp = my_data_tmp[which(my_data_tmp$Concentration == 0),]
  LOD_tmp = 3 * sd(blanks_tmp$Intensity)
  LOQ_tmp  = 10 * sd(blanks_tmp$Intensity)
  picewise_linear_data = piecewise.linear(my_data_tmp$Concentration, my_data_tmp$Intensity)
  LOL_tmp = picewise_linear_data[["change.point"]]

  picewise_linear_data_x = picewise_linear_data$x
  picewise_linear_data_y = picewise_linear_data$y
  picewise_linear_data = data.frame(matrix(NaN, nrow = length(picewise_linear_data$x), ncol = 2))
  picewise_linear_data[,1] = picewise_linear_data_x
  picewise_linear_data[,2] = picewise_linear_data_y
  picewise_linear_fit[[protein_tmp]] = picewise_linear_data

  quan_limits$LOD[which(quan_limits$Protein == protein_tmp)] = LOD_tmp
  quan_limits$LOQ[which(quan_limits$Protein == protein_tmp)]  = LOQ_tmp
  quan_limits$LOL[which(quan_limits$Protein == protein_tmp)] = LOL_tmp

  my_data_tmp$LOD_criteria = (my_data_tmp$Intensity < LOD_tmp)
  my_data_tmp$LOQ_criteria = (my_data_tmp$Intensity < LOQ_tmp)
  LOD_criteria_inx = min(which(as.numeric(aggregate(LOD_criteria  ~ Concentration, my_data_tmp, sum)$LOD_criteria) == 0))
  LOQ_criteria_inx = min(which(as.numeric(aggregate(LOQ_criteria  ~ Concentration, my_data_tmp, sum)$LOQ_criteria) == 0))

  my_data_tmp$LOL_criteria = (my_data_tmp$Concentration < LOL_tmp)
  LOL_criteria_inx = max(which(as.numeric(aggregate(LOL_criteria  ~ Concentration, my_data_tmp, sum)$LOL_criteria) != 0))

  quan_limit_conc$LOD[which(quan_limit_conc$Protein == protein_tmp)] = ((sort(unique(my_data_tmp$Concentration))[LOD_criteria_inx] / 5) * 31.5)
  quan_limit_conc$LOQ[which(quan_limit_conc$Protein == protein_tmp)] = ((sort(unique(my_data_tmp$Concentration))[LOQ_criteria_inx] / 5) * 31.5)
  quan_limit_conc$LOL[which(quan_limit_conc$Protein == protein_tmp)] = ((sort(unique(my_data_tmp$Concentration))[LOL_criteria_inx] / 5) * 31.5)
}

for (i_protein in 1:nrow(my_data)) {
  protein_tmp = my_data$'Protein'[i_protein]

  concentrations = grouping_data$Condition
  my_data_tmp = data.frame(matrix(NaN, nrow = length(concentrations), ncol = 2))
  colnames(my_data_tmp) = c("Concentration", "Intensity")
  my_data_tmp$Concentration = (as.numeric(concentrations) / 31.5) * 5
  my_data_tmp$Intensity = as.numeric(my_data[i_protein, 2:ncol(my_data)])

  concentration_LOD = quan_limits$LOD[which(quan_limits$Protein == protein_tmp)]
  concentration_LOQ = quan_limits$LOQ[which(quan_limits$Protein == protein_tmp)]
  concentration_LOL = quan_limits$LOL[which(quan_limits$Protein == protein_tmp)]

  concentration_lower_limit = quan_limit_conc$LOQ[which(quan_limit_conc$Protein == protein_tmp)]
  concentration_upper_limit = quan_limit_conc$LOL[which(quan_limit_conc$Protein == protein_tmp)]

  which_linear_lower = which(my_data_tmp$Concentration >= ((concentration_lower_limit / 31.5) * 5))
  which_linear_upper = which(my_data_tmp$Concentration <= ((concentration_upper_limit / 31.5) * 5))
  which_linear = intersect(which_linear_lower, which_linear_upper)

  which_not_linear = setdiff((1:length(my_data_tmp$Concentration)), which_linear)
  my_data_tmp_not_linear = my_data_tmp[which_not_linear,]
  my_data_tmp_linear = my_data_tmp[which_linear,]

  base_p = ggplot(data = my_data_tmp_linear, aes(x = Concentration, y = Intensity)) +
           geom_point() + ylab("Intensity") + xlab("Protein loading [ng]") +
           ggtitle(paste("Calibration curve: ", protein_tmp, sep="")) +
           theme(plot.title = element_text(size=14, face="bold", hjust = 0.5), axis.text.x = element_text(color="black", size=12),
                 axis.ticks = element_line(color = "black"), axis.text.y = element_text(color="black", size=12),
                 axis.title = element_text(size=12))

  base_p = base_p + geom_point(data = my_data_tmp[setdiff(which_linear_upper, which_linear),],
                          aes(x = Concentration, y = Intensity), shape=16, color="red", size = 1.5)
  base_p = base_p + geom_point(data = my_data_tmp[setdiff(which_linear_lower, which_linear),],
                     aes(x = Concentration, y = Intensity), shape=16, color="black", size = 1.5)
  base_p = base_p + geom_hline(yintercept = concentration_LOD, linetype="longdash", color = "darkred") +
    geom_hline(yintercept = concentration_LOQ, linetype="longdash", color = "darkorange")
  p_log = base_p + scale_y_log10(labels = label_log(), limits = c(9.8e3,3e7)) +
                   scale_x_log10(labels = label_log()) +
                   ylab("Intensity") + xlab("Protein loading [ng]")
  p_upper_lim = base_p + geom_point(data = my_data_tmp[setdiff(which_linear_lower, which_linear),],
                                    aes(x = Concentration, y = Intensity), shape=16, color="red", size = 1.5)
  p_upper_lim = p_upper_lim + geom_point(data = my_data_tmp[setdiff(which_linear_lower, which_linear),],
                                         aes(x = Concentration, y = Intensity), shape=16, color="red", size = 1.5)
  p_upper_lim = p_upper_lim +  geom_hline(yintercept = picewise_linear_fit[[protein_tmp]]$X2[which.min(abs(picewise_linear_fit[[protein_tmp]]$X1 - concentration_LOL))],
                                          linetype="longdash", color = "red")

  p_upper_lim = p_upper_lim + scale_y_log10(labels = label_log(), limits = c(9.8e3,3e7)) + scale_x_log10(labels = label_log())
  p_upper_lim

  ggsave(plot = p_upper_lim, filename = paste("figures/", "calibration_", strsplit(protein_tmp, "\\|")[[1]][3], ".tiff", sep=""),
         width = 15.3, height = 15.3, units = "cm", dpi = "retina")
  ggsave(plot = p_upper_lim, filename = paste("figures/", "calibration_", strsplit(protein_tmp, "\\|")[[1]][3], ".emf", sep=""),
         width = 15.3, height = 15.3, units = "cm")
}

# ######################################################################################
# # Figure 2 / 3 - Missing data imputation - use reduced dataset (2.5 - 25ng per 31.5uL)
# ######################################################################################
NAs_proportion = c(0, 0.01, 0.05, 0.1)
random_NAs_sampling = 100
set_kNN_k_parameter = 3
TP_proteins = c("sp|P00330|ADH1_YEAST", "sp|P04807|HXKB_YEAST", "sp|P00924|ENO1_YEAST")
keep_concentrations_range = c(2.5, 25)

my_data = as.data.frame(read_excel(data_path))
grouping_data = get_grouping_data(my_data)
my_data = my_data[my_data$Protein %in% TP_proteins, ]
grouping_data = grouping_data[grouping_data$Condition >= keep_concentrations_range[1] & grouping_data$Condition <= keep_concentrations_range[2],]
columns_to_keep = colnames(my_data) %in% grouping_data$Replicate
columns_to_keep[1:9] = TRUE
my_data = my_data[,columns_to_keep]

my_data_prepared = my_data_input_prepare(my_data, grouping_data, n_first = 3, n_last = 3,
                                         rm_precursors = TRUE,
                                         rm_zero_na = FALSE,
                                         rm_blanks = FALSE)
my_data = my_data_prepared[[1]]
my_data_additional = my_data_prepared[[2]]
grouping_data = my_data_prepared[[3]]
labels = my_data_prepared[[4]]
metabolite_names = my_data_prepared[[5]]

my_data = my_data[,-which(colnames(my_data) == "Group")]
rownames_my_data = rownames(my_data)
my_data = data.frame(sapply(my_data, function(x) as.numeric(as.character(x))))
rownames(my_data) = rownames_my_data
my_data = cbind(my_data_additional,  t(my_data))
my_data_full = my_data

for (i_NAs_proportion in 1:length(NAs_proportion)) {
  final_Pearsons_r = list("no_MDI_Pearsons_r" = c(), "kNN_Pearsons_r" = c(), "DD_Pearsons_r" = c(), "AI_Pearsons_r" = c(), "MSstats_Pearsons_r" = c())
  final_R2 = list("no_MDI_R2" = c(), "kNN_R2" = c(), "DD_R2" = c(), "AI_R2" = c(), "MSstats_R2" = c())
  final_Lins_CCC = list("no_MDI_Lins_CCC" = c(), "kNN_Lins_CCC" = c(), "DD_Lins_CCC" = c(), "AI_Lins_CCC" = c(), "MSstats_Lins_CCC" = c())
  final_MAE = list("no_MDI_MAE" = c(), "kNN_MAE" = c(), "DD_MAE" = c(), "AI_MAE" = c(), "MSstats_MAE" = c())
  final_RMSE = list("no_MDI_RMSE" = c(), "kNN_RMSE" = c(), "DD_RMSE" = c(), "AI_RMSE" = c(), "MSstats_RMSE" = c())
  final_MAPE = list("no_MDI_MAPE" = c(), "kNN_MAPE" = c(), "DD_MAPE" = c(), "AI_MAPE" = c(), "MSstats_MAPE" = c())

  for (i_loop in 1:random_NAs_sampling) {
    writeLines(paste("Loop ", as.character(i_loop), " / ", as.character(random_NAs_sampling), sep = ""))

    my_data = my_data_full
    set.seed(1243 + i_loop)
    my_data[,((ncol(my_data_additional)+1):ncol(my_data))] = missForest::prodNA(my_data[,((ncol(my_data_additional)+1):ncol(my_data))], noNA = NAs_proportion[i_NAs_proportion])

    my_data_no_MDI = my_data[,((ncol(my_data_additional)+1):ncol(my_data))]
    my_data_no_MDI[is.na(my_data_no_MDI)] = 0
    my_data_no_MDI = cbind(my_data_additional, my_data_no_MDI)
    out_data_combination = data_combination_summing(my_data_no_MDI, grouping_data)
    no_MDI_out = print_calibration_curves(out_data_combination, grouping_data, "linear", skip_print_plot = TRUE)

    my_data_new_kNN_return = MDI_kNN(my_data, grouping_data, set_kNN_k_parameter)
    my_data_new_kNN_return = my_data_new_kNN_return["my_data_new"][[1]]
    out_data_combination = data_combination_summing(my_data_new_kNN_return, grouping_data)
    kNN_out = print_calibration_curves(out_data_combination, grouping_data, "linear", skip_print_plot = TRUE)

    my_data_new_DD_return = MDI_data_driven(my_data, grouping_data, show_pb = FALSE)
    my_data_new_DD_return = my_data_new_DD_return["my_data_new"][[1]]
    out_data_combination = data_combination_summing(my_data_new_DD_return, grouping_data)
    DD_out = print_calibration_curves(out_data_combination, grouping_data, "linear", skip_print_plot = TRUE)

    my_data_new_AI_return = MDI_AI(my_data, grouping_data, choose_ref = "CVs", set_Koina_model = "ms2pip_timsTOF2024")
    my_data_new_AI_return = my_data_new_AI_return["my_data_new"][[1]]
    out_data_combination = data_combination_summing(my_data_new_AI_return, grouping_data)
    AI_out = print_calibration_curves(out_data_combination, grouping_data, "linear", skip_print_plot = TRUE)

    my_data_prepared_for_MSstats = MSstats_prepare_own(my_data, grouping_data)
    invisible(capture.output({ my_data_new_MSstats_return = dataProcess(raw = my_data_prepared_for_MSstats, use_log_file = FALSE, normalization = FALSE,
                                                                        featureSubset = "all", MBimpute = TRUE)}))
    my_data_new_MSstats_return_new = MSstats_return_feature_own(my_data_new_MSstats_return, my_data)
    out_data_combination_MSstats = data_combination_summing(my_data_new_MSstats_return_new, grouping_data)
    MSstats_out = print_calibration_curves(out_data_combination_MSstats, grouping_data, "linear", skip_print_plot = TRUE)

    final_Pearsons_r$no_MDI_Pearsons_r = rbind(final_Pearsons_r$no_MDI_Pearsons_r, no_MDI_out$Pearsons_r)
    final_Pearsons_r$kNN_Pearsons_r = rbind(final_Pearsons_r$kNN_Pearsons_r, kNN_out$Pearsons_r)
    final_Pearsons_r$DD_Pearsons_r = rbind(final_Pearsons_r$DD_Pearsons_r, DD_out$Pearsons_r)
    final_Pearsons_r$AI_Pearsons_r = rbind(final_Pearsons_r$AI_Pearsons_r, AI_out$Pearsons_r)
    final_Pearsons_r$MSstats_Pearsons_r = rbind(final_Pearsons_r$MSstats_Pearsons_r, MSstats_out$Pearsons_r)

    final_R2$no_MDI_R2 = rbind(final_R2$no_MDI_R2, no_MDI_out$R2)
    final_R2$kNN_R2 = rbind(final_R2$kNN_R2, kNN_out$R2)
    final_R2$DD_R2 = rbind(final_R2$DD_R2, DD_out$R2)
    final_R2$AI_R2 = rbind(final_R2$AI_R2, AI_out$R2)
    final_R2$MSstats_R2 = rbind(final_R2$MSstats_R2, MSstats_out$R2)

    final_Lins_CCC$no_MDI_Lins_CCC = rbind(final_Lins_CCC$no_MDI_Lins_CCC, no_MDI_out$Lins_CCC)
    final_Lins_CCC$kNN_Lins_CCC = rbind(final_Lins_CCC$kNN_Lins_CCC, kNN_out$Lins_CCC)
    final_Lins_CCC$DD_Lins_CCC = rbind(final_Lins_CCC$DD_Lins_CCC, DD_out$Lins_CCC)
    final_Lins_CCC$AI_Lins_CCC = rbind(final_Lins_CCC$AI_Lins_CCC, AI_out$Lins_CCC)
    final_Lins_CCC$MSstats_Lins_CCC = rbind(final_Lins_CCC$MSstats_Lins_CCC, MSstats_out$Lins_CCC)

    final_MAE$no_MDI_MAE = rbind(final_MAE$no_MDI_MAE, no_MDI_out$MAE)
    final_MAE$kNN_MAE = rbind(final_MAE$kNN_MAE, kNN_out$MAE)
    final_MAE$DD_MAE = rbind(final_MAE$DD_MAE, DD_out$MAE)
    final_MAE$AI_MAE = rbind(final_MAE$AI_MAE, AI_out$MAE)
    final_MAE$MSstats_MAE = rbind(final_MAE$MSstats_MAE, MSstats_out$MAE)

    final_RMSE$no_MDI_RMSE = rbind(final_RMSE$no_MDI_RMSE, no_MDI_out$RMSE)
    final_RMSE$kNN_RMSE = rbind(final_RMSE$kNN_RMSE, kNN_out$RMSE)
    final_RMSE$DD_RMSE = rbind(final_RMSE$DD_RMSE, DD_out$RMSE)
    final_RMSE$AI_RMSE = rbind(final_RMSE$AI_RMSE, AI_out$RMSE)
    final_RMSE$MSstats_RMSE = rbind(final_RMSE$MSstats_RMSE, MSstats_out$RMSE)

    final_MAPE$no_MDI_MAPE = rbind(final_MAPE$no_MDI_MAPE, no_MDI_out$MAPE)
    final_MAPE$kNN_MAPE = rbind(final_MAPE$kNN_MAPE, kNN_out$MAPE)
    final_MAPE$DD_MAPE = rbind(final_MAPE$DD_MAPE, DD_out$MAPE)
    final_MAPE$AI_MAPE = rbind(final_MAPE$AI_MAPE, AI_out$MAPE)
    final_MAPE$MSstats_MAPE = rbind(final_MAPE$MSstats_MAPE, MSstats_out$MAPE)
  }

  list_of_all_lists = c(final_Pearsons_r, final_R2, final_Lins_CCC, final_MAE, final_RMSE, final_MAPE)
  writexl::write_xlsx(lapply(list_of_all_lists, as.data.frame), paste(path, "figures/", "Fig2_data_MDI_",
                                                                      as.character(round(NAs_proportion[i_NAs_proportion] * 100, 2)), ".xlsx", sep=""))

  print_plot_Fig2(final_Pearsons_r, "Pearson's r", 0.965, 0.995, path, NAs_proportion[i_NAs_proportion])
  print_plot_Fig2(final_R2, "R²", 0.935, 0.988, path, NAs_proportion[i_NAs_proportion])
  print_plot_Fig2(final_Lins_CCC, "Lin's CCC", 0.97, 0.992, path, NAs_proportion[i_NAs_proportion])
  print_plot_Fig2(final_MAE, "MAE", 0.05, 0.135, path, NAs_proportion[i_NAs_proportion])
  print_plot_Fig2(final_RMSE, "RMSE", 0.08, 0.16, path, NAs_proportion[i_NAs_proportion])
  print_plot_Fig2(final_MAPE, "MAPE", 0.05, 0.2, path, NAs_proportion[i_NAs_proportion])
}

#####################################################################################
# Figure 5 - Data combination and testing
#####################################################################################
TP_proteins = c("sp|P00330|ADH1_YEAST", "sp|P04807|HXKB_YEAST", "sp|P00924|ENO1_YEAST")
set_kNN_k_parameter = 3
keep_concentrations_range = c(2.5, 25)

my_data = as.data.frame(read_excel(data_path))
grouping_data = get_grouping_data(my_data)
grouping_data = grouping_data[grouping_data$Condition >= keep_concentrations_range[1] & grouping_data$Condition <= keep_concentrations_range[2],]
columns_to_keep = colnames(my_data) %in% grouping_data$Replicate
columns_to_keep[1:9] = TRUE
my_data = my_data[,columns_to_keep]

my_data_prepared = my_data_input_prepare(my_data, grouping_data, n_first = 3, n_last = 3,
                                         rm_precursors = TRUE,
                                         rm_zero_na = TRUE,
                                         rm_blanks = FALSE)
my_data = my_data_prepared[[1]]
my_data_additional = my_data_prepared[[2]]
grouping_data = my_data_prepared[[3]]
labels = my_data_prepared[[4]]
metabolite_names = my_data_prepared[[5]]
my_data = my_data[,-which(colnames(my_data) == "Group")]
save_rownames = rownames(my_data)
my_data = data.frame(sapply(my_data, function(x) as.numeric(as.character(x))))
rownames(my_data) = save_rownames
my_data = cbind(my_data_additional,  t(my_data))

my_data_saved_i_MDI = my_data
grouping_data_saved_i_MDI = grouping_data
for (i_MDI in 1:4) {
  my_data = my_data_saved_i_MDI
  grouping_data = grouping_data_saved_i_MDI
  if (i_MDI == 1) {
    my_data = MDI_kNN(my_data, grouping_data, set_kNN_k_parameter)["my_data_new"][[1]]
    MDI_name = "Classical kNN"
  } else if (i_MDI == 2) {
    my_data = MDI_data_driven(my_data, grouping_data, show_pb = FALSE)["my_data_new"][[1]]
    MDI_name = "Data-driven MDI"
  } else if (i_MDI == 3) {
    my_data_MSstats = MSstats_prepare_own(my_data, grouping_data)
    invisible(capture.output({my_data_new_MSstats_return = dataProcess(raw = my_data_MSstats, use_log_file = FALSE, normalization = FALSE, 
                                                                       featureSubset = "all", MBimpute = TRUE)}))
    my_data = MSstats_return_feature_own(my_data_new_MSstats_return, my_data)
    MDI_name = "MSstats MDI"
  } else if (i_MDI == 4) {
    my_data = MDI_AI(my_data, grouping_data, set_Koina_model = "ms2pip_timsTOF2024")["my_data_new"][[1]]
    MDI_name = "AI MDI"
    my_data_additional = my_data_additional[-which(my_data_additional$`Product Charge` == 2),]
  }

  all_group_combinations = combn(unique(grouping_data$Condition), m=2)
  all_group_combinations = all_group_combinations[, order(as.numeric(abs(all_group_combinations[1,] - all_group_combinations[2,]))), drop = FALSE]

  out_test = data.frame(matrix(data = NA, ncol = 2+ncol(all_group_combinations), nrow = (length(my_data$Peptide))))
  colnames(out_test) = c("Protein", "Peptide", paste(all_group_combinations[1,], "_", all_group_combinations[2,], sep=""))
  out_test$Protein = my_data$`Protein Name`
  out_test$Peptide = my_data$Peptide
  out_test = distinct(out_test)
  out_test = list("summing" = out_test, "best_peak_global_max" = out_test, "best_peak_local_max" = out_test,
                  "AI_highest" = out_test, "AI_heaviest" = out_test,
                  "p_value_integration" = out_test, "multivariate_testing" = out_test,
                  "MSstats_lm" = out_test, "MSstats_TMP" = out_test)

  out_contab = data.frame(matrix(data = NA, ncol = ncol(all_group_combinations), nrow = 4))
  colnames(out_contab) = c(paste(all_group_combinations[1,], "_", all_group_combinations[2,], sep=""))
  rownames(out_contab) = c("TP", "FN", "FP", "TN")
  out_contab = list("summing" = out_contab, "best_peak_global_max" = out_contab, "best_peak_local_max" = out_contab,
                    "AI_highest" = out_contab, "AI_heaviest" = out_contab,
                    "p_value_integration" = out_contab, "multivariate_testing" = out_contab,
                    "MSstats_lm" = out_contab, "MSstats_TMP" = out_contab)

  grouping_data_saved = grouping_data

  for (i_iter in 1:ncol(all_group_combinations)) {
    writeLines(paste("Loop ", as.character(i_iter), " / ", as.character(ncol(all_group_combinations)) ,sep = ""))
    writeLines("")

    grouping_data = grouping_data_saved

    my_data_core = my_data[,which(colnames(my_data) %in% grouping_data$Replicate)]
    my_data_core_1 = my_data_core[, which(grouping_data$Condition == all_group_combinations[1,i_iter])]
    my_data_core_2 = my_data_core[, which(grouping_data$Condition == all_group_combinations[2,i_iter])]
    my_data_tmp = cbind(my_data_additional, my_data_core_1, my_data_core_2)

    grouping_data_1 = grouping_data[which(grouping_data$Condition == all_group_combinations[1,i_iter]),]
    grouping_data_2 = grouping_data[which(grouping_data$Condition == all_group_combinations[2,i_iter]),]
    grouping_data = rbind(grouping_data_1, grouping_data_2)

    out_data_sum = data_combination_summing(my_data_tmp, grouping_data, sum_level = "Peptide")
    out_test_sum_tmp = U_Mann_Whitney_peptide_test(out_data_sum, grouping_data)

    out_data_best_global = data_combination_choose_best(my_data_tmp, grouping_data, sum_level = "Peptide", max_type = "global")
    out_test_best_global = U_Mann_Whitney_peptide_test(out_data_best_global, grouping_data)

    out_data_best_local = data_combination_choose_best(my_data_tmp, grouping_data, sum_level = "Peptide", max_type = "local")
    out_test_best_local = U_Mann_Whitney_peptide_test(out_data_best_local, grouping_data)

    out_data_AI_highest = data_rescaling_AI_without_change(my_data_tmp, grouping_data, show_double_charge_warn = FALSE,
                                                           choose_ref = "highest", consolidation_method = "Peptide")
    out_data_AI_highest = data_combination_summing(out_data_AI_highest, grouping_data, sum_level = "Peptide")
    out_test_AI_highest_tmp = U_Mann_Whitney_peptide_test(out_data_AI_highest, grouping_data)

    out_data_AI_heaviest = data_rescaling_AI_without_change(my_data_tmp, grouping_data, show_double_charge_warn = FALSE,
                                                            choose_ref = "heaviest", consolidation_method = "Peptide")
    out_data_AI_heaviest = data_combination_summing(out_data_AI_heaviest, grouping_data, sum_level = "Peptide")
    out_test_AI_heaviest_tmp = U_Mann_Whitney_peptide_test(out_data_AI_heaviest, grouping_data)

    out_data_int = data_combination_p_int_Browns(my_data_tmp, grouping_data, int_level = "Peptide")

    out_data_multivariate = data_combination_PERMANOVA(my_data_tmp, grouping_data, comb_level = "Peptide")

    my_data_tmp_MSstats = MSstats_prepare_own(my_data_tmp, grouping_data, level = "Peptide")
    invisible(capture.output({ out_data_MSstats_lm = dataProcess(raw = my_data_tmp_MSstats, use_log_file = FALSE, summaryMethod = 'linear', normalization = FALSE, 
                                                                 featureSubset = "all", MBimpute = FALSE)}))
    out_data_MSstats_lm = MSstats_return_protein_own(out_data_MSstats_lm)
    colnames(out_data_MSstats_lm)[1] = "Peptide"
    out_data_MSstats_lm_tmp = U_Mann_Whitney_peptide_test(out_data_MSstats_lm, grouping_data)
    out_data_MSstats_lm_tmp = out_data_MSstats_lm_tmp[match(out_test_sum_tmp$Peptide, out_data_MSstats_lm_tmp$Peptide),]

    invisible(capture.output({out_data_MSstats_TMP = dataProcess(raw = my_data_tmp_MSstats, use_log_file = FALSE, summaryMethod = 'TMP', normalization = FALSE, 
                                                                 featureSubset = "all", MBimpute = FALSE)}))
    out_data_MSstats_TMP = MSstats_return_protein_own(out_data_MSstats_TMP)
    colnames(out_data_MSstats_TMP)[1] = "Peptide"
    out_data_MSstats_TMP_tmp = U_Mann_Whitney_peptide_test(out_data_MSstats_TMP, grouping_data)
    out_data_MSstats_TMP_tmp = out_data_MSstats_TMP_tmp[match(out_test_sum_tmp$Peptide, out_data_MSstats_TMP_tmp$Peptide),]

    out_test$summing[,which(colnames(out_test$summing) %in% paste(all_group_combinations[1,i_iter], "_", all_group_combinations[2,i_iter], sep=""))] = out_test_sum_tmp$signif
    out_test$best_peak_global_max[,which(colnames(out_test$summing) %in% paste(all_group_combinations[1,i_iter], "_", all_group_combinations[2,i_iter], sep=""))] = out_test_best_global$signif
    out_test$best_peak_local_max[,which(colnames(out_test$summing) %in% paste(all_group_combinations[1,i_iter], "_", all_group_combinations[2,i_iter], sep=""))] = out_test_best_local$signif
    out_test$AI_highest[,which(colnames(out_test$summing) %in% paste(all_group_combinations[1,i_iter], "_", all_group_combinations[2,i_iter], sep=""))] = out_test_AI_highest_tmp$signif
    out_test$AI_heaviest[,which(colnames(out_test$summing) %in% paste(all_group_combinations[1,i_iter], "_", all_group_combinations[2,i_iter], sep=""))] = out_test_AI_heaviest_tmp$signif
    out_test$p_value_integration[,which(colnames(out_test$summing) %in% paste(all_group_combinations[1,i_iter], "_", all_group_combinations[2,i_iter], sep=""))] = out_data_int$Signif
    out_test$multivariate_testing[,which(colnames(out_test$summing) %in% paste(all_group_combinations[1,i_iter], "_", all_group_combinations[2,i_iter], sep=""))] = out_data_multivariate$Signif
    out_test$MSstats_lm[,which(colnames(out_test$summing) %in% paste(all_group_combinations[1,i_iter], "_", all_group_combinations[2,i_iter], sep=""))] = out_data_MSstats_lm_tmp$signif
    out_test$MSstats_TMP[,which(colnames(out_test$summing) %in% paste(all_group_combinations[1,i_iter], "_", all_group_combinations[2,i_iter], sep=""))] = out_data_MSstats_TMP_tmp$signif
  }

  out_contab_summary_tmp = data.frame(matrix(data = NA, ncol = ncol(all_group_combinations), nrow = length(out_contab)))
  colnames(out_contab_summary_tmp) = c(paste(all_group_combinations[1,], "_", all_group_combinations[2,], sep=""))
  rownames(out_contab_summary_tmp) = c("Mathematical sum", "Best peak - global max", "Best peak - local max", "AI scaling - highest fragment",
                                   "AI scaling - heaviest fragment", "p-value integration", "Multivariate testing", "MSstats - linear", "MSstats - TMP")

  out_contab_summary = list("FDR" = out_contab_summary_tmp, "Sensitivity" = out_contab_summary_tmp, "Specificity" = out_contab_summary_tmp,
                            "Accuracy" = out_contab_summary_tmp, "Precision" = out_contab_summary_tmp)

  for (i_contab in 1:length(out_contab)) {
    tmp_table = out_test[[i_contab]]
    for (i_col in 3:ncol(tmp_table)) {
      signif_pos = which(tmp_table[,i_col] != "ns")
      non_signif_pos = which(tmp_table[,i_col] == "ns")

      actual_signif_pos = which(tmp_table[,1] %in% TP_proteins)
      actual_non_signif_pos = (1:length(tmp_table[,i_col]))[! (1:length(tmp_table[,i_col])) %in% actual_signif_pos]

      out_contab[[i_contab]][1, (i_col - 2)] = length(which(signif_pos %in% actual_signif_pos)) # TP
      out_contab[[i_contab]][2, (i_col - 2)] = length(which(non_signif_pos %in% actual_signif_pos)) # FN
      out_contab[[i_contab]][3, (i_col - 2)] = length(which(signif_pos %in% actual_non_signif_pos)) # FP
      out_contab[[i_contab]][4, (i_col - 2)] = length(which(non_signif_pos %in% actual_non_signif_pos)) # TN

      out_contab_summary$FDR[i_contab, (i_col - 2)] = out_contab[[i_contab]][3, (i_col - 2)] / (out_contab[[i_contab]][3, (i_col - 2)] + out_contab[[i_contab]][1, (i_col - 2)])
      out_contab_summary$Sensitivity[i_contab, (i_col - 2)] = out_contab[[i_contab]][1, (i_col - 2)] / (out_contab[[i_contab]][1, (i_col - 2)] + out_contab[[i_contab]][2, (i_col - 2)])
      out_contab_summary$Specificity[i_contab, (i_col - 2)] = out_contab[[i_contab]][4, (i_col - 2)] / (out_contab[[i_contab]][3, (i_col - 2)] + out_contab[[i_contab]][4, (i_col - 2)])
      out_contab_summary$Accuracy[i_contab, (i_col - 2)] = (out_contab[[i_contab]][1, (i_col - 2)] + out_contab[[i_contab]][4, (i_col - 2)]) / (out_contab[[i_contab]][1, (i_col - 2)] +
                                                            out_contab[[i_contab]][4, (i_col - 2)] + out_contab[[i_contab]][3, (i_col - 2)] + out_contab[[i_contab]][2, (i_col - 2)])
      out_contab_summary$Precision[i_contab, (i_col - 2)] = out_contab[[i_contab]][1, (i_col - 2)] / (out_contab[[i_contab]][1, (i_col - 2)] + out_contab[[i_contab]][3, (i_col - 2)])
    }
  }

  writexl::write_xlsx(lapply(out_contab_summary, as.data.frame), paste(path, "figures/", "Fig3_data_", clean_filename(MDI_name), ".xlsx", sep=""))

  out_contab_summary_tmp = data.frame(matrix(data = NA, ncol = 3, nrow = length(out_contab)))
  colnames(out_contab_summary_tmp) = c("Lower_CI", "Mean", "Upper_CI")
  rownames(out_contab_summary_tmp) = c("Mathematical sum", "Best peak - global max", "Best peak - local max", "AI scaling - highest fragment",
                                       "AI scaling - heaviest fragment", "p-value integration", "Multivariate testing", "MSstats - linear", "MSstats - TMP")
  out_contab_summary_CI = list("PPV" = out_contab_summary_tmp, "Sensitivity" = out_contab_summary_tmp, "Specificity" = out_contab_summary_tmp,
                               "Accuracy" = out_contab_summary_tmp, "Precision" = out_contab_summary_tmp)

  for (i_list in 1:length(out_contab_summary)) {
    for (i_method in 1:nrow(out_contab_summary[[i_list]])) {
      if (i_list == 1) {
        out_contab_summary_CI[[i_list]][i_method,] = CI(na.omit(1 - as.double(out_contab_summary[[i_list]][i_method,])))
      }
      else {
        out_contab_summary_CI[[i_list]][i_method,] = CI(na.omit(as.double(out_contab_summary[[i_list]][i_method,])))
      }
    }
  }

  out_contab_summary_CI["Sensitivity"] = NULL

  plot_data_radar = data.frame(matrix(nrow = nrow(out_contab_summary[[1]]), ncol = length(out_contab_summary_CI)))
  colnames(plot_data_radar) = names(out_contab_summary_CI)
  rownames(plot_data_radar) = rownames(out_contab_summary_CI[[1]])
  for (i_list in 1:length(out_contab_summary_CI)) {
    for (i_method in 1:nrow(out_contab_summary_CI[[i_list]])) {
      plot_data_radar[i_method, i_list] = out_contab_summary_CI[[i_list]][i_method,2]
    }
  }
  plot_data_radar = rbind(rep(1,length(out_contab_summary_CI)) ,
                          rep(0.88,length(out_contab_summary_CI)),
                          plot_data_radar)
  plot_data_radar = plot_data_radar[c(1,2,8,3,9,4,5,6,7,10,11),]
  colnames(plot_data_radar) = c("PPV", "Specificity  ", "Accuracy", "  Precision")

  # Radarplot for *.tiff
  tiff(paste(path, "figures/radarchart_", MDI_name, ".tiff", sep=""), units="in", width=8, height=8, res=1200)
  colors_border = c(
    "#0072FF",  # p-value integration
    "#00C2A8",  # Mathematical sum
    "#7CFF6B",  # Multivariate testing
    "#FFD166",  # Best peak - global max
    "#F4A261",  # Best peak - local max
    "#9B5DE5",  # AI scaling - highest fragment
    "#C77DFF",  # AI scaling - heaviest fragment
    "#E63946",  # MSstats - linear
    "#B22222"   # MSstats - TMP
  )
  colors_in = alpha(colors_border,0.1)
  radarchart(plot_data_radar, pcol=colors_border, pfcol=colors_in, plwd=1, plty=1, cglcol = "grey", cglty = 1,
             vlcex = 1.26, axistype = 1, caxislabels = c("0.88  ", "", "", "", "1.0   "), axislabcol = "black",
             title = MDI_name)
  # legend(x=0.7, y=1.3, legend = rownames(plot_data_radar[-c(1,2),]), bty = "n", pch=20,
  #        col=colors_border , text.col = "black", cex=0.8, pt.cex=2)
  dev.off()

  # Radarplot for *.emf
  devEMF::emf(paste(path, "figures/radarchart_", MDI_name, ".emf", sep=""), units="in", width=8, height=8)
  colors_border = c(
    "#0072FF",  # p-value integration
    "#00C2A8",  # Mathematical sum
    "#7CFF6B",  # Multivariate testing
    "#FFD166",  # Best peak - global max
    "#F4A261",  # Best peak - local max
    "#9B5DE5",  # AI scaling - highest fragment
    "#C77DFF",  # AI scaling - heaviest fragment
    "#E63946",  # MSstats - linear
    "#B22222"   # MSstats - TMP
  )
  colors_in = alpha(colors_border,0.1)
  radarchart(plot_data_radar, pcol=colors_border, pfcol=colors_in, plwd=1, plty=1, cglcol = "grey", cglty = 1,
             vlcex = 1.26, axistype = 1, caxislabels = c("0.85  ", "", "", "", "1.0   "), axislabcol = "black",
             title = MDI_name)
  # legend(x=0.7, y=1.3, legend = rownames(plot_data_radar[-c(1,2),]), bty = "n", pch=20,
  #        col=colors_border , text.col = "black", cex=0.8, pt.cex=2)
  dev.off()

  #Print legend only
  tiff(paste(path, "figures/radarplot_legend.tiff", sep=""), units="in", width=4, height=4, res=1200)
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend(x=0.08, y=0.85, legend = rownames(plot_data_radar[-c(1,2),]), bty = "n", pch=20,
         col=colors_border , text.col = "black", cex=0.8, pt.cex=2)
  dev.off()
  devEMF::emf(paste(path, "figures/radarplot_legend.emf", sep=""), units="in", width=4, height=4)
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend(x=0.08, y=0.85, legend = rownames(plot_data_radar[-c(1,2),]), bty = "n", pch=20,
         col=colors_border , text.col = "black", cex=0.8, pt.cex=2)
  dev.off()

  # # Legacy full plots
  # for (i_param in 1:length(out_contab_summary)) {
  #   param_name = names(out_contab_summary)[i_param]
  #   out_contab_summary_tmp = out_contab_summary[[i_param]]
  # 
  #   plot_data = data.frame(matrix(data=NA, ncol = 3, nrow = 0))
  #   colnames(plot_data) = c("Method", "Grouping", "Parameter")
  #   for (i_method in 1:nrow(out_contab_summary_tmp)) {
  #     for (i_col in 1:ncol(out_contab_summary_tmp)) {
  #       plot_data[nrow(plot_data)+1, 1] = rownames(out_contab_summary_tmp)[i_method]
  # 
  #       conc = as.numeric(strsplit(colnames(out_contab_summary_tmp)[i_col], "_")[[1]])
  #       conc = format(round(((conc / 31.5) * 5), digits=2), nsmall = 2)
  # 
  #       plot_data[nrow(plot_data), 2] = paste(conc[1], " vs ", conc[2], sep="")
  #       plot_data[nrow(plot_data), 3] = out_contab_summary_tmp[i_method, i_col]
  #     }
  #   }
  # 
  #   all_group_comb_conc_1 = format(round(((as.numeric(all_group_combinations[1,]) / 31.5) * 5), digits=2), nsmall = 2)
  #   all_group_comb_conc_2 = format(round(((as.numeric(all_group_combinations[2,]) / 31.5) * 5), digits=2), nsmall = 2)
  #   plot_data$Grouping = factor(plot_data$Grouping, levels = paste(all_group_comb_conc_1, " vs ", all_group_comb_conc_2, sep=""))
  # 
  #   p = ggplot(data = plot_data, aes(x = Grouping, y = Parameter, group = Method, color = Method)) +
  #     geom_line(size=1) + ylab(param_name) + xlab("Grouping") + ggtitle(paste(param_name, " plot: ", MDI_name, sep="")) +
  #     scale_y_continuous(labels = scales::percent) +
  #     theme(plot.title = element_text(size = 12, face="bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 11),
  #           axis.ticks = element_line(color = "black"), axis.text.y = element_text(color="black", size = 11),
  #           axis.title = element_text(size = 12)) + theme(legend.position = "none")
  # 
  #   ggsave(filename = paste(path, "figures/", param_name, "_", MDI_name, ".tiff", sep=""),
  #          width = 18.2, height = 13.6, units = "cm", dpi = "retina")
  #   ggsave(filename = paste(path, "figures/", param_name, "_", MDI_name, ".emf", sep=""),
  #          width = 18.2, height = 13.6, units = "cm", dpi = "retina")
  # 
  #   if (i_param*i_MDI == 15) {
  #     p = p + theme(legend.position = "right")
  #     p
  #     ggsave(filename = paste(path, "figures/", param_name, "_", MDI_name, "_LEGEND.tiff", sep=""),
  #            width = 22.2, height = 13.6, units = "cm", dpi = "retina")
  #     ggsave(filename = paste(path, "figures/", param_name, "_", MDI_name, "_LEGEND.emf", sep=""),
  #            width = 22.2, height = 13.6, units = "cm", dpi = "retina")
  #   }
  # }
}

######################################################################################



























