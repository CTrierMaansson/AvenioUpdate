
zip_path1 <- "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AbpTmAnr_h9DPL1skMvyje5U"
zip_path2 <- "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ABE9wiNEoilBeYhlV1UO3VRC"
zip_path3 <- "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AMw9wjlnV4lHRanHC8-jUnoh"
test_path <- "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ANei_Hsi__lGmaJ1k_Ji_d3O"



library(readxl)
library(dplyr)
library(stringr)
`%ni%` <- Negate(`%in%`)

add_samples <- function(Directory){
  Run_ID <- gsub("//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-","",Directory)
  zip_file <- paste0(Directory,"/","all_CSVs-",Run_ID,".zip")
  sample_file <- paste0(Directory,"/SampleMetrics.csv")
  run_ID_short <- substr(Run_ID,1,8)
  unzip(zipfile = zip_file,
        exdir = Directory)
  Filter_variants <- paste0(Directory,"/","VariantMetrics-UserFilterSet.csv")
  variants <- read.csv(Filter_variants)
  AVENIO_runs <- read_xlsx("//Synology_m1/Synology_folder/AVENIO/AVENIO_runs.xlsx",col_types = c(rep("guess",4),"date",rep("guess",6)))
  AVENIO_runs_select <- AVENIO_runs %>% 
    filter(grepl(run_ID_short,Run_ID))
  variants_select <- variants %>% 
    filter(Sample.ID %in% AVENIO_runs_select$Sample_name)
  samples_with_variants <- unique(variants_select$Sample.ID)
  samples_of_interest <- unique(AVENIO_runs_select$Sample_name)
  if(length(samples_with_variants) != length(samples_of_interest)){
    sample_metrics <- read.csv(sample_file) %>% 
      filter(Sample.ID %in% samples_of_interest) %>% 
      filter(Sample.ID %ni% samples_with_variants)
    sample_metrics <- sample_metrics %>% 
      select(colnames(sample_metrics)[colnames(sample_metrics) %in% colnames(variants_select)])
    sele_cols <- colnames(variants_select)[colnames(variants_select) %in% colnames(sample_metrics)]
    for(i in 1:length(sele_cols)){
      if(class(sample_metrics[[sele_cols[i]]]) != class(variants_select[[sele_cols[i]]])){
        classes <- c(class(sample_metrics[[sele_cols[i]]]),class(variants_select[[sele_cols[i]]]))
        if("character" %in% classes){
          sample_metrics[[sele_cols[i]]] <- as.character(sample_metrics[[sele_cols[i]]])
          variants_select[[sele_cols[i]]] <- as.character(variants_select[[sele_cols[i]]])
        }
        if("numeric" %in% classes){
          sample_metrics[[sele_cols[i]]] <- as.numeric(sample_metrics[[sele_cols[i]]])
          variants_select[[sele_cols[i]]] <- as.numeric(variants_select[[sele_cols[i]]])
        }
      }
    }
    combined_samples <- full_join(variants_select,sample_metrics,copy = TRUE)
  }
  else{
    combined_samples <- variants_select
  }
  if(!"sample_index" %in% colnames(AVENIO_runs_select)){
    AVENIO_runs_select <- AVENIO_runs_select %>% 
      mutate(sample_index = paste0(Project,"_",
                                   Name_in_project,"_",
                                   substr(str_split_i(as.character(Sample_date),"-",1),3,4),
                                   str_split_i(as.character(Sample_date),"-",2),
                                   str_split_i(as.character(Sample_date),"-",3))) %>% 
      mutate(sample_index = ifelse(Material != "cfDNA",paste0(sample_index,"_",Material),sample_index))
  }
  AVENIO_runs_select_merge <- AVENIO_runs_select %>% 
    dplyr::select(Sample_name,sample_index)
  combined_samples <- combined_samples %>% 
    left_join(AVENIO_runs_select_merge,by = c("Sample.ID"="Sample_name")) %>% 
    relocate(sample_index,.before = Sample.ID) %>% 
    arrange(sample_index,Gene)
  return(list(combined_samples,
              AVENIO_runs_select))
}
sample_res_1 <- add_samples(zip_path1)
sample_res_2 <- add_samples(zip_path2)
sample_res_3 <- add_samples(zip_path3)
combined_samples1 <- AvenioUpdate:::add_DNAfusion_res(sample_res_1[[1]],
                                                      repo_path = zip_path1,
                                                      sample_info = sample_res_1[[2]])
combined_samples2 <- AvenioUpdate:::add_DNAfusion_res(sample_res_2[[1]],
                                                      repo_path = zip_path2,
                                                      sample_info = sample_res_2[[2]])
combined_samples3 <- AvenioUpdate:::add_DNAfusion_res(sample_res_3[[1]],
                                                      repo_path = zip_path3,
                                                      sample_info = sample_res_3[[2]])
create_df_list <- function(sample_df){
  samples <- unique(sample_df$Sample.ID)
  AVENIO_runs <- read_xlsx("//Synology_m1/Synology_folder/AVENIO/AVENIO_runs.xlsx",col_types = c(rep("guess",4),"date",rep("guess",6)))
  AVENIO_runs_select <- AVENIO_runs %>% 
    filter(Sample_name %in% samples)
  unique_patients <- unique(AVENIO_runs_select$CPR)
  df_list <- list()
  for (i in 1:length(unique_patients)){
    CPR_patient <- unique_patients[i]
    patient_df <- AVENIO_runs_select %>% 
      filter(CPR == CPR_patient)
    patient_variants <- sample_df %>% 
      filter(Sample.ID %in% patient_df$Sample_name)
    df_list[[i]] <- patient_variants 
    names(df_list)[i] <- CPR_patient
  }
  return(df_list)
}
df_list1 <- create_df_list(combined_samples1)
df_list2 <- create_df_list(combined_samples2)
df_list3 <- create_df_list(combined_samples3)
master_list <- append(df_list1,df_list2)
length(master_list)
reanalyze_samples <- function(master_list, df_list){
  existing <- names(df_list)[names(df_list) %in% names(master_list)]
  master_sele <- master_list[names(master_list) %in% existing]
  df_list_sele <- df_list[names(df_list) %in% existing]
  test_df <- do.call(rbind,append(master_sele,df_list_sele)) %>% 
    arrange(sample_index,Gene) %>% 
    unique()
  rownames(test_df) <- NULL
  validated_samples <- add_variants_from_bam_files(test_df)
  validated_list <- create_df_list(validated_samples)
  for(i in 1:length(validated_list)){
    sample_df <- validated_list[[i]]
    if(any(grepl("_BC",sample_df$sample_index))){
      patient_sample <- paste0(str_split_i(unique(sample_df$sample_index)[1],"_",1),"_",str_split_i(unique(sample_df$sample_index)[1],"_",2))
      print(paste0("There is a BC sample for ", patient_sample, " and variants are marked accordingly"))
      sample_BC <- sample_df %>% 
        filter(grepl("_BC", sample_index))
      sample_other <- sample_df %>% 
        filter(!grepl("_BC", sample_index))
      if(nrow(sample_other) > 0){
        sample_BC <- sample_BC %>% 
          mutate(identifier = paste0(Gene,"_",Coding.Change,"_",Genomic.Position))
        sample_other <- sample_other %>% 
          mutate(identifier = paste0(Gene,"_",Coding.Change,"_",Genomic.Position)) %>% 
          mutate(Flags = ifelse(identifier %in% sample_BC$identifier,
                                ifelse(nchar(Flags)>0,
                                       paste0(Flags,", BC_mut"),
                                       paste0("BC_mut")),
                                Flags)) %>% 
          select(-identifier)
        sample_BC <- sample_BC %>% 
          select(-identifier)
        sample_df <- rbind(sample_other,sample_BC) %>% 
          arrange(sample_index,Gene)
      }
    }
    validated_list[[i]] <- sample_df
  }
  master_list <- master_list[names(master_list) %ni% names(validated_list)]
  master_list <- append(master_list,validated_list)
  df_list_new <- df_list[names(df_list) %ni% existing]
  master_list <- append(master_list,df_list_new)
  master_list <- master_list[sort(names(master_list))]
  return(master_list)
}

master_list <- reanalyze_samples(master_list = master_list, df_list = df_list3)

add_run_to_list <- function(master_list, Directory){
  samples <- add_samples(Directory)
  n_patients_before <- length(master_list)
  unlisted_before <- do.call(rbind,master_list) %>% 
      select(Analysis.Name,Sample.ID) %>% 
      unique()
  n_runs_before <- nrow(unlisted_before)
  df_list <- create_df_list(samples)
  reanalyzed <- reanalyze_samples(master_list,df_list)
  if(is.list(reanalyzed)){
      if(length(reanalyzed) >= length(master_list)){
          n_patients_after <- length(reanalyzed)
          unlisted_after <- do.call(rbind,reanalyzed) %>% 
              select(Analysis.Name,Sample.ID) %>% 
              unique()
          n_runs_after <- nrow(unlisted_after)
          print("Before the dataset consisted of ",
                n_patients_before, 
                " patients and ",
                n_runs_before, 
                " samples analyzed")
          print("Now the dataset consists of ",
                n_patients_after, 
                " patients and ",
                n_runs_after, 
                " samples analyzed")
          saveRDS(reanalyzed,file = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
      }
  }
  return(reanalyzed)
}


for(i in 1:length(master_list)){
    print(class(master_list[[i]]))
}

master_list <- add_run_to_list(master_list = master_list,Directory = test_path)


create_simple_output <- function(df_list,CPR,synonymous = TRUE){
  df <- df_list[[CPR]]
  if(is.null(df)){
    stop("CPR number is not in the dataset")
  }
  print("Reading AVENIO_runs.xlsx file")
  AVENIO_runs <- read_xlsx("//Synology_m1/Synology_folder/AVENIO/AVENIO_runs.xlsx",
                           col_types = c(rep("guess",4),"date",rep("guess",6)))
  print("Selecting relevant sample")
  AVENIO_sele <- AVENIO_runs[AVENIO_runs$CPR == CPR,]
  AVENIO_sele <- AVENIO_sele[!is.na(AVENIO_sele$CPR),] %>% 
    select(Sample_name, Sample_note, Material, Notes)
  print("Marging run and variant information")
  df <- df %>% 
    left_join(AVENIO_sele, by = c("Sample.ID" = "Sample_name")) %>% 
    select(sample_index,
           Mutation.Class,Gene,Amino.Acid.Change,
           Variant.Description,Flags,Allele.Fraction,
           Variant.Depth,Unique.Depth,Analysis.Name,Sample.ID,
           Sample_note, Material, Notes) %>%
    mutate(Variant.Description = gsub(" variant","",Variant.Description)) %>% 
    tibble()
  colnames(df) <- c("sample_index","Class","Gene","AA","Description",
                    "Flags","MAF","Variant_depth","Unique_depth","Analysis",
                    "Sample.ID","Sample_note","Material","Notes")
  if(synonymous){
    df <- df
  }
  else{
      print("Removing Synonymous variants")
    df <- df %>% 
      filter(!grepl("Synonymous",Description))
  }
  return(list("Patient CPR:" = CPR,"Results" = df))
}
create_simple_output(master_list,"2603530058",synonymous = FALSE)
length(master_list)


#### Creating package ####
setwd("C:/Users/chris/OneDrive/1PhD/AvenioUpdate")
devtools::document()
devtools::build()

devtools::install_github("CTrierMaansson/AvenioUpdate",force = T)

AvenioUpdate::explore_AVENIO_runs()
AvenioUpdate::explore_AVENIO_runs_Info()
AvenioUpdate::result_stats()

master_list <- AvenioUpdate::add_run_to_list(master_list = master_list,
                              Directory = test_path)
names(master_list)
AvenioUpdate::create_simple_output(master_list,"0406531473",synonymous = FALSE)

names(analyses[[2]])
AvenioUpdate:::wrapper_tjek_bam_snv

