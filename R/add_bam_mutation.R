#' @noRd
#' @importFrom dplyr filter `%>%` any_of select left_join bind_cols
#' @importFrom stringr str_glue
add_bam_mutation<-function(df_mut,Sample1,Sample2,GenomicPosition,
                           VariantDepth,UniqueDepth){
    #lav række med informationer fra sample 1 (informationer fra kørsel):
    df_sample1<-df_mut%>%dplyr::filter(sample_index==Sample1)
    list_sample1<-df_sample1[1,]
    #lav række med informationer fra sample 2 (informationer fra varianten):
    list_sample2<-df_mut%>%
        dplyr::filter(sample_index==Sample2&Genomic.Position==GenomicPosition)
    
    
    redundant_cols <- redundant_columns()
    list_sample1_red <- list_sample1 %>% 
        dplyr::select(dplyr::any_of(redundant_cols))
    list_sample2 <- list_sample2 %>% 
        dplyr::select(-dplyr::any_of(redundant_cols))
    merged_list <- bind_cols(list_sample1_red,list_sample2) 
    merged_list$Variant.Depth <- as.character(VariantDepth)
    merged_list$Unique.Depth <- as.character(UniqueDepth)
    merged_list$Allele.Fraction <- stringr::str_glue(
        "{round((VariantDepth/UniqueDepth)*100,digits=2)}%")
    merged_list$No..Mutant.Molecules.per.mL <- as.character(
        round(
            (as.numeric(
                merged_list$Isolated.DNA.Mass..ng.)*330*(VariantDepth/UniqueDepth))
            /as.numeric(merged_list$Plasma.Volume..mL.),2))
    merged_list$Flags<-"BAM"
    #Listen der bestemmer,¨
    #hvilke index der skal komme fra hvilken liste (0 = list1, 1 = list2)
    # index_list<-c(rep(0,15),rep(1,31),rep(0,39))
    # # Funktion til at vælge elementer baseret på index_list
    # merge_lists <- function(list1, list2, index_list) {
    #     mapply(function(x, y, z) if (z == 0) x else y, list1, list2, index_list)
    # }
    # # Brug af funktionen
    # merged_list <- merge_lists(list_sample1, list_sample2, index_list)
    # 
    # #tilføj dybder fra BAM fil:
    # merged_list[25]<-as.character(VariantDepth)
    # merged_list[26]<-as.character(UniqueDepth)
    # #Beregn AF og MM:
    # merged_list[22]<-stringr::str_glue(
    #     "{round((VariantDepth/UniqueDepth)*100,digits=2)}%")
    # merged_list[23]<-as.character(
    #     round(
    #         (as.numeric(
    #             merged_list[9])*330*(VariantDepth/UniqueDepth))
    #         /as.numeric(merged_list[10]),2))
    
    # indsæt"BAM" under flags:
    # merged_list[3]<-"BAM"
    
    #Hvis der ikke er fundet en mutation i sample 1 tidligere,
    #så udskiftes den gamle række med den nye:
    if (is.na(list_sample1$Mutation.Class)){
        df_mut<-df_mut%>%dplyr::filter(sample_index!=list_sample1[[1]])
    }
    merged_list <- merged_list %>% 
        unique()
    #tilføj den nye række:
    df_mut<-rbind(df_mut,merged_list)
    return(df_mut)
}