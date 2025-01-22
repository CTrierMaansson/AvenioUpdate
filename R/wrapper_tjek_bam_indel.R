#' @noRd
#' @importFrom dplyr tibble `%>%` filter
#' @importFrom stringr str_glue
wrapper_tjek_bam_indel<-function(df){
    
    df_varianter_fundet_i_bam<-dplyr::tibble(
        sample_1=c(0),
        sample_2=c(0),
        genomic_position=c(0),
        variant_depth=c(0),
        unique_depth=c(0))
    df_nested<-lav_nested_df_indel(df)
    path_AVENIO_results<-"//Synology_m1/Synology_folder/AVENIO/AVENIO_results"
    
    for(i in 1:nrow(df_nested)){
        #get sample1 name:
        sample1<-df_nested$Sample[i]
        print(paste0("Processing ", sample1))
        
        #path til Bam fil sample i:
        analyseID<-df_nested$Analysis.ID[i]
        sampleID<-df_nested$Sample.ID[i]
        
        files <- list.files(stringr::str_glue(
          "{path_AVENIO_results}/Plasma-{analyseID}/"),
          pattern = stringr::str_glue("{sampleID}-"),
          full.names = TRUE)
        path_sample_mappe <- files[grepl("-Surveillance","-Expanded"),files]
        
        
        #Andre prøver fra samme patient:
        df_patient_other_samples<-df_nested%>%
            dplyr::filter(Patient==df_nested$Patient[i]&Sample!=sample1)
        
      #hvis der er kørt flere prøver på patienten, så tjekkes disse en efter en:
        if (nrow(df_patient_other_samples) == 0)
            {print(paste0("No other samples for the patient than ", sample1))}
        else{
            #loop over andre prøver:
            for (j in 1:nrow(df_patient_other_samples)){
                sample2<-df_patient_other_samples$Sample[j]
                varianter<-df_patient_other_samples$variants[j][[1]]
                #Hvis indel fundet, tjekkes disse i BAM:
                if (nrow(varianter)>0){print(paste0("INDEL in ", sample2))
                    
                    for (v in 1:nrow(varianter)){
                        variant<-varianter[v,]
                        
                        #tjek om variant allerede er fundet i sample1:
                        x<-df_nested$variants[[i]][,1:2]
                        chrom<-variant$chrom
                        pos<-variant$pos
                        
                        if (sum(chrom==x$chrom&pos==x$pos)==0){
                            
    reads_variant_pos<-get_reads_at_variant_position_using_nested_in_posdeduped(
                  variant,sampleID,path_sample_mappe)
    variant_depth_at_pos<-get_indel_depth_at_variant_position_using_nested(
        variant,reads_variant_pos)
    liste_variant<-c(c(sample1,sample2),variant_depth_at_pos)
    #hvis varianten er fundet med mininum 1 read, tilføjes denne:
    if (liste_variant[4]!="0"){
          df_varianter_fundet_i_bam<-rbind(
              df_varianter_fundet_i_bam,liste_variant)
                            }}
                    }}
            }}
    }
    #slet den tomme række
    df_varianter_fundet_i_bam<-df_varianter_fundet_i_bam[-1,]
    
    return(df_varianter_fundet_i_bam)
}