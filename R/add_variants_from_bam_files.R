#' @noRd
add_variants_from_bam_files<- function(df,synology_path){
    #tjek bamtiler for snv og indels
    snv_fundet_i_bam<-wrapper_tjek_bam_snv(df,synology_path)
    indel_findet_i_bam<-wrapper_tjek_bam_indel(df,synology_path)
    
    #tilføj varianter til df:
    if(nrow(snv_fundet_i_bam)>0){
        df<-wrapper_add_variants(df,snv_fundet_i_bam)}
    
    if(nrow(indel_findet_i_bam)>0){
        df<-wrapper_add_variants(df,indel_findet_i_bam)}
    
    return(df)
}