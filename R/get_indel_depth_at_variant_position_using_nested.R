#' @noRd
#' @importFrom stringr str_glue
get_indel_depth_at_variant_position_using_nested<-function(Variants,Reads){
    cigar<-Reads[[1]]$cigar
    indel<-Variants$type
    #find ct og wt deoth:
    ct_depth <- sum(grepl(indel,cigar))
    wt_depth <- sum(!grepl(indel,cigar))
    #beregn de forskellige ting:  
    total_depth<-ct_depth+wt_depth
    chrom<-Variants$chrom
    pos<-as.numeric(Variants$pos)
    position<-stringr::str_glue("{chrom}:{pos}")
    return(c(position,ct_depth,total_depth))
}