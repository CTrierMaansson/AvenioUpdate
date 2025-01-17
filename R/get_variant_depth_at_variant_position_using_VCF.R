#' @noRd
#' @importFrom stringr str_glue
get_variant_depth_at_variant_position_using_VCF<-function(VCF,Reads){
    ct_depth <- as.integer(unlist(mapply(VariantDepth, 
                                         Chr=VCF@fix[,1], 
                                         Pos=as.numeric(VCF@fix[,2]), 
                                         Base=VCF@fix[,5],
                                         MoreArgs=list(Bam=Reads))))
    wt_depth <- as.integer(unlist(mapply(VariantDepth, 
                                         Chr=VCF@fix[,1], 
                                         Pos=as.numeric(VCF@fix[,2]), 
                                         Base=VCF@fix[,4],
                                         MoreArgs=list(Bam=Reads))))
    total_depth<-ct_depth+wt_depth
    chrom<-VCF@fix[,1]
    pos<-as.numeric(VCF@fix[,2])
    position<-string::str_glue("{chrom}:{pos}")
    return(c(position,ct_depth,total_depth))
}