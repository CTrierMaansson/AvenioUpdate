#' @noRd
#' @importFrom utils glob2rx
#' @importFrom Rsamtools ScanBamParam scanBam scanBamWhat BamFile
#' @importFrom stringr str_glue
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
get_reads_at_variant_position_using_VCF_in_deduped <- function(VCF, SampleID , Path) {
    #tjek om der er fundet en variant i prøven:
    if (nrow(VCF) == 0) {print(paste0("No variants in ", SampleID))}
    else {
        #find navn på Deduped bam fil:
        bamfile <- list.files(Path,pattern=utils::glob2rx("Deduped-*.bam"))
        bamfileP <- bamfile[grep(SampleID, bamfile)]
        #find navn på Deduped bai fil:
        indexx <- list.files(Path,pattern=utils::glob2rx("Deduped-*.bai"))
        indexx <- indexx[grep(SampleID, indexx)]
        #tjek om der er fundet to filnavne:
        if (length(bamfileP) + length(indexx) != 2) 
        {print(paste0("BAM file ", 
                      SampleID, 
                      " not found, proceeding without."))}
        else {
            #lav reference til bamfilen og dens index med BamFile():
            bamfile <- Rsamtools::BamFile(
                stringr::str_glue("{Path}/{bamfileP}"), yieldSize=5e8)
            
            print(paste0("Processing ", SampleID))
            #VCF$fix: Column 1 is the chromosome number, 
            #2 is the chromosome position
            #Lav parameter object med positionerne for variaterne:
            param <- Rsamtools::ScanBamParam(
                which=GenomicRanges::GRanges(
                    VCF@fix[,1],
                    IRanges::IRanges(start=as.numeric(VCF@fix[,2]),
                                     end=as.numeric(VCF@fix[,2]))),
                                  what=Rsamtools::scanBamWhat())
            #Importer info fra alle reads der dækker positionerne i param
            bam <- Rsamtools::scanBam(bamfile, param=param)
            
            return(bam)
        }
    }
}