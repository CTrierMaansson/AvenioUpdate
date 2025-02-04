#' @noRd
#' @importFrom utils glob2rx
#' @importFrom Rsamtools BamFile ScanBamParam scanBam scanBamWhat
#' @importFrom stringr str_glue
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
get_reads_at_variant_position_using_nested_in_posdeduped <- function(Variants, SampleID , Path) {
    #find navn på PosDeduped bam fil:
    bamfile <- list.files(Path,pattern=utils::glob2rx("PosDeduped-*.bam"))
    bamfileP <- bamfile[grep(SampleID, bamfile)]
    #find navn på PosDeduped bai fil:
    indexx <- list.files(Path,pattern=utils::glob2rx("PosDeduped-*.bai"))
    indexx <- indexx[grep(SampleID, indexx)]
    #tjek om der er fundet to filnavne:
    if (length(bamfileP) + length(indexx) != 2) 
    {message(paste0("BAM file ",
                  SampleID,
                  " not found, proceeding without."))}
    else {
        #lav reference til bamfilen og dens index med BamFile():
        bamfile <- Rsamtools::BamFile(
            stringr::str_glue("{Path}/{bamfileP}"), yieldSize=5e8)
    
        message(paste0("Processing ", SampleID))
        #VCF$fix: Column 1 is the chromosome number, 
        #2 is the chromosome position
        #Lav parameter object med positionerne for variaterne:
        param <- Rsamtools::ScanBamParam(
            which=GenomicRanges::GRanges(
                Variants$chrom,
                IRanges::IRanges(start=as.numeric(Variants$pos),
                                 end=as.numeric(Variants$pos))), 
            what=Rsamtools::scanBamWhat())
        #Importer info fra alle reads der dækker positionerne i param
        bam <- Rsamtools::scanBam(bamfile, param=param)
        return(bam)
    }
}