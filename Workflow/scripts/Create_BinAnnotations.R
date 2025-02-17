###### https://github.com/asntech/QDNAseq.hg38/blob/main/README.md #######
##### Genmap: https://github.com/cpockrandt/genmap #####
###### bigWigAverageOverBed: https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ #######
###### chrom.sizes: https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/
##### Blacklist: https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/ ######


library(Biobase)
library(BSgenome.Hsapiens.UCSC.hg38)
library(QDNAseq)
library(future)
resourcePath = "/media/sunny/WorkDisk/SV/HCsig/resources/hg38/"
#set virtual mem
options(future.globals.maxSize= 8912896000)
plan(multisession,workers=35)
# change the current plan to access parallelization
#future::plan("multiprocess", workers = 4)

#for (binsize in c(1000, 500, 30, 15, 50, 10, 5, 1)) {
for (binsize in c(30, 50, 100, 150, 200, 250, 500, 1000)) {
  
  bins <- createBins(bsgenome=BSgenome.Hsapiens.UCSC.hg38, binSize=binsize)
  bins$mappability <- calculateMappability(bins,
                                           bigWigFile=paste0(resourcePath,"mappability.hg38.genmap.100mer.bigwig"),
                                           bigWigAverageOverBed="/media/sunny/WorkDisk/SV/HCsig/workflow/scripts/bigWigAverageOverBed")
  
  bins$blacklist <- calculateBlacklist(bins, bedFiles=c(paste0(resourcePath,"wgEncodeDukeMapabilityRegionsExcludable.bed"), paste0(resourcePath,"wgEncodeDacMapabilityConsensusExcludable.bed")))
  
  bins$residual <- NA
  bins$use <- bins$chromosome %in% as.character(1:22) & bins$bases > 0
  
  #
  tg <- binReadCounts(bins,
                      path="/media/sunny/WorkDisk/SV/HCsig/results/02_alignment", cache=TRUE)
  
  bins$residual <- iterateResiduals(tg)
  
  bins <- AnnotatedDataFrame(bins,
                             varMetadata=data.frame(labelDescription=c(
                               "Chromosome name",
                               "Base pair start position",
                               "Base pair end position",
                               "Percentage of non-N nucleotides (of full bin size)",
                               "Percentage of C and G nucleotides (of non-N nucleotides)",
                               "Average mappability of 100mers with a maximum of 2 mismatches",
                               "Percent overlap with ENCODE blacklisted regions",
                               "Median loess residual from 1000 Genomes (100mers)",
                               "Whether the bin should be used in subsequent analysis steps"),
                               row.names=colnames(bins)))
  
  QDNAseqInfo <- list(
    author="Srinivas Veerla",
    date=Sys.time(),
    organism='Hsapiens',
    build='hg38',
    version=packageVersion("QDNAseq"),
    url=paste0(
      "https://github.com/asntech/QDNAseq.hg38/raw/master/data/hg38.",
      binsize, "kbp.SR100.rds"),
    md5=digest::digest(bins@data),
    sessionInfo=sessionInfo())
  
  attr(bins, "QDNAseq") <- QDNAseqInfo
  saveRDS(bins, file=paste0(resourcePath,"BinAnnotations/","hg38.", binsize, "kbp.SR100.rds"))
}

