# summarization portion of existing 190131 DE analysis pipeline

require(Rsubread)
params = commandArgs(trailingOnly = T) #grabs parameters from bash call
setwd(params[1]) #root directory of dataset
stranding = as.numeric(params[2])
threads = as.numeric(params[3])
gtfPath = params[4]
minFrag = as.numeric(params[5])
includeMM = as.numeric(params[6])
outPrefix = params[7]

#perform summarization with given parameters
files = list.files(path = "./alignment", pattern = "Aligned\\.out\\.sam$", full.names = T)
counts = Rsubread::featureCounts(files = files, annot.ext = gtfPath, minFragLength = minFrag, 
                                 nthreads = threads, isGTFAnnotationFile = T, countMultiMappingReads = includeMM)

#output processed counts
outName = paste0(getwd(), "/counts/", outPrefix, "_subread_counts.rds")
saveRDS(counts, outName)