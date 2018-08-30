# Exactly like the homologene package, but uses most recent NCBI homologene dataset
tempLoc = tempfile()
raw = download.file("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data", tempLoc)
hg = read.delim(tempLoc, header = F)
colnames(hg) = c("HID", "taxid", "geneID", "symbol", "proteinID", "proteinAccession")
hg$symbol = as.character(hg$symbol) #remove factorization

hgConvert = function(geneSymbols, inTax, outTax) #analogous to hgConvert()
{
  df_in = subset(hg, taxid == inTax)
  df_out = subset(hg, taxid == outTax)
  outData = data.frame(matrix(NA, nrow = length(geneSymbols), ncol = 4))
  colnames(outData) = c(inTax, outTax, paste0(inTax, "_ID"), paste0(outTax, "_ID"))
  
  for(i in 1:length(geneSymbols))
  {
    outData[i, 1] = geneSymbols[i] #inTax gene symbol
    inSet = subset(df_in, symbol == geneSymbols[i])
    curHID = inSet$HID #same between both taxa
    inGeneID = inSet$geneID
    outData[i, 3] = ifelse(length(inSet$geneID) == 1, inGeneID, NA) #inTax gene ID, make sure it's in the dataset first
    
    outSet = subset(df_out, HID == curHID)
    outSymbol = outSet$symbol[1] #hoping that this normally just returns one value
    outData[i, 2] = ifelse(length(outSymbol) == 1, outSymbol, NA) #outTax gene symbol
    outGeneID = outSet$geneID[1]
    outData[i, 4] = ifelse(length(outGeneID) == 1, outGeneID, NA) #outTax geneID
  }
  return(na.omit(outData))
}