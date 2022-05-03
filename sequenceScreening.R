library(Biostrings)
library(app)
DNAset = readBStringSet(file.choose(), format="fasta")
dnaLength = length(DNAset[[1]])
genLength = dnaLength - 18
genSet = (successiveViews(DNAset[[1]], width = rep(19, genLength), gapwidth = -18, from = 1))
sink("result.txt") 
cnt = 1
repeat {
   if(cnt > genLength) {
      break
   }
   DNA.str = genSet[[cnt]]
   if(!(countPattern("AAAA", DNA.str) || countPattern("TTTT", DNA.str) || countPattern("CCCC", DNA.str) || countPattern("GGGG", DNA.str))){
      allGC = letterFrequency(DNA.str, "GC", as.prob = TRUE) 
      headGC = letterFrequencyInSlidingView(head(DNA.str, 6), view.width = 6, letter = 'GC')
      tailGC = letterFrequencyInSlidingView(tail(DNA.str, 6), view.width = 6, letter = 'GC')	
      if(allGC >= 0.4 && allGC <= 0.65 && headGC > tailGC ) { 
         site = cnt - dnaLength - 1
         print(site)
         print(DNA.str)
      }
   } 
   cnt <- cnt+1
}
sink()


