###########################Threshold Spp influence###########################

#Check sampling effort effect in threshold value determination. 


#Require

library(spider)
library(ggplot2)

#Preparing data
#require A matrix of DNA sequences as csv file:
##Very important## csvfile= Must have the first column with numeric values for species, one number per specie.
dat <- read.dna("name.fas", format="fasta")    
dat<-  as.character(dat)
write.csv(dat, file = "MyData.csv")

#Threshold spp influence

#nspp= Number of species in the matrix
#csvfile= "Name of" csv file 
#it= Number of iterations

##At least requires more than 5 species for threshold calculations


species <- function(csvfile, it) {
  x1=read.csv(csvfile, header= TRUE, sep=",", comment.char="")
  j= as.matrix(x1)
  threshold=numeric()
  index=numeric()
      for (iteration in 1:it){
      c= j[,1]
      A= sample(6:65, 1)
      x <- sample(1:65,A, replace =FALSE)
      r= subset(j, c %in% c(x))
      h= r[,-1]
      p=as.DNAbin(h)
      featherDist <- dist.dna(p, pairwise.deletion = TRUE)
      capture.output(p<-localMinima(featherDist))
      threshold=c(threshold, p$localMinima[1] *100) 
      index=c(index, A)
      indthresholds<- cbind(index, threshold)
      o<- as.matrix(indthresholds)
      output= o
      write.csv(output, file = "output.csv")
  
}
  p=qplot(index, threshold, geom=c("point", "smooth"), method="loess")
  p + theme(axis.title=element_text(face="bold.italic",size="12", color="brown"), legend.position="top")
  p + theme_classic()
}
