require(ape)
require(dplyr)

#require muscle aligment to be locally installed.

Ciclando <- function(ind1,ind2,ind3,ind4) { 
  
  data_list<- c(ind1,ind2,ind3,ind4)
  data_cicles<- c(35,40,45,50,55)
  save<- numeric()
  
  for (i in 1:length(data_list)){
    individual<- data_list[i]
    cons_name1<- paste(individual,"consensus", "R1.fasta", sep="-")
    cons_name2<- paste(individual,"consensus", "R2.fasta", sep="-")
    
    cons1<- read.dna(file=cons_name1, format="fasta")
    cons2<- read.dna(file=cons_name2, format="fasta")
    
    for (i in 1:length(data_cicles)){
      
      data_name1<- paste(individual, data_cicles[i], "R1 .fasta", sep="-")
      print(data_name1)
      data_name2<- paste(individual, data_cicles[i], "R2 .fasta", sep="-")
      data1<-read.dna(data_name1, format= "fasta")
      data2<-read.dna(data_name2, format= "fasta")
      data_cicles_def<-data_cicles[i] 
      
      for (i in 1:nrow(data1)){
        print(i)
        print(data_name1)
        sequence<- data1[i,]
        if(ncol(sequence)==ncol(cons1)){
          
          matrix<- rbind(cons1,sequence)
          dist<-dist.dna(matrix, model="raw")[1]
          dist<- dist * 50
          save<- rbind(save,c(i,individual, data_cicles_def, dist,"R1"))
          
        } else {
          
          sequence<-gsub("-", "1", sequence[1,])
          sequence<-sequence[, sequence != "1"]
          
          if(length(sequence)==ncol(cons1)){
            sequence<- as.DNAbin(sequence)
            matrix<- rbind(cons1,sequence)
            dist<-(dist.dna(matrix, model="raw")[1])+0.02
            dist<- as.integer(dist * 50)
            save<- rbind(save,c(i,individual, data_cicles_def, dist,"R1_gap"))
            
          } else {
            sequence<- data1[i,]
            matrix<- append(as.list(cons1),as.list(sequence))
            matrix<- muscle(matrix)
            dist<-(dist.dna(matrix, model="raw")[1])+0.02
            dist<-as.integer(dist * (length(matrix)/2))
            save<- rbind(save,c(i,individual, data_cicles_def, dist,"R1_indel"))
          }
        }
      }  
      for (i in 1:nrow(data2)){
          sequence<- data2[i,]
          if(ncol(sequence)==ncol(cons2)){
            matrix<- rbind(cons2,sequence)
            dist<-dist.dna(matrix, model="raw")[1]
            dist<- dist * 50
            save<- rbind(save,c(i,individual, data_cicles_def, dist,"R2"))
            
          } else {
            sequence<-gsub("-", "1", sequence[1,])
            sequence<-sequence[, sequence != "1"]
            if(length(sequence)==ncol(cons1)){
              sequence<- as.DNAbin(sequence)
              
              matrix<- rbind(cons2,sequence)
              dist<-(dist.dna(matrix, model="raw")[1])+0.02
              dist<- as.integer(dist * 50)
              save<- rbind(save,c(i,individual, data_cicles_def, dist,"R2_gap"))
              
            } else {
              sequence<- data1[i,]
              matrix<- append(as.list(cons2),as.list(sequence))
              matrix<- muscle(matrix)
              dist<-(dist.dna(matrix, model="raw")[1])+0.02
              dist<-as.integer( dist * (length(matrix)/2))
              save<- rbind(save,c(i,individual, data_cicles_def, dist,"R2_indel"))
            }
            
            
          }
          
        }
        
      }
      
    }
    write.csv(save, "data_ouput_mut.csv")
  }  
    
#Ciclando("CEP006","CEP007", "CEP016","CEP023") Use example.
