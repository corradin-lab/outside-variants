Args <- commandArgs()
heatmap_output<-read.delim(Args[3], head=F)
values<-heatmap_output[c(-1,-2)]
chr_coord<-paste(heatmap_output$V1,heatmap_output$V2, sep='_')
values<-apply(values,2,as.numeric)
mean<-mean(values, na.rm=T)
std<-sd(as.vector(values), na.rm=T)
zscores<-(values-mean)/std
zscoresfinal<-cbind(chr_coord, zscores)
write.table(zscoresfinal, "temp.txt",col.names=F, row.names=F, sep="\t",quote=F)
q()