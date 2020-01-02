
args<-commandArgs(TRUE)


input<-read.delim(args[1],head=T)

aaa <- function(mydata) {
	ID<-apply(mydata[,1:3],2,as.character)

	#rowid <- mydata[,1]
	#geno <-  mydata[,2]
	nr <- nrow(mydata)
	
	final.out <- matrix(nrow=0,ncol=11)
	i=1
	while (i < nr ) {
		temp.out <- matrix(nrow=0, ncol=0)
		temp.out <- matrix(nrow=1, ncol=11)
		frst_A <- apply(mydata[i,4:ncol(mydata[i,])],2,as.numeric)
		scnd_A <- apply(mydata[i+1,4:ncol(mydata[i+1,])],2,as.numeric)
		#frst_M <- apply(mydata[i,12:35],2,as.numeric)
		#scnd_M <- apply(mydata[i+1,12:35],2,as.numeric)
		
		temp.out[1,1] <- ID[i,2]	
		temp.out[1,2] <- ID[i+1,2]
		temp.out[1,3] <- ID[i,3]
		temp.out[1,4] <- ID[i+1,3]
		
		temp.out[1,5] <- mean(frst_A,na.rm=TRUE)
		temp.out[1,6] <- mean(scnd_A,na.rm=TRUE)
		temp.out[1,7] <- sd(frst_A,na.rm=TRUE)
		temp.out[1,8] <- sd(scnd_A,na.rm=TRUE)
		
		temp.out[1,9] <- tryCatch(ks.test(frst_A,scnd_A,alternative = c("less"),exact=NULL)$p.value, error = function(e) rep (NA, 1)); 
		temp.out[1,10] <- ID[i,1]

		final.out <- rbind(final.out, temp.out)
		i = i+2
	}
	
	return(final.out)
}

nn<-paste(args[2],"aggregate_ksonesided_paired_greater_out",sep="_")
input_out<-aaa(input)
write.table(input_out,file=nn, sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
