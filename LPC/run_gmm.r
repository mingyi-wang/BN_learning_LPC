library(corpcor);

library("GeneNet");
setwd("h:\\projects\\LPC\\simulation");
measurements<-c(10,20,50,100,200,500,1000);
for (i in 1:10){
 for (j in 1:7)
{
   fname<-paste("..\\data\\data_100genes\\100g_10rec_x_",i,"_",measurements[j],"measurements.txt",sep="")

   fname1<-paste("..\\results\\ggm\\100g_10rec_x_",i,"_",measurements[j],"_ggm.txt",sep="")
   fname2<-paste("..\\results\\ggm\\100g_10rec_x_",i,"_",measurements[j],"_ggm_net.txt",sep="")
   fname3<-paste("..\\results\\ggm\\100g_10rec_x_",i,"_",measurements[j],"_ggm_edg.txt",sep="")

   t<-read.table(fname,header=FALSE)

   
   #x<-ggm.estimate.pcor(t);

     
   
   # below is using shrinkage method to estimate partial correlation
   pc<-ggm.estimate.pcor(t);
   t1.edges <- network.test.edges(pc)
   

   t1.net <- extract.network(t1.edges, cutoff.ggm=0.8)  #use local fdr cutoff 0.2
   
   write.table(t1.net,file=fname2,sep="\t",row.names = FALSE,col.names =FALSE)
   write.table(pc,file=fname1,sep="\t",row.names = FALSE,col.names =FALSE)
   write.table(t1.edges,file=fname3,sep="\t",row.names = FALSE,col.names =FALSE)

}
}
