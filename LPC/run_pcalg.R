library("pcalg")
setwd("h:\\projects\\LPC\\LPC")
n_genes<-100;
measurements<-c(10,20,50,100,200,500,1000);
for (i in 1:10){
 for (j in 1:7)
{
    fname<-paste("..\\data\\data_100genes\\100g_10rec_x_",i,"_",measurements[j],"measurements.txt",sep="")
   
    fname1<-paste("..\\results\\pc\\",n_genes,"g_10rec_x_",i,"_",measurements[j],"_pc.txt",sep="")
    fname2<-paste("..\\results\\pc\\",n_genes,"g_10rec_x_",i,"_",measurements[j],"_cpdag.txt",sep="")
 

   t<-read.table(fname,header=FALSE)

   #t<-t(t);
   #x<-pcAlgo(t,alpha=0.005);

   #r<-ggm.test.edges(x);
   #net<-extract.network(r,cutoff.ggm=0.9)
   #write.table(net,file=fname1,sep="\t")

   resD <- pcAlgo(t, alpha = 0.05, corMethod ="standard",directed=FALSE,0)
   resD <- udag2pdag(resD,verbose=1)
  # resD <- pcAlgo(t, alpha = 0.05, corMethod ="standard",directed=FALSE)

   edges<-edges(resD@graph)

   edges$'100'[1]

   plot(resD,zvalue.lwd=TRUE)
   m<-wgtMatrix(resD@graph)
   write.table(resD@zMin,file=fname1,sep="\t",row.names = FALSE,col.names =FALSE)
   write.table(m,file=fname2,sep="\t",row.names = FALSE,col.names =FALSE)

}
} 
