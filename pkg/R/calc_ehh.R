calc_ehh<-function(haplohh,mrk,limhaplo=2,limehh=0.05,plotehh=TRUE,main_leg="EHH plot"){
  if(!(is.haplohh(haplohh))){stop("Data oject is not of valid haplohh object... (see data2haplohh() function)")}
  if(mrk<1 | mrk>haplohh@nsnp){stop(paste("Focal snp index must be between",1,"and",haplohh@nsnp))}
  if(limhaplo<2){stop("limhaplo must be >1")}
  if(limehh<0 | limehh>1){stop("limehh must be between 0 and 1")}
  
  ehh<-matrix(0,nrow=haplohh@nsnp,ncol=2) ; nhaplo_eval<-matrix(0,nrow=haplohh@nsnp,ncol=2) ; ihh<-rep(0,2)
  res.ehh<-.Fortran(name="r_ehh",
                    mrk=as.integer(mrk),
                    nmrk=as.integer(haplohh@nsnp),
                    nhap=as.integer(haplohh@nhap),
                    haplo=as.integer(haplohh@haplo),
                    map_pos=as.double(haplohh@position),
                    ehh=as.double(ehh),
                    out_nhaplo_eval=as.double(nhaplo_eval),
                    ihh=as.double(ihh),
                    limhaplo=as.integer(limhaplo),
                    limehh=as.double(limehh))
  ehh=matrix(res.ehh$ehh,2,haplohh@nsnp,byrow=T)
  nhaplo_eval=matrix(res.ehh$out_nhaplo_eval,2,haplohh@nsnp,byrow=T)
  rownames(ehh)=rownames(nhaplo_eval)=c("Anc. Allele","Der. Allele")
  colnames(ehh)=colnames(nhaplo_eval)=haplohh@snp.name
  ihh=res.ehh$ihh ; names(ihh)=c("Anc. Allele","Der. Allele")

 if(plotehh){
   sel_reg<-(colSums(nhaplo_eval)>0)
   if(sum(sel_reg)>0){
     matplot(haplohh@position[sel_reg],t(ehh[,sel_reg]),col=c("red","blue"),lty=1,
             type="l",main=main_leg,bty="n",xlab="Position",ylab="EHH")
     abline(v=haplohh@position[mrk],lty=2)
     smartlegend("left","top",c("Anc. Allele","Der. Allele"),col=c("red","blue"),bty="n",lty=1)
   }
   }

  return(list(ehh=ehh,nhaplo_eval=nhaplo_eval,freq_all1=nhaplo_eval[1,mrk]/sum(nhaplo_eval[,mrk]),ihh=ihh))

}
