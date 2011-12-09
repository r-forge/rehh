scan_hh<-function(haplohh,limhaplo=2,limehh=0.05,limehhs=0.05){
  if(!(is.haplohh(haplohh))){stop("Data oject is not of valid haplohh object... (see data2haplohh() function)")}
  if(limhaplo<2){stop("limhaplo must be >1")}
  if(limehh<0 | limehh>1){stop("limehh must be between 0 and 1")}
  if(limehhs<0 | limehhs>1){stop("limehhs must be between 0 and 1")}

    out_fort=matrix(0,haplohh@nsnp,4)
   
    res_scan<-.Fortran(name="r_scan_hh",
                       nmrk=as.integer(haplohh@nsnp),
                       nhap=as.integer(haplohh@nhap),
                       haplo=as.integer(haplohh@haplo),
                       map_pos=as.double(haplohh@position),
                       limhaplo=as.integer(limhaplo),
                       limehh=as.double(limehh),
                       limehhs=as.double(limehhs),
                       out_res_all=as.double(out_fort))

    RES_ALL=matrix(res_scan$out_res_all,haplohh@nsnp,4)
    RES_ALL=cbind(rep(haplohh@chr.name,haplohh@nsnp),haplohh@position,RES_ALL)
    rownames(RES_ALL)=haplohh@snp.name ; colnames(RES_ALL)=c("CHR","POSITION","FREQ_a","IHHa","IHHd","IES")
  return(RES_ALL)

}
