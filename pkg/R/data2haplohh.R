setClass(Class = "haplohh",
	representation(haplo = "matrix",position = "numeric",snp.name="character",chr.name = "numeric",nhap = "numeric",nsnp = "numeric")
)

is.haplohh <- function (x) 
{
    res <- (is(x,"haplohh") & validObject(x))
    return(res)
}

data2haplohh<-function(haplo,map){
   res<-new("haplohh")
   res@haplo<-as.matrix(haplo)

   if(sum(!(res@haplo==0 | res@haplo==1 | res@haplo==2 ))>0){
    stop("Alleles in haplotype must be coded as 0 (missing data), 1 (ancestral allele) or 2 (derived allele)")
    }

   res@nhap<-nrow(res@haplo) ; res@nsnp<-ncol(res@haplo)

   tmp_nom<-rownames(map)
   if(length(tmp_nom)!=res@nsnp){stop("Length of map and haplo matrices differ")}
   tmp_chr=unique(as.numeric(map[,1]))
   if(length(tmp_chr)>1){stop("Several chromosome names in map matrix")}
   res@chr.name<-tmp_chr[1]
   res@snp.name<-tmp_nom

   tmp_pos<-as.numeric(map[,2])
   if(sum(diff(tmp_pos)<0)>0){
     stop("SNP should be ordered on the map, check also that both haplotypes (column of haplo) and map (row of map) are ordered in the same way")}
   if(sum(diff(tmp_pos)==0)>0){
     warning("Some SNPs map to the same position")}

   res@position<-tmp_pos

  return(res)
}

