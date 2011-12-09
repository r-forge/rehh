module ehh_utils
!---------------  EHH  -------------------------------! 
contains

subroutine ehh_calc(haplos_in,dir,all_test,tmp_ehh,tmp_neval)
implicit none
!dir=1 (1:gauche, 2:droite)
!tmp_neval=tmp_neval: vecteur pour le nombre d'haplos associé a all 1 et all 2
integer, intent(in)  :: dir,all_test
integer, intent(in), dimension(:,:) :: haplos_in
integer, intent(out) :: tmp_neval
real (kind=8), intent(out)    :: tmp_ehh
integer,allocatable  :: listehap(:,:),cnthapdiff(:) 
integer              :: nhapdiff, haptrouve,tmpnzero, tmp_nmrk, tmp_nhap,hapteste,nalldiff,i1,j1,tmp_posmrk

! initialisation
 tmp_nhap=size(haplos_in,1) ; tmp_nmrk=size(haplos_in,2)
 tmp_neval=0 ; tmp_ehh=0. ; nhapdiff=0
allocate(listehap(tmp_nhap,tmp_nmrk),cnthapdiff(tmp_nhap))
 cnthapdiff(:)=0 ; listehap(:,:)=0

if(dir==1) then 
  tmp_posmrk=1
else
  tmp_posmrk=tmp_nmrk
end if


!!LISTE HAPLOS
   do i1=1,tmp_nhap !2
     if(haplos_in(i1,tmp_posmrk)==all_test) then
      tmpnzero=0
!recherche de 0 dans l'haplo
      do j1=1,tmp_nmrk !3
       if(haplos_in(i1,j1)==0) tmpnzero=tmpnzero+1
      end do      !3
      if(tmpnzero==0) then
        if(nhapdiff==0) then
          nhapdiff=nhapdiff+1 ;  cnthapdiff(nhapdiff)=1
          do j1=1,tmp_nmrk !3
            listehap(nhapdiff,j1)=haplos_in(i1,j1)
          end do !3
        else
          haptrouve=0 ; hapteste=0
          do while (haptrouve==0 .and. hapteste<nhapdiff) !3                 
            hapteste=hapteste+1 ;  nalldiff=0 ;  j1=0
             do while (nalldiff==0 .and. j1<tmp_nmrk) !4
              j1=j1+1
              if(haplos_in(i1,j1)/=listehap(hapteste,j1)) nalldiff=1
             end do !4
            if(nalldiff==0) then 
             haptrouve=1 ;  cnthapdiff(hapteste)=cnthapdiff(hapteste)+1
            end if
          end do !3
          if(haptrouve==0) then
            nhapdiff=nhapdiff+1 ;  cnthapdiff(nhapdiff)=1 
            do j1=1,tmp_nmrk !3
             listehap(nhapdiff,j1)=haplos_in(i1,j1)
            end do !3
          end if        
        end if
     end if 
  end if
 end do !2
!!EHH
if(nhapdiff/=0) then
 do i1=1,nhapdiff
   tmp_neval=tmp_neval+cnthapdiff(i1)
   tmp_ehh=tmp_ehh+cnthapdiff(i1)*(cnthapdiff(i1)-1.)!cnthapdiff(i1)**2
 end do
   if(tmp_neval>1) then
    tmp_ehh=tmp_ehh/(tmp_neval*(tmp_neval-1.))
   else
    tmp_ehh=0.
   end if
end if
deallocate(listehap,cnthapdiff)

end subroutine ehh_calc

!---------------AIRE--------------------------------------! 

subroutine calc_aire(ordonne,abscisse,aire)
implicit none
real (kind=8), intent(in), dimension(:) :: ordonne,abscisse
real (kind=8), intent(out) :: aire
real (kind=8), allocatable :: hauteurs(:),longueurs(:)
integer :: i,tmp_nmrk

if(size(ordonne)/=size(abscisse)) then
print *,'size (ord)=',size(ordonne),' abs=',size(abscisse)
stop 'ERR: probleme calcul Aire:' 
end if
tmp_nmrk=size(ordonne)
allocate(hauteurs(tmp_nmrk-1),longueurs(tmp_nmrk-1))
do i=1,(tmp_nmrk-1)
  hauteurs(i)=abs(ordonne(i+1)-ordonne(i))
  longueurs(i)=abscisse(i+1)-abscisse(i)
end do   

aire=0.
do i=1,(tmp_nmrk-1)
 aire=aire + (longueurs(i)*hauteurs(i)/2) + longueurs(i)*ordonne(i) !Methode des trapezes!!!
end do
deallocate(hauteurs,longueurs) 
end subroutine calc_aire

!------------------------------------------------------------
!----------EHHS---------------------------------------------! 
!Selon la methode dans Tang et al. (2007): n'est plus allele specifique:
subroutine ehhs_calc(haplos_in,dir,tmp_ehhs,tmp_neval)
implicit none
!dir=1 (1:gauche, 2:droite) !pour avoir les frequences alleliques
integer, intent(in)  :: dir
integer, intent(in), dimension(:,:) :: haplos_in
integer, intent(out) :: tmp_neval
real (kind=8), intent(out)    :: tmp_ehhs
integer,allocatable  :: listehap(:,:),cnthapdiff(:) 
integer              :: nhapdiff, haptrouve,tmpnzero, tmp_nmrk, tmp_nhap,hapteste,nalldiff,i1,j1,tmp_posmrk,cnt_all(2)
real (kind=8)                :: hetero_hap,hetero_all,homo_all,homo_hap

! initialisation
 tmp_nhap=size(haplos_in,1) ; tmp_nmrk=size(haplos_in,2) ; nhapdiff=0
 tmp_neval=0 ; tmp_ehhs=0. ; homo_hap=0. ; homo_all=0. ; hetero_hap=0. ; hetero_all=0.
allocate(listehap(tmp_nhap,tmp_nmrk),cnthapdiff(tmp_nhap))
 cnthapdiff(:)=0 ; listehap(:,:)=0 ; cnt_all(:)=0

if(dir==1) then 
  tmp_posmrk=1
else
  tmp_posmrk=tmp_nmrk
end if
!!LISTE HAPLOS
   do i1=1,tmp_nhap !2
!recherche de 0 dans l'haplo
      tmpnzero=0
      do j1=1,tmp_nmrk !3
       if(haplos_in(i1,j1)==0) tmpnzero=tmpnzero+1
      end do      !3

      if(tmpnzero==0) then
        if(nhapdiff==0) then
          nhapdiff=nhapdiff+1 ;  cnthapdiff(nhapdiff)=1
          do j1=1,tmp_nmrk !3
            listehap(nhapdiff,j1)=haplos_in(i1,j1)
          end do !3
        else
          haptrouve=0 ; hapteste=0
          do while (haptrouve==0 .and. hapteste<nhapdiff) !3                 
            hapteste=hapteste+1 ;  nalldiff=0 ; j1=0
             do while (nalldiff==0 .and. j1<tmp_nmrk) !4
              j1=j1+1
              if(haplos_in(i1,j1)/=listehap(hapteste,j1)) nalldiff=1
             end do !4
            if(nalldiff==0) then 
             haptrouve=1 ; cnthapdiff(hapteste)=cnthapdiff(hapteste)+1
            end if
          end do !3
          if(haptrouve==0) then
            nhapdiff=nhapdiff+1 ; cnthapdiff(nhapdiff)=1 
            do j1=1,tmp_nmrk !3
             listehap(nhapdiff,j1)=haplos_in(i1,j1)
            end do !3
          end if        
        end if
     end if 
 end do !2

!!EHHS
if(nhapdiff/=0) then
 do i1=1,nhapdiff
   tmp_neval=tmp_neval+cnthapdiff(i1)
   homo_hap=homo_hap+cnthapdiff(i1)**2
!decompte des 2 alleles au marqueur
     if(listehap(i1,tmp_posmrk)==1) cnt_all(1)=cnt_all(1)+cnthapdiff(i1)
     if(listehap(i1,tmp_posmrk)==2) cnt_all(2)=cnt_all(2)+cnthapdiff(i1)    
 end do
   homo_all=cnt_all(1)**2 + cnt_all(2)**2
   homo_hap=homo_hap/(tmp_neval**2)
   homo_all=homo_all/(tmp_neval**2)
   hetero_hap=(1-homo_hap)*tmp_neval/(tmp_neval-1)
   hetero_all=( 1- homo_all) * tmp_neval/(tmp_neval-1)
   tmp_ehhs=(1-hetero_hap)/(1-hetero_all)
end if
deallocate(listehap,cnthapdiff)

end subroutine ehhs_calc


!---------------------------------------------------------!

end module ehh_utils
