!version du 24 02 2011 (calul selon approche Sabeti) et quelques modfis de forme
subroutine r_scan_hh(nmrk,nhap,haplo,map_pos,limhaplo,limehh,limehhs,out_res_all)
use ehh_utils
implicit none

integer:: lim_mrk=0,nmrk, nhap, mrk,dum_int, tmp,tmp1,tmp_neval,&
          i,j,limhaplo,nmrk_g,nmrk_d,tmp_all, haplo(nhap,nmrk)
integer,allocatable :: haplo_test(:,:),nhaplo_eval(:,:)
real (kind=8) :: limehh,limehhs,tmp_ehh,tmp_aire,map_pos(nmrk),dum_real,out_res_all(nmrk,4)
real (kind=8),allocatable ::tmp_abs(:),ehh(:,:)

allocate(ehh(nmrk,3),nhaplo_eval(nmrk,3),tmp_abs(nmrk))

do mrk=1,nmrk
 if(lim_mrk==0) then
 nmrk_g=nmrk ; nmrk_d=1
 else
 nmrk_g=min(nmrk,mrk+lim_mrk) ; nmrk_d=max(1,mrk-lim_mrk)
 end if
 
 ehh(:,:)=0. ; nhaplo_eval(:,:)=0 ; ehh(mrk,1:3)=1.
 do j=1,nhap
  if(haplo(j,mrk)==1) nhaplo_eval(mrk,1)=nhaplo_eval(mrk,1)+1
  if(haplo(j,mrk)==2) nhaplo_eval(mrk,2)=nhaplo_eval(mrk,2)+1
  nhaplo_eval(mrk,3)=nhaplo_eval(mrk,1)+nhaplo_eval(mrk,2)
 end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!EHH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!calcul des EHH à gauche 
 if(mrk<nmrk_g) then
  do tmp_all=1,2
  tmp=0 ! variable de controle pour l'allele tmp/=0 si on est inferieur a limhaplo ou limehh 
  if(nhaplo_eval(mrk,tmp_all)<limhaplo) then
      tmp=1 !le marqueur est monomorphe: on passe
      ehh(mrk,tmp_all)=0. !eviter de biaiser le calcul des aires
  end if
  j=0
  do while((mrk+j)<nmrk_g .and. tmp==0)
    j=j+1
    allocate(haplo_test(nhap,j+1))
    haplo_test=haplo(1:nhap,mrk:(mrk+j))
      call ehh_calc(haplo_test,1,tmp_all,tmp_ehh,tmp_neval) 
      if(tmp_neval>limhaplo .and. tmp_ehh>limehh) then
       ehh((mrk+j),tmp_all)=tmp_ehh ; nhaplo_eval((mrk+j),tmp_all)=tmp_neval
      else
       tmp=tmp+1
      end if
   deallocate(haplo_test)
  end do
  end do
 end if

!calcul des EHH à droite
if(mrk>nmrk_d) then
 do tmp_all=1,2
  tmp=0 ! variable de controle pour l'allele 1 (resp 2) tmp/=0 si on est inferieur a limhaplo ou limehh sur l'allele 1
  if(nhaplo_eval(mrk,tmp_all)<limhaplo) then
   tmp1=1 !le marqueur est monomorphe: on passe
   ehh(mrk,tmp_all)=0. !eviter de biaiser le calcul des aires
  end if
  j=0
  do while((mrk-j)>nmrk_d .and. tmp==0)
   j=j+1
   allocate(haplo_test(nhap,j+1))
   haplo_test=haplo(1:nhap,(mrk-j):mrk)
    call ehh_calc(haplo_test,2,tmp_all,tmp_ehh,tmp_neval) 
    if(tmp_neval>limhaplo .and. tmp_ehh>limehh) then
     ehh((mrk-j),tmp_all)=tmp_ehh ; nhaplo_eval((mrk-j),tmp_all)=tmp_neval
    else
     tmp=tmp+1
    end if
   deallocate(haplo_test)
  end do
 end do
end if

!!Calcul des iHH        !!!!!!!!!!!!!!!!!!!!!!!!
out_res_all(mrk,1)=(nhaplo_eval(mrk,1)+0.)/sum(nhaplo_eval(mrk,1:2))
do tmp_all=1,2
 tmp_abs=ehh(:,tmp_all)
 call calc_aire(tmp_abs,map_pos(:),tmp_aire)
 out_res_all(mrk,tmp_all+1)=tmp_aire
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!EHHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! a gauche
 if(mrk<nmrk_g) then
  tmp=0
  if(nhaplo_eval(mrk,3)<limhaplo) then
      ehh(mrk,3)=0. ; tmp=1
  end if
  j=0
  do while((mrk+j)<nmrk_g .and. tmp==0)
    j=j+1
    allocate(haplo_test(nhap,j+1))
    haplo_test=haplo(1:nhap,mrk:(mrk+j))
     call ehhs_calc(haplo_test,1,dum_real,dum_int)
    deallocate(haplo_test)
      if(dum_int>limhaplo .and. dum_real>limehhs) then
       ehh(mrk+j,3)=dum_real ; nhaplo_eval(mrk+j,3)=dum_int
      else
       tmp=1
      end if
  end do
end if

! a droite
if(mrk>nmrk_d) then
  tmp=0
  if(nhaplo_eval(mrk,3)<limhaplo) then
      ehh(mrk,3)=0. ; tmp=1
  end if
  j=0
  do while((mrk-j)>nmrk_d .and. tmp==0)
    j=j+1
    allocate(haplo_test(nhap,j+1))
    haplo_test=haplo(1:nhap,(mrk-j):mrk)
    call ehhs_calc(haplo_test,2,dum_real,dum_int)
    deallocate(haplo_test)
      if(dum_int>limhaplo .and. dum_real>limehhs) then
       ehh(mrk-j,3)=dum_real ; nhaplo_eval(mrk-j,3)=dum_int
      else
       tmp=1
      end if
  end do
end if
!!Calcul des iES        !!!!!!!!!!!!!!!!!!!!!!!!
tmp_abs=ehh(:,3)
 call calc_aire(tmp_abs,map_pos(:),tmp_aire)
out_res_all(mrk,4)=tmp_aire
!---------------------------------------------------------!

!write(*,'(A,1x,i7,1x,A,1x,i7,1x,A,f6.4,1x,3(f16.4),A)'),'SNP ',mrk,' out of ',nmrk,'(',out_res_all(mrk,:),')'
!if(mod(mrk,100)==0) write(*,'(A,1x,i7,1x,A,1x,i7)'),'SNP ',mrk,' out of ',nmrk
write(*,'(A,1x,i7,1x,A,1x,i7)'),'SNP ',mrk,' out of ',nmrk
end do

end subroutine r_scan_hh