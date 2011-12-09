!version du 24 02 2011 (calul selon approche Sabeti) et quelques modfis de forme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calcul des ehh (fait appel a ehh_calc)
!!renvoie: ehh(nmrk,2) !ehh en chaque position pour les deux alleles et nhaplo_eval(1:nmrk,2)
!!renvoie ihh_1 et ihh_2
subroutine r_ehh(mrk,nmrk,nhap,haplo,map_pos,ehh,nhaplo_eval,ihh,limhaplo,limehh)
use ehh_utils
implicit none

integer::lim_mrk=0,nmrk, nhap, mrk, tmp,tmp1,tmp_neval,i,j,limhaplo,nmrk_g,nmrk_d,tmp_all, haplo(nhap,nmrk),nhaplo_eval(nmrk,2)
integer,allocatable :: haplo_test(:,:)
real (kind=8) :: limehh,tmp_ehh,tmp_aire,ihh(2),ehh(nmrk,2),map_pos(nmrk)
real (kind=8),allocatable ::tmp_abs(:)

 if(lim_mrk==0) then
 nmrk_g=nmrk ; nmrk_d=1
 else
 nmrk_g=min(nmrk,mrk+lim_mrk) ; nmrk_d=max(1,mrk-lim_mrk)
 end if
 
 ehh(:,:)=0. ; nhaplo_eval(:,:)=0 ; ehh(mrk,1:2)=1.
 do j=1,nhap
  if(haplo(j,mrk)==1) nhaplo_eval(mrk,1)=nhaplo_eval(mrk,1)+1
  if(haplo(j,mrk)==2) nhaplo_eval(mrk,2)=nhaplo_eval(mrk,2)+1
 end do

!!!!!!!!!!!!!!!!!!!!!!!!!!
!calcul des EHH à gauche 
!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!
!calcul des EHH à droite
!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Calcul des iHH        !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do tmp_all=1,2
 tmp_abs=ehh(:,tmp_all)
 call calc_aire(tmp_abs,map_pos(:),tmp_aire)
 ihh(tmp_all)=tmp_aire
end do

return

end subroutine r_ehh

