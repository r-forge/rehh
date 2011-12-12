!version du 24 02 2011 (calul selon approche Sabeti) et quelques modfis de forme
!patch pour windows ou on rentre nhaplo_eval en double
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calcul des ehh (fait appel a ehh_calc)
!!renvoie: ehh(nmrk,2) !ehh en chaque position pour les deux alleles et nhaplo_eval(1:nmrk,2)
!!renvoie ihh_1 et ihh_2
subroutine r_ehhs(mrk,nmrk,nhap,haplo,map_pos,ehhs,out_nhaplo_eval,ies,limhaplo,limehhs)
!haplo=as.integer => vecteur qu'il faut ensuite decomposer
!map_pos=as.single =>vecteur de longueur nmrk
!nmrk et nhap =>as.single
use ehh_utils

implicit none

integer::nmrk, nhap, mrk, tmp,tmp1,direction,dum_int=100,i,j,limhaplo,lim_mrk=0,&
         nmrk_g,nmrk_d,tmp_all, haplo(nhap,nmrk)
! si lim_nmrk=0, on fait comme avnat (on explore tout),sinon ion ne considere qe lim_mrk à agauche et a droite
integer,allocatable :: haplo_test(:,:),nhaplo_eval(:)
real (kind=8) :: limehhs,tmp_ehh=0.,tmp_aire,ies,dum_real,ehhs(nmrk),map_pos(nmrk),out_nhaplo_eval(nmrk)
real (kind=8),allocatable ::tmp_abs(:)!,map_pos(:)

allocate(nhaplo_eval(nmrk))

 if(lim_mrk==0) then
 nmrk_g=nmrk ; nmrk_d=1
 else
 nmrk_g=min(nmrk,mrk+lim_mrk) ; nmrk_d=max(1,mrk-lim_mrk)
 end if
 
 ehhs(:)=0. ; nhaplo_eval(:)=0 ; ehhs(mrk)=1.

 do j=1,nhap
  if(haplo(j,mrk)==1 .or. haplo(j,mrk)==2) nhaplo_eval(mrk)=nhaplo_eval(mrk)+1
 end do

! a gauche
 if(mrk<nmrk_g) then
  tmp=0
  if(nhaplo_eval(mrk)<limhaplo) then
      ehhs(mrk)=0. ; tmp=1
  end if
  j=0
  do while((mrk+j)<nmrk_g .and. tmp==0)
    j=j+1
    allocate(haplo_test(nhap,j+1))
    haplo_test=haplo(1:nhap,mrk:(mrk+j))
     call ehhs_calc(haplo_test,1,dum_real,dum_int)
    deallocate(haplo_test)
      if(dum_int>limhaplo .and. dum_real>limehhs) then
       ehhs(mrk+j)=dum_real ; nhaplo_eval(mrk+j)=dum_int
      else
       tmp=1
      end if
  end do
end if

! a droite
if(mrk>nmrk_d) then
  tmp=0
  if(nhaplo_eval(mrk)<limhaplo) then
      ehhs(mrk)=0. ; tmp=1
  end if
  j=0
  do while((mrk-j)>nmrk_d .and. tmp==0)
    j=j+1
    allocate(haplo_test(nhap,j+1))
    haplo_test=haplo(1:nhap,(mrk-j):mrk)
    call ehhs_calc(haplo_test,2,dum_real,dum_int)
    deallocate(haplo_test)
      if(dum_int>limhaplo .and. dum_real>limehhs) then
       ehhs(mrk-j)=dum_real ; nhaplo_eval(mrk-j)=dum_int
      else
       tmp=1
      end if
  end do
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Calcul des iES        !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
tmp_abs=ehhs(:)
 call calc_aire(tmp_abs,map_pos(:),tmp_aire)
ies=tmp_aire

out_nhaplo_eval(:)=nhaplo_eval(:) +0.

!print *,map_pos !ies,tmp_aire,ehhs(mrk)

!---------------------------------------------------------!
end subroutine r_ehhs

