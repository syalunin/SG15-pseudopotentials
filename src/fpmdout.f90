!
! Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! interpolates various arrays onto linear radial mesh to create file
! for PWSCF input using the fpmd file format

 subroutine fpmdout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&                   zz,zion,mmax,iexc,icmod,nrl,drl,atsym,epstot, &
&                   na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&                   fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact, &
&                   epsh1,epsh2,depsh,rlmax,psfile)


!lmax  maximum angular momentum
!lloc  l for local potential
!rc  core radii
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
!nproj  number of vkb projectors for each l
!rr  log radial grid
!vpuns  unscreened semi-local pseudopotentials (vp(:,5) is local potential 
!  if linear combination is used)
!rho  valence pseudocharge
!rhomod  model core charge
!zz  atomic number
!zion  at this point, total valence charge (becomes psuedoion charge)
!mmax  size of log radial grid
!iexc  type of exchange-correlation
!icmod  1 if model core charge is used, otherwise 0
!nrl size of linear radial grid
!drl spacing of linear radial grid
!atsym  atomic symbol
!epstot  pseudoatom total energy
!psfile  should be 'fpmd'

 implicit none
 integer, parameter :: dp=kind(1.0d0)

 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp

!Input variables
 integer :: lmax,lloc,iexc,mmax,nrl,icmod
 integer :: nproj(6)
 real(dp) :: drl,fcfact,zz,zion,epstot
 real(dp) :: rr(mmax),vpuns(mmax,5),rho(mmax),vkb(mmax,2,4)
 real(dp) :: rhomod(mmax,5)
 real(dp):: rc(6),evkb(2,4)
 character*2 :: atsym

!additional input for fpmd output to echo input file, all as defined
! in the main progam
 integer :: na(30),la(30),ncon(6),nbas(6)
 integer :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 integer :: nc,nv,lpopt,ncnf
 real(dp) :: fa(30),rc0(6),ep(6),qcut(6),debl(6),facnf(30,5)
 real(dp) :: dvloc0,epsh1,epsh2,depsh,rlmax
 character*4 :: psfile

!Output variables - printing only

!Local variables
 integer :: ii,jj,ll,l1,iproj
 integer :: dtime(8)
 real(dp) :: dmat(2,2,lmax+1)
 real(dp), allocatable :: rhomodl(:,:)
 real(dp),allocatable :: rhol(:),rl(:),vkbl(:,:,:),vpl(:,:)
 real(dp), allocatable :: vkborl(:,:,:)
 character*2 :: pspd(3)


 allocate(rhol(nrl),rl(nrl),vkbl(nrl,2,4),vpl(nrl,5),rhomodl(nrl,5))
 allocate(vkborl(mmax,2,4))

! divide Kleinman-Bylander / Vanderbilt potentials by their r -> 0 dependence

 do l1 = 1, lmax + 1
   ll = l1 - 1
   if (ll .ne. lloc) then
     vkborl(:,1,l1)=vkb(:,1,l1)/rr(:)**(ll+1)
     vkborl(:,2,l1)=vkb(:,2,l1)/rr(:)**(ll+1)
   end if
 end do
!
! interpolation of everything onto linear output mesh
! cubic extrapolation to zero from nearest 4 points of linear mesh
!
 do  ii=1,nrl
   rl(ii)=drl*dble(ii-1)
 end do
!
 do l1=1,max(lmax+1,lloc+1)
   call dp3int(rr,vpuns(1,l1),mmax,rl,vpl(1,l1),nrl)
   vpl(1,l1)=4.0d0*vpl(2,l1)-6.0d0*vpl(3,l1)  &
&           +4.0d0*vpl(4,l1)-      vpl(5,l1)
   if(l1 .ne. lloc + 1) then

     call dp3int(rr,vkborl(1,1,l1),mmax,rl,vkbl(1,1,l1),nrl)
     call dp3int(rr,vkborl(1,2,l1),mmax,rl,vkbl(1,2,l1),nrl)

     vkbl(1,1,l1)=4.0d0*vkbl(2,1,l1)-6.0d0*vkbl(3,1,l1)  &
&                +4.0d0*vkbl(4,1,l1)-      vkbl(5,1,l1)
     vkbl(1,2,l1)=4.0d0*vkbl(2,2,l1)-6.0d0*vkbl(3,2,l1)  &
&                +4.0d0*vkbl(4,2,l1)-      vkbl(5,2,l1)
   end if
 end do
 
! rescale linear-mesh vkbl to restore actual r dependence

 do l1=1,lmax+1
   ll=l1-1
   if(ll.ne.lloc)then
     do ii = 1, nrl
       ! UPF and fpmd differ by factor r
       vkbl(ii,1,l1)=vkbl(ii,1,l1)*rl(ii)**ll
       vkbl(ii,2,l1)=vkbl(ii,2,l1)*rl(ii)**ll
     end do
   end if
 end do

 call dp3int(rr,rho,mmax,rl,rhol,nrl)

 rhol(1)= 4.0d0*rhol(2)-6.0d0*rhol(3) &
&        +4.0d0*rhol(4)-      rhol(5)

 do jj=1,5
   call dp3int(rr,rhomod(1,jj),mmax,rl,rhomodl(1,jj),nrl)

   rhomodl(1,jj)=4.0d0*rhomodl(2,jj)-6.0d0*rhomodl(3,jj) &
&               +4.0d0*rhomodl(4,jj)-      rhomodl(5,jj)
 end do


 call date_and_time(VALUES=dtime)
 ii=dtime(1)-2000
 if(ii<10) then
   write(pspd(1),'(a,i1)') '0',ii
 else
   write(pspd(1),'(i2)') ii
 end if
 ii=dtime(2)
 if(ii<10) then
   write(pspd(2),'(a,i1)') '0',ii
 else
   write(pspd(2),'(i2)') ii
 end if
 ii=dtime(3)
 if(ii<10) then
   write(pspd(3),'(a,i1)') '0',ii
 else
   write(pspd(3),'(i2)') ii
 end if

! section for fpmd output for pwscf

 write(6,'(/a)') 'Begin PSP_FPMD'
 write(6,'(a,/a,/a,/a,/a)') &
&  '<?xml version="1.0" encoding="UTF-8"?>', &
&  '<fpmd:species xmlns:fpmd="http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0"', &
&  '  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"', &
&  '  xsi:schemaLocation="http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0', &
&  '  species.xsd">', &
&  '<description>'
 write(6,'(/t2,a/t2,a/t2,a/t2,a/t2,a/t2,a//)') &
       'This pseudopotential file has been produced using the code', &
&      'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&      'scalar-relativistic version 2.1.1, 03/26/2014 by D. R. Hamann', &
&      'The code is available through a link at URL www.mat-simresearch.com.', &
&      'Documentation with the package provides a full discription of the', &
&      'input data below.'

 write(6,'(t2,a/t2,a/t2,a/t2,a//)') &
&      'While it is not required under the terms of the GNU GPL, it is',&
&      'suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)',&
&      'in any publication using these pseudopotentials.' 

 write(6,'(/a/)') &
&      'Input file for PP generation:'

! output printing (echos input data)

 write(6,'(a)') '# ATOM AND REFERENCE CONFIGURATION'
 write(6,'(a)') '# atsym  z    nc    nv    iexc   psfile'
 write(6,'(a,a,f6.2,3i6,2a)') '  ',trim(atsym),zz,nc,nv,iexc, &
&      '      ',trim(psfile)
 write(6,'(a/a)') '#','#   n    l    f        energy (Ha)'
 do ii=1,nc+nv
   write(6,'(2i5,f8.2)') na(ii),la(ii),fa(ii)
 end do

 write(6,'(a/a/a)') '#','# PSEUDOPOTENTIAL AND OPTIMIZATION','# lmax'
 write(6,'(i5)')  lmax
 write(6,'(a/a)') '#','#   l,   rc,     ep,   ncon, nbas, qcut'
 do l1=1,lmax+1
   write(6,'(i5,2f10.5,2i5,f10.5)') l1-1,rc(l1),ep(l1),ncon(l1),nbas(l1),qcut(l1)
 end do

 write(6,'(a/a/a,a)') '#','# LOCAL POTENTIAL','# lloc, lpopt,  rc(5),', &
&      '   dvloc0'
 write(6,'(2i5,f10.5,a,f10.5)') lloc,lpopt,rc(5),'   ',dvloc0

 write(6,'(a/a/a)') '#','# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
&      '# l, nproj, debl'
 do l1=1,lmax+1
   write(6,'(2i5,f10.5)') l1-1,nproj(l1),debl(l1)
 end do

 write(6,'(a/a/a)') '#','# MODEL CORE CHARGE', &
&      '# icmod, fcfact'
 write(6,'(i5,f10.5,2i5,f10.5)') icmod,fcfact

 write(6,'(a/a/a)') '#','# LOG DERIVATIVE ANALYSIS', &
&      '# epsh1, epsh2, depsh'
 write(6,'(3f8.2)') epsh1,epsh2,depsh

 write(6,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl'
 write(6,'(2f8.2)') rlmax,drl

 write(6,'(a/a/a)') '#','# TEST CONFIGURATIONS','# ncnf'
 write(6,'(i5)') ncnf
 write(6,'(a/a)') '# nvcnf','#   n    l    f'
 do jj=2,ncnf+1
   write(6,'(i5)') nvcnf(jj)
   do ii=nc+1,nc+nvcnf(jj)
     write(6,'(2i5,f8.2)') nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
   end do
   write(6,'(a)') '#'
 end do
 write(6,'(a)') &
&      '</description>'

 ! general information about the atom
 write(6,'(a,a,a)') '<symbol>',atsym,'</symbol>'
 write(6,'(a,i3,a)') '<atomic_number>',nint(zz),'</atomic_number>'
 write(6,'(a,f15.8,a)') '<mass>',mass(nint(zz)),'</mass>'
 write(6,'(a)') '<norm_conserving_semilocal_pseudopotential>'
 write(6,'(a,i6,a)') '<valence_charge>',nint(zion),'</valence_charge>'
 write(6,'(a,f10.4,a)') '<mesh_spacing>',drl,'</mesh_spacing>'
 ! nonlinear core correction
 if (icmod == 1) then
   write(6,'(a,i4,a)') '<core_density size="',nrl,'">'
   write(6,'(1p,4e20.10)') (rhomodl(ii,1)/(4.0d0*pi),ii=1,nrl)
   write(6,'(a)') '</core_density>'
 end if
 ! local potential
 write(6,'(a,i4,a)') '<local_potential size="',nrl,'">'
 ! write local potential in Hartree units
 write(6,'(1p,4e20.10)') (vpl(ii,lloc+1),ii=1,nrl)
 write(6,'(a)') '</local_potential>'
 ! loop on angular mommentum for projector outputs
 dmat = 0
 do l1=1,lmax+1
   if(l1==lloc+1) cycle
   do jj=1,nproj(l1)
     iproj=iproj+1
     dmat(jj,jj,l1)=evkb(jj,l1)
     write(6,'(a,i1,a,i1,a,i4,a)') '<projector l="',l1-1,'" i="',jj,'" size="',nrl,'">'
     write(6,'(1p,4e20.10)') (vkbl(ii,jj,l1),ii=1,nrl)
     write(6,'(a)') '</projector>'
   end do
 end do
 ! D_ij matrix per angular moment
 do l1=1,lmax+1
   do ii=1,nproj(l1)
     do jj=1,nproj(l1)
       write(6,'(a,i1,a,i1,a,i1,a,e20.10,a)') &
& '<d_ij l="',l1-1,'" i="',ii,'" j="',jj,'">',dmat(ii,jj,l1),'</d_ij>'
     end do
   end do
 end do

 write(6,'(a)') '</norm_conserving_semilocal_pseudopotential>'
 write(6,'(a)') '</fpmd:species>'

 write(6,'(a)') 'End PSP_FPMD'
 deallocate(rhol,rl,vkbl,vpl,rhomodl)
 deallocate(vkborl)

 return

 contains

 ! mass = mass of ion
 real(dp) function mass(zzi)

 implicit none

 integer, intent(in) :: zzi

 select case (zzi)
 case(  1) ! H
   mass = 1.008
 case(  2) ! He
   mass = 4.0026
 case(  3) ! Li
   mass = 6.94
 case(  4) ! Be
   mass = 9.0122
 case(  5) ! B
   mass = 10.81
 case(  6) ! C
   mass = 12.011
 case(  7) ! N
   mass = 14.007
 case(  8) ! O
   mass = 15.999
 case(  9) ! F
   mass = 18.998
 case( 10) ! Ne
   mass = 20.180
 case( 11) ! Na
   mass = 22.990
 case( 12) ! Mg
   mass = 24.305
 case( 13) ! Al
   mass = 26.982
 case( 14) ! Si
   mass = 28.085
 case( 15) ! P
   mass = 30.974
 case( 16) ! S
   mass = 32.06
 case( 17) ! Cl
   mass = 35.45
 case( 18) ! Ar
   mass = 39.948
 case( 19) ! K
   mass = 39.098
 case( 20) ! Ca
   mass = 40.078
 case( 21) ! Sc
   mass = 44.956
 case( 22) ! Ti
   mass = 47.867
 case( 23) ! V
   mass = 50.942
 case( 24) ! Cr
   mass = 51.996
 case( 25) ! Mn
   mass = 54.938
 case( 26) ! Fe
   mass = 55.845
 case( 27) ! Co
   mass = 58.933
 case( 28) ! Ni
   mass = 58.693
 case( 29) ! Cu
   mass = 63.546
 case( 30) ! Zn
   mass = 65.38
 case( 31) ! Ga
   mass = 69.723
 case( 32) ! Ge
   mass = 72.63
 case( 33) ! As
   mass = 74.922
 case( 34) ! Se
   mass = 78.96
 case( 35) ! Br
   mass = 79.904
 case( 36) ! Kr
   mass = 83.798
 case( 37) ! Rb
   mass = 85.468
 case( 38) ! Sr
   mass = 87.62
 case( 39) ! Y
   mass = 88.906
 case( 40) ! Zr
   mass = 91.224
 case( 41) ! Nb
   mass = 92.906
 case( 42) ! Mo
   mass = 95.96
 case( 43) ! Tc
   mass = 97.91
 case( 44) ! Ru
   mass = 101.07
 case( 45) ! Rh
   mass = 102.91
 case( 46) ! Pd
   mass = 106.42
 case( 47) ! Ag
   mass = 107.87
 case( 48) ! Cd
   mass = 112.41
 case( 49) ! In
   mass = 114.82
 case( 50) ! Sn
   mass = 118.71
 case( 51) ! Sb
   mass = 121.76
 case( 52) ! Te
   mass = 127.60
 case( 53) ! I
   mass = 126.90
 case( 54) ! Xe
   mass = 131.29
 case( 55) ! Cs
   mass = 132.91
 case( 56) ! Ba
   mass = 137.33
 case( 71) ! Lu
   mass = 174.97
 case( 72) ! Hf
   mass = 178.49
 case( 73) ! Ta
   mass = 180.95
 case( 74) ! W
   mass = 183.84
 case( 75) ! Re
   mass = 186.21
 case( 76) ! Os
   mass = 190.23
 case( 77) ! Ir
   mass = 192.22
 case( 78) ! Pt
   mass = 195.08
 case( 79) ! Au
   mass = 196.97
 case( 80) ! Hg
   mass = 200.59
 case( 81) ! Tl
   mass = 204.38
 case( 82) ! Pb
   mass = 207.2
 case( 83) ! Bi
   mass = 208.98
 case( 84) ! Po
   mass = 208.98
 case( 85) ! At
   mass = 209.99
 case( 86) ! Rn
   mass = 222.02
 case( 87) ! Fr
   mass = 223.02
 case( 88) ! Ra
   mass = 226.03
 case(103) ! Lr
   mass = 262.11
 case(104) ! Rf
   mass = 265.12
 case(105) ! Db
   mass = 268.13
 case(106) ! Sg
   mass = 271.13
 case(107) ! Bh
   mass = 270
 case(108) ! Hs
   mass = 277.15
 case(109) ! Mt
   mass = 276.15
 case(110) ! Ds
   mass = 281.16
 case(111) ! Rg
   mass = 280.16
 case(112) ! Cn
   mass = 285.17
 case(113) ! Uut
   mass = 284.18
 case(114) ! Fl
   mass = 289.19
 case(115) ! Uup
   mass = 288.19
 case(116) ! Lv
   mass = 293
 case(117) ! Uus
   mass = 294
 case(118) ! Uuo
   mass = 294
 ! Lanthanoids
 case( 57) ! La
   mass = 138.91
 case( 58) ! Ce
   mass = 140.12
 case( 59) ! Pr
   mass = 140.91
 case( 60) ! Nd
   mass = 144.24
 case( 61) ! Pm
   mass = 144.91
 case( 62) ! Sm
   mass = 150.36
 case( 63) ! Eu
   mass = 151.96
 case( 64) ! Gd
   mass = 157.25
 case( 65) ! Tb
   mass = 158.93
 case( 66) ! Dy
   mass = 162.50
 case( 67) ! Ho
   mass = 164.93
 case( 68) ! Er
   mass = 167.26
 case( 69) ! Tm
   mass = 168.93
 case( 70) ! Yb
   mass = 173.05
 ! Actinoids
 case( 89) ! Ac
   mass = 227.03
 case( 90) ! Th
   mass = 232.04
 case( 91) ! Pa
   mass = 231.04
 case( 92) ! U
   mass = 238.03
 case( 93) ! Np
   mass = 237.05
 case( 94) ! Pu
   mass = 244.06
 case( 95) ! Am
   mass = 243.06
 case( 96) ! Cm
   mass = 247.07
 case( 97) ! Bk
   mass = 247.07
 case( 98) ! Cf
   mass = 251.08
 case( 99) ! Es
   mass = 252.08
 case(100) ! Fm
   mass = 257.10
 case(101) ! Md
   mass = 258.10
 case(102) ! No
   mass = 259.10
 case default
   mass = 0
 end select

 end function mass

 end subroutine fpmdout
