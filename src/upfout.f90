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
! for PWSCF input using the UPF file format

 subroutine upfout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod,rhoc,uuc, &
&                  zz,zion,mmax,iexc,icmod,nrl,drl,atsym,epstot, &
&                  na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&                  fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact, &
&                  epsh1,epsh2,depsh,rlmax,psfile)


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
!psfile  should be 'upf'
!rhoc    all-electron core density

 implicit none
 integer, parameter :: dp=kind(1.0d0)

 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp

!Input variables
 integer :: lmax,lloc,iexc,mmax,nrl,icmod
 integer :: nproj(6)
 real(dp) :: drl,fcfact,zz,zion,epstot
 real(dp) :: rr(mmax),vpuns(mmax,5),rho(mmax),vkb(mmax,2,4)
 real(dp) :: rhomod(mmax,5),rhoc(mmax),uuc(mmax,nc)
 real(dp):: rc(6),evkb(2,4)
 character*2 :: atsym

!additional input for upf output to echo input file, all as defined
! in the main progam
 integer :: na(30),la(30),ncon(6),nbas(6)
 integer :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 integer :: nc,nv,lpopt,ncnf
 real(dp) :: fa(30),rc0(6),ep(6),qcut(6),debl(6),facnf(30,5)
 real(dp) :: dvloc0,epsh1,epsh2,depsh,rlmax
 character*4 :: psfile

!Output variables - printing only

!Local variables
 integer :: ii,jj,ll,l1,iproj,ntotproj,nrlproj
 integer :: dtime(8)
 real(dp) :: dmat(8,8)
 real(dp), allocatable :: rhomodl(:,:),rhocl(:)
 real(dp), allocatable :: rhol(:),rl(:),vkbl(:,:,:),vpl(:,:)
 real(dp), allocatable :: vkborl(:,:,:),uucl(:,:)
 character*5 :: lnames
 character*2 :: pspd(3)

 lnames = "SPDFG"

 allocate(rhol(nrl),rl(nrl),vkbl(nrl,2,4),vpl(nrl,5),rhomodl(nrl,5),rhocl(nrl))
 allocate(vkborl(mmax,2,4),uucl(nrl,nc))

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
       vkbl(ii,1,l1)=vkbl(ii,1,l1)*rl(ii)**(ll+1)
       vkbl(ii,2,l1)=vkbl(ii,2,l1)*rl(ii)**(ll+1)
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

! section for upf output for pwscf

 write(6,'(/a)') 'Begin PSP_UPF'

 write(6,'(a,/a)') &
&      '<UPF version="2.0.1">', &
&      '  <PP_INFO>'
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

 write(6,'(a)') &
&      '    <PP_INPUTFILE>'

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
&      '    </PP_INPUTFILE>'

 write(6,'(a,3(/a))') &
&      '  </PP_INFO>', &
&      '  <!--                               -->', &
&      '  <!-- END OF HUMAN READABLE SECTION -->', &
&      '  <!--                               -->'

 write(6,'(a)') &
&      '    <PP_HEADER'

   write(6,'(t8,a)') &
&        'generated="Generated using ONCVPSP code by D. R. Hamann"'
   write(6,'(t8,a)') &
&        'author="anonymous"'
   write(6,'(t8,5a)') &
&        'date="',pspd(:),'"'
   write(6,'(t8,a)') &
&        'comment=""'
   write(6,'(t8,a,a2,a)') &
&        'element="',atsym,'"'
   write(6,'(t8,a)') &
&        'pseudo_type="NC"'
   write(6,'(t8,a)') &
&        'relativistic="scalar"'
   write(6,'(t8,a)') &
&        'is_ultrasoft="F"'
   write(6,'(t8,a)') &
&        'is_paw="F"'
   write(6,'(t8,a)') &
&        'is_coulomb="F"'
   write(6,'(t8,a)') &
&        'has_so="F"'
   write(6,'(t8,a)') &
&        'has_wfc="F"'
   write(6,'(t8,a)') &
&        'has_gipaw="T"'

   if(icmod==1) then
     write(6,'(t8,a)') &
&        'core_correction="T"'
   else
     write(6,'(t8,a)') &
&        'core_correction="F"'
   end if

   if(iexc==3) then
     write(6,'(t8,a)') &
&        'functional="LDA"'
   else if(iexc==4) then
     write(6,'(t8,a)') &
&        'functional="PBE"'
   else
     write(6,'(t8,a)') &
&        'upfout:  only iexc=3 or iexc=4 are supported for UPF output'
     stop
   end if

   write(6,'(t8,a,f8.2,a)') &
&        'z_valence="',zion,'"'

   write(6,'(t8,a,1p,e20.11,a)') &
&        'total_psenergy="',epstot,'"'

   write(6,'(t8,a,1p,e20.11,a)') &
&        'rho_cutoff="',rl(nrl),'"'

   write(6,'(t8,a,i1,a)') &
&        'l_max="',lmax,'"'

   if(lloc==4) then
     write(6,'(t8,a)') &
&        'l_local="-1"'
   else
     write(6,'(t8,a,i1,a)') &
&        'l_local="',lloc,'"'
   end if

! calculate total number of projectors and maximum linear mesh point for
! projectors

   ntotproj=0
   nrlproj=0
   do l1=1,lmax+1
     if(l1/=lloc+1) then
       ntotproj=ntotproj+nproj(l1)
     end if
     nrlproj=max(nrlproj,4+int(rc(l1)/drl))
   end do
   if(mod(nrlproj,2)/=0) nrlproj=nrlproj+1

   write(6,'(t8,a,i6,a)') &
&        'mesh_size="',nrl,'"'

   write(6,'(t8,a)') &
&        'number_of_wfc="0"'

   write(6,'(t8,a,i1,a)') &
&        'number_of_proj="',ntotproj,'"/>' !end of PP_HEADER

   write(6,'(t2,a)') &
&        '<PP_MESH>'

 write(6,'(t4,a,i4,a)') &
&      '<PP_R type="real"  size="',nrl,'" columns="8">'

 write(6,'(8f10.4)') (rl(ii),ii=1,nrl)

 write(6,'(t4a)') &
&      '</PP_R>'

 write(6,'(t4a,i4,a)') &
&      '<PP_RAB type="real"  size="',nrl,'" columns="8">'

 write(6,'(8f10.4)') (drl,ii=1,nrl)

 write(6,'(t4a)') &
&      '</PP_RAB>'

   write(6,'(t2,a)') &
&        '</PP_MESH>'

! write local potential with factor of 2 for Rydberg units
 write(6,'(a,i4,a)') &
&      '  <PP_LOCAL type="real"  size="',nrl,'" columns="4">'

 write(6,'(1p,4e20.10)') (2.0d0*vpl(ii,lloc+1),ii=1,nrl)

 write(6,'(a)') &
&      '  </PP_LOCAL>'

! loop on angular mommentum for projector outputs

 write(6,'(t2,a)') &
&      '<PP_NONLOCAL>'

 dmat(:,:)=0.0d0
 iproj=0

 do l1=1,lmax+1
   if(l1==lloc+1) cycle
   do jj=1,nproj(l1)
     iproj=iproj+1
     dmat(iproj,iproj)=2.0d0*evkb(jj,l1) !2 for Rydbergs
     write(6,'(t4,a,i1)') &
&          '<PP_BETA.',iproj
     write(6,'(t8,a)') &
&          'type="real"'
     write(6,'(t8,a,i4,a)') &
&          'size="',nrl,'"'
     write(6,'(t8,a)') &
&          'columns="4"'
     write(6,'(t8,a,i1,a)') &
&          'index="',iproj,'"'
!     write(6,'(t8,a)') &
!&          'label="3P"'
     write(6,'(t8,a,i1,a)') &
&          'angular_momentum="',l1-1,'"'
     write(6,'(t8,a,i4,a)') &
&          'cutoff_radius_index="',nrlproj,'"'
     write(6,'(t8,a,1p,e20.10,a)') &
&          'cutoff_radius="',(nrlproj-1)*drl,'" >'

     write(6,'(1p,4e20.10)') (vkbl(ii,jj,l1),ii=1,nrl)

     write(6,'(t4,a,i1,a)') &
&          '</PP_BETA.',iproj,'>'
   end do
 end do

 write(6,'(t4,a,i4,a)') &
&      '<PP_DIJ type="real"  size="',ntotproj**2,'" columns="4">'

 write(6,'(1p,4e20.10)') ((dmat(ii,jj),ii=1,ntotproj),jj=1,ntotproj)

 write(6,'(t4,a)') &
&      '</PP_DIJ>'

 write(6,'(t2,a)') &
&      '</PP_NONLOCAL>'

 write(6,'(t2,a)') &
&      '<PP_PSWFC>'
 write(6,'(t2,a)') &
&      '</PP_PSWFC>'

 if(icmod==1) then
   write(6,'(t2,a,i4,a)') &
&        '<PP_NLCC type="real"  size="',nrl,'" columns="4">'

   write(6,'(1p,4e20.10)') (rhomodl(ii,1)/(4.0d0*pi),ii=1,nrl)

   write(6,'(t2,a)') &
&        '</PP_NLCC>'
 end if

 write(6,'(t2,a,i4,a)') &
&      '<PP_RHOATOM type="real"  size="',nrl,'" columns="4">'

 write(6,'(1p,4e20.10)') ((rl(ii)**2)*rhol(ii),ii=1,nrl)

 write(6,'(t2,a)') &
&      '</PP_RHOATOM>'


! here begins the GIPAW section
 write(6,'(t2,a,i0,a)') '<PP_GIPAW gipaw_data_format="',2,'">'
 write(6,'(t2,a,i0,a)') '<PP_GIPAW_CORE_ORBITALS number_of_core_orbitals="',nc,'">'

 do ii=1,nc
    l1 = la(ii)
    call dp3int(rr,uuc(1,ii),mmax,rl,uucl(1,ii),nrl)
    write(6,'(t4,a,i0)') &
&       '<PP_GIPAW_CORE_ORBITAL.',ii
    write(6,'(t8,a)') &
&       'type="real"'
    write(6,'(t8,a,i0,a)') &
&       'size="',nrl,'"'
    write(6,'(t8,a)') &
&       'columns="4"'
    write(6,'(t8,a,i0,a)') &
&       'index="',ii,'"'
    write(6,'(t8,a,i0,a,a)') &
&       'label="',na(ii),lnames(l1+1:l1+1),'"'
    write(6,'(t8,a,i0,a,i0,a)') &
&       'n="',na(ii),'" l="',l1,'">'
    write(6,'(1p,4e20.10)') &
&       (rl(jj) * uucl(jj,ii),jj=1,nrl)
    write(6,'(t4,a,i0,a)') &
&       '</PP_GIPAW_CORE_ORBITAL.',ii,'>'

 end do
 write(6,'(t2,a)') '</PP_GIPAW_CORE_ORBITALS>'
 write(6,'(t2,a)') '</PP_GIPAW>'


 write(6,'(a)') &
&      '</UPF>'
  
 write(6,'(a)') 'End PSP_UPF'

 deallocate(rhol,rl,vkbl,vpl,rhomodl)
 deallocate(vkborl)

 return
 end subroutine upfout
