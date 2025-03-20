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
 subroutine excwig(rho,vxc,exc,mmax)

! Wigner ingerpolation formula for exchange-correlation potential and
! energy density

!rho  charge density
!vxc  exchange-correlation potential
!exc  exchange-correlation energy density
!mmax  dimension of log radial grid

 implicit none

 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: pi4=4.0d0*pi
 real(dp), parameter :: pi4i=1.0d0/pi4

!Input variables
 real(dp) ::  rho(mmax)
 integer :: mmax

!Output variables
 real(dp) :: vxc(mmax),exc(mmax)

!Local vafriables
 real(dp) ::  rh,xx
 integer :: ii

! Wigner interpolation formula, E. Wigner, Phys. Rev. 46, 1002 (1934).

 do ii=1,mmax
   rh=rho(ii)*pi4i
   xx=rh**0.33333333333333333d0
   vxc(ii)=-xx*(0.984 + (0.943656+8.8963*xx)/(1.0+12.57*xx)**2)
   exc(ii)=-0.738*xx*(1.0+0.959/(1.0+12.57*xx))
 end do
 return
 end subroutine excwig

