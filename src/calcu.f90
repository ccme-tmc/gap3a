!BOP
!
! !ROUTINE: calcu
!
! !INTERFACE:
      subroutine calcu(iat,zz,l,el,uv,duv,nodes,isp)

! !DESCRIPTION:
!
! This subroutine calculates the radial wave function $u_l(r,E_l)$ and 
! normalizes it
!
! The result is stored in the global variables \verb"a" and \verb"b".
!
! !USES:

      use constants, only: cfein1, cfein2
      use radwf,     only: a, b,vr, rel
      use struk,     only: ro, nrpt, dh

! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: iat ! inequivalent atom
      integer(4), intent(in) :: isp ! index for spin
      integer(4), intent(in) :: l   ! azimutal quantum number of u
      real(8),    intent(in) :: el  ! linearization energy in Hartree
      real(8),    intent(in) :: zz  ! nuclear charge

! !OUTPUT PARAMETERS:

      real(8),    intent(out) :: duv   ! value of  d u_l(r,E)/d r at the  muffin tin radius
      real(8),    intent(out) :: uv    ! value of u_l(r,E_l) at the muffin tin radius.
      integer(4), intent(out) :: nodes ! number of nodes of u_l(r,E_l)


! !LOCAL VARIABLES:
      
      integer(4) :: npt ! Number of radial mesh points 
      integer(4) :: irp ! runs over the spherical grid points
    
      real(8) :: fl   ! Angular momentum l
      real(8) :: dx   ! Step size of the logarithmic radial mesh 
      real(8) :: ovlp ! square of the norm of u_l; |u_l(r,E_1=E_l-Delta E)|^2
      real(8) :: rnot ! first radial mesh point (for each atom the corresponding ro(iat) is stored here).
      real(8) :: trx  ! normalization factor of u_l; |u_l(r,E_1=E_l-Delta E)|^-1
 
! !EXTERNAL ROUTINES: 

      external outwin
      external rint13

! !INTRINSIC ROUTINES: 


      intrinsic sqrt

! !REVISION HISTORY:
! 
! Created Feb. 6th. 2004 by RGA
! Last modified July 20th 2004 by RGA
!
!EOP
!
!BOC
      fl = dble(l)
      rnot = ro(iat)
      npt = nrpt(iat)
      dx = dh(iat)

      call outwin(rel,vr(1,iat,isp),rnot,dx,npt,el,fl,uv,duv,nodes,zz)
!
!     Calculate |u_l(r,E_l)|^2:
!
      call rint13(cfein1,cfein2,a,b,a,b,ovlp,iat)
!
!     Calculate the normalization factor:
!
      trx = 1.0d+0/sqrt(ovlp)
!
!     Store the normalized values at the boundaries in p and dp:
!
      uv = trx * uv
      duv = trx * duv

#ifndef MPI
!     debug use - zmy
!      write(*, "(A6,I3,A6,F6.1,A15,F11.8)") &
!     &    "iat = ", iat, "; l = ", fl, "; el(Ha) = ", el
!      write(*, *) "nodes = ", nodes
!      do irp = 1, npt
!        !write(*, "(2ES19.17)") a(irp), b(irp)
!        write(*, *) rnot*exp(dx*(irp-1)), a(irp), b(irp)
!      enddo

!     debug use - zmy
#endif

!
!     Normalize the function u_l:
!
      do irp = 1, npt
        a(irp) = trx * a(irp)
        b(irp) = trx * b(irp)
      enddo

      end subroutine calcu

!EOC
