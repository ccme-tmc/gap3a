!BOP
!
! !ROUTINE: int1ipw
!
! !INTERFACE:
      subroutine int1ipw(ig,integral)
      
! !DESCRIPTION:

!This function calculates the integral of a plane wave with wave vector
!belonging to the reciprocal Bravais lattice in the
!interstitial region by the difference between the integral over the whole
!unit cell and the Muffin Tin spheres:
! !USES:
      
      use constants, only: czero,pi
      use struk,     only: nat,mult, pos, rmt, vmt

! !INPUT PARAMETERS:

      implicit none
      integer(4) :: ig(3)
      
! !OUTPUT PARAMETERS:
      
      complex(8) :: integral

! !LOCAL VARIABLES:      
      
      integer(4) :: i, iat,ieq,idf
      
      real(8), dimension(3) :: g  ! The reciprocal lattice vector.

      real(8) :: glen ! Length of g
      real(8) :: gr   ! glen times MT radius
      real(8) :: cgr  ! cos(gr)
      real(8) :: sgr  ! sin(gr)
      real(8) :: intmod ! Modulus of the MT integral
      real(8) :: phase  ! the phase of the MT integral
      real(8) :: mtintr ! The integral over the MT Sphere (real part)
      real(8) :: mtinti ! The integral over the MT Sphere (imaginary part)

      real(8) :: integr
      real(8) :: integi

! !EXTERNAL ROUTINES:      
      
      external k2cart

! !REVISION HISTORY:
!
! Created 30. April 2004 by RGA
!
!EOP
!BOC
      integr=0.0d0
      integi=0.0d0
      call k2cart(ig,1,g)
      glen=dsqrt(g(1)*g(1)+g(2)*g(2)+g(3)*g(3))
      if(glen.gt.1.0d-10)then
        idf=0
        do iat=1,nat
          mtintr=0.0d0
          mtinti=0.0d0
          gr=glen*rmt(iat)
          sgr=dsin(gr)
          cgr=dcos(gr)
          intmod = 3.0d0*vmt(iat)*(sgr/gr - cgr)/gr/gr
          do ieq=1,mult(iat)
            idf=idf+1
            phase=0.0d0
            do i=1,3
              phase=phase+dble(ig(i))*pos(i,idf)
            enddo
            mtintr=mtintr+dcos(2.0d0*pi*phase)
            mtinti=mtinti+dsin(2.0d0*pi*phase)
          enddo
          integr=integr-intmod*mtintr
          integi=integi-intmod*mtinti
        enddo
        integral =cmplx(integr,integi,8)  
      else  
        integral=cmplx(1.0d0,0.0d0,8)
        do iat=1,nat
          integral=integral-cmplx(vmt(iat)*dble(mult(iat)),0.0d0,8)
        enddo
      endif  

      return
      end subroutine int1ipw
!EOC      
