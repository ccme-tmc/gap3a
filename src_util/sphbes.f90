!BOP
!
! !ROUTINE: sphbes
!
! !INTERFACE:
      subroutine sphbes(n,x,sj,sjp)

! !DESCRIPTION:
!
!   Calculates the spherical bessel function $j_l(x)$ and its first derivative $j'_l(x)$ 
!   at the point $x\le 0$, for l=0,1,...n       
!
! !INPUT PARAMETERS:
     
      implicit none
      
      integer(4), intent(in) :: n ! Maximum l to be calculated
      real(8), intent(in) :: x    ! The argument at which the bessel functions are calculated

! !OUTPUT PARAMETERS:

      real(8), intent(out) :: sj(*) ! values of spherical bessel function
      real(8), intent(out) :: sjp(*)! values of first derivative of the spherical bessel function

! !LOCAL VARIABLES:

      integer(4) :: l
      real(8) :: factor  ! Multiplying factor between the spherical bessel func and the regular bessel func of half integer order
      real(8) :: rj(n+1) ! Array containing the regular bessel function of order i-1/2
      real(8) :: rjp(n+1)! Array containing the first derivative of the regular bessel function of order i-1/2

! !DEFINED PARAMETERS:
      real(8), parameter :: rtpio2 = 1.25331413731550025120788d+0
      


!
! !EXTERNAL ROUTINES: 
!


      external bessj


!
! !INTRINSIC ROUTINES: 
!

      
      intrinsic sqrt      

! !REVISION HISTORY:
!
! Original version: sphbes.for from numerical recipies.
!  (c) copr. 1986-92 numerical recipes software &124i..
!
! Modified: 30. Nov. 2004 by RGA: output vector with values for l=0...nl
!            1st. Dic. 2004 by RGA: values at x=0 included
!
!EOP
!BOC 
            
      if(n.lt.0.or.x.lt.0.)then
        write(6,*) "ERROR in sphbes:bad arguments"
        write(6,*)'  input: x = ',x
        write(6,*)'  input: n = ',n
        stop 'bad arguments in sphbes'
      endif   
      
      if(x.gt.0.0d0)then
        call bessj(x,n,rj,rjp)

        factor=rtpio2/sqrt(x)
        do l=1,n+1
          sj(l)=factor*rj(l)
          sjp(l)=factor*rjp(l)-sj(l)/(2.0d0*x)
        enddo  
      else
        do l=1,n+1
          sj(l)=0.0d0
          sjp(l)=0.0d0
        enddo 
        sj(1)=1.0d0
        sjp(2)=1.0d0/3.0d0 
      endif
      
      end subroutine sphbes
!EOC
