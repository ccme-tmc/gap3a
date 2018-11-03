! BOP

! 
! !ROUTINE: trilinear
!
! !INTERFACE:

      subroutine trilinear(nods,coordina,resul)

! !DESCRIPTION:
! This subroutine use the trilinear method to get the interpolated energy.

! !INPUT PARAMETERS:

      real(8), intent(in) :: nods(8) ! energy on the vertices

      real(8), intent(in) :: coordina(3) ! the coordinates of the point to
!                                          be drawn

! !OUTPUT PARAMETERS:

      real(8), intent(out) :: resul

! !LOCAL VARIABLES:

      real(8) :: x,y,z,resul0

! !REVISION HISTORY:
! created on July 1st, 2004 by XZL
!EOP

!BOC
      x=coordina(1)
      y=coordina(2)
      z=coordina(3)

      resul0=0.0d0
!      print*, 'nods:', nods(:)
!      print*, 'coordinates', coordina(:) 

!      resul0=nods(1)+(nods(5)-nods(1))*x+ &
!     &       (nods(4)-nods(1))*y+(nods(2)-nods(1))*z

      resul0=nods(1)*(1-x)*(1-y)*(1-z)+  &
     &       nods(2)*(1-x)*(1-y)*z+      &
     &       nods(3)*(1-x)*y*z+          &
     &       nods(4)*(1-x)*y*(1-z)+      &
     &       nods(5)*x*(1-y)*(1-z)+      &
     &       nods(6)*x*(1-y)*z+          &
     &       nods(7)*x*y*z+              &
     &       nods(8)*x*y*(1-z)

      resul=resul0
!      print*, 'result: ', resul
      return

      end subroutine trilinear

!EOC
