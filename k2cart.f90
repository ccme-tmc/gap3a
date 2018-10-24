!BOP
!
! !ROUTINE: k2cart
!
! !INTERFACE:
      subroutine k2cart(ik,idiv,rk)

! !DESCRIPTION:
!
!Given the integer submesh coordinates of the kpoint, this subroutine
!returns its real cartesian coordinates
!
 
! !USES:

      use struk, only: pia,br2,ortho
      
! !INPUT PARAMETERS:

      implicit none
      integer(4), intent(in) :: ik(3) ! submesh coordinates of the k-point
      integer(4), intent(in) :: idiv  ! divisor of the submesh coords.
      
! !OUTPUT PARAMETERS:
      real(8), intent(out) :: rk(3)   ! real cartesian coords of the k-point

! !LOCAL VARIABLES:
      
      integer(4) :: i

 
! !REVISION HISTORY:
!
! Created: 24th. March 2004, by RGA
!
!EOP
!BOC

      if(ortho)then
        do i=1,3
          rk(i)=dble(ik(i))*pia(i)/dble(idiv)
        enddo     
      else
        do i=1,3
          rk(i)=( dble(ik(1))*br2(i,1)+dble(ik(2))*br2(i,2)             &
     &           +dble(ik(3))*br2(i,3))/dble(idiv)
        enddo
      endif
    
      end subroutine k2cart
!EOC
