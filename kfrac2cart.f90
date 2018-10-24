!BOP
!
! !ROUTINE: kfrac2cart
!
! !INTERFACE:
      subroutine kfrac2cart(xk,rk)

! !DESCRIPTION:
!
!Given the integer submesh coordinates of the kpoint, this subroutine
!returns its real cartesian coordinates
!
 
! !USES:

      use struk, only: pia,br2,ortho
      
! !INPUT PARAMETERS:

      implicit none
      
! !OUTPUT PARAMETERS:
      real(8), intent(in) ::  xk(3)  ! the k-point in the fractioal coordinates
      real(8), intent(out) :: rk(3)  ! real cartesian coords of the k-point

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
          rk(i)=xk(i)*pia(i)
        enddo     
      else
        do i=1,3
          rk(i)= xk(1)*br2(i,1)+xk(2)*br2(i,2)+ xk(3)*br2(i,3) 
        enddo
      endif
    
      end subroutine kfrac2cart
!EOC
