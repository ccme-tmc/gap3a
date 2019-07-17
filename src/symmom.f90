!BOP
!
! !ROUTINE: symmom
!
! !INTERFACE:
      subroutine symmom(p,psq)
!
! !DESCRIPTION:
!
! This subroutine symetrizes one momentum matrix elements according to the
! symmetries of the lattice      
!
! !USES:

      use constants, only: czero
      use struk,     only: nsym, imat
      
! !INPUT/OUTPUT PARAMETERS:

      implicit none
      
      complex(8), intent(in) :: p(3)
      
      real(8), intent(out)   :: psq(3)
      
! !LOCAL VARIABLES:
      
      integer(4) :: isym, i, j
      
      real(8) :: normfac
      real(8) :: ps2(3)
      
      complex(8) :: p1(3)
 
! !REVISION HISTORY:
!
! Created 07.09.05 by RGA
!
!EOP
!BOC
      ps2(1:3)=0.0d0
      normfac=1.0d0/dble(nsym)
      do isym=1,nsym
        p1(1:3)=czero      
        do i=1,3
          do j=1,3
            p1(i)=p1(i)+dble(imat(j,i,isym))*p(j)
          enddo
        enddo
        do i=1,3
          ps2(i)=ps2(i)+real(p1(i)*conjg(p1(i)))
        enddo  
      enddo
      do i=1,3
        psq(i)=ps2(i)*normfac
      enddo
      
      end subroutine symmom
!EOC               
