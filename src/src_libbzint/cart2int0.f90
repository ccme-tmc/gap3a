!BOP
! 
! !ROUTINE: cart2int0
!
! !INTERFACE:
      subroutine cart2int0(nkp,rbas,aaa,klist,idvk)

! !DESCRIPTION:
!
! This subroutine transform the submesh coordinates of the kpoints into
! cartesian coordinates in the reciprocal space
!
! !INPUT PARAMETERS:
      implicit none
      integer, intent(in) :: nkp,idvk
      real(8), intent(in) :: rbas(3,3)  ! Basis vectors of the direct lattice
      real(8), intent(in) :: aaa(3)  
      
! !INPUT/OUTPUT PARAMETERS:
      integer, intent(inout) :: klist(3,nkp)! integer coordinates of  the kpoints

! !LOCAL VARIABLES:
      integer :: kpi,i,j,idiv,klst(3,nkp) 
      real(8) :: ak(3) 

!EOP      
!
!BOC
      do kpi=1,nkp
        do i=1,3
          ak(i)=0.0d0
          do j=1,3
            ak(i)=ak(i)+rbas(j,i)*dble(klist(j,kpi))
          enddo
        enddo
        do i=1,3
          klst(i,kpi)=nint(ak(i)/aaa(i))
        enddo
      enddo
      call divisi(nkp,idvk,klst)
      klist(1:3,1:nkp)=klst(1:3,1:nkp)

      end subroutine cart2int0

!EOC
