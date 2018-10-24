!BOP
!
! !ROUTINE: w2k_writeklist
!
! !INTERFACE:
      subroutine w2k_writeklist(fid,nikp,ndiv,idiv,weight,klist)

! !DESCRIPTION:
!
! This subroutine writes the klist to the file case.klist, conforming to
! the WIEN2k standards
!
! !INPUT PARAMETERS:
      implicit none
      integer(4), intent(in) :: fid 
      integer(4), intent(in) :: nikp    ! Number of irreducible k-points
      integer(4), intent(in) :: ndiv(3) ! Number of divisions of the submesh in each direction.
      integer(4), intent(in) :: idiv    ! Common divisor of the integeral     kpoints
      integer(4), intent(in) :: weight(nikp) ! The geometrical weight of each k-point
      integer(4), intent(in) :: klist(3,nikp) !list of integer irreducible k-points

! !LOCAL VARIABLES:
    
      integer(4) :: ikp,j    
    
!EOP
!BOC
      write(fid,100) 1,(klist(j,1),j=1,3),idiv,dble(weight(1)),   &
     &      ndiv(1)*ndiv(2)*ndiv(3),ndiv
      do ikp =2,nikp
        write(fid,101) ikp,(klist(j,ikp),j=1,3),idiv,dble(weight(ikp))
      enddo
      write(fid,'("END")')
          
  100 format(i10,4i5,f5.1,5x,i5,' k, div: (',3i3,')')          
  101 format(i10,4i5,f5.1)
  
      end subroutine w2k_writeklist
!EOC
