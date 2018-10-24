!BOP
!
! !ROUTINE: writesym
!
! !INTERFACE:
      subroutine writesym(nsym,iio)

! !INPUT PARAMETERS:      
      implicit none
      
      integer(4), intent(in) :: nsym       ! Number of symmetry operations
      
      integer(4), intent(in) :: iio(3,3,*) ! Symmetry operation matrices

! !LOCAL VARIABLES:      
      
      integer(4) :: i,i1,i2,i3,i4,j,k     ! Several counters
      
! !REVISION HISTORY:
!
! Created: Feb. 2004 by RGA
! Last Modified: 8th. March 2004 by RGA
!
!EOP
!BOC
      do i=1,int(float(nsym)/4.+.9)                                 
        i1=4*i-3                                                          
        i2=4*i-2                                                          
        i3=4*i-1                                                          
        i4=4*i                                                            
        write(27,120) i1,i2,i3,i4                                        
        do j=1,3                                                      
          write(27,140) (iio(j,k,i1),k=1,3),(iio(j,k,i2),k=1,3),&
     &                  (iio(j,k,i3),k=1,3),(iio(j,k,i4),k=1,3)
        enddo
      enddo
  120 format (t5,'SYMMETRY MATRIX NR.',i3,t30,'SYMMETRY MATRIX NR.' &    
     & ,i3,t55,'SYMMETRY MATRIX NR.',i3,t80,'SYMMETRY MATRIX NR.',i3)   
  140 format (t5,3i5,t30,3i5,t55,3i5,t80,3i5)                           
      end subroutine writesym
!EOC
