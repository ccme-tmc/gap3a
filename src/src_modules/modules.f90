!BOP
!
! !MODULE: q0barc
      module q0barc

! !DESCRIPTION:
!
!This module declares the variables used in barcq0 and calcsing
!
! !PUBLIC VARIABLES:

        integer(4) :: nglen
        real(8), allocatable :: glen0(:)
        real(8), allocatable :: sing(:,:,:)
        complex(8), allocatable :: phase(:,:,:)
        
      end module q0barc
!EOP               
!---------------------------------------------------------------
!BOP
!



!------------------------------------------------------------------------
!BOP
!
! !MODULE: rotylm
      module rotylm
      
! !PUBLIC VARIABLES:
      
      complex(8), allocatable :: djmm(:,:) ! rotation matrix for ylm's
      
      end module rotylm
      
!EOP                        
!--------------------------------------------------------------
