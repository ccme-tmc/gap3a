      
!EOP
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: rspevec
      module rspevec
      
! !PUBLIC VARIABLES:

      integer(4) :: igmin,igmax,at1,at2,ibmin,ibmax
      
      integer(4) :: iik,jjk
      
      integer(4) :: ibmin2,ibmax2
      
      complex(8), allocatable :: evecmt(:) 
      
      complex(8), allocatable :: evecmts(:) 
      
      real(8), allocatable :: rg(:)
      
      real(8), dimension(3) :: rd
      
      end module rspevec
