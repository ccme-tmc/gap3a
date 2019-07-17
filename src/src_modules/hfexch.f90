!BOI
!
! !MODULE: hfexch
      module hfexch
      implicit none 
!EOI
!BOP      
! !PUBLIC VARIABLES:
                                                   
      real(8)  :: exhf                       ! HF exchange energy
      complex(8), allocatable :: exq(:)      ! q-dependent exchange energy 
      complex(8) :: exq0s                    ! the singular part exq at q=0  

!EOP
!BOC
      integer,private::ierr
      contains
      
      subroutine init_hfexch(nirkp)
      integer(4) nirkp
      
      allocate(exq(nirkp),stat=ierr)
      if(ierr.ne.0) then 
        write(6,*) "initselfe: Fail to allocate memory"
        stop
      endif       

      end subroutine 

      subroutine end_hfexch
        deallocate(exq) 
      end subroutine

      end module 
!EOC
