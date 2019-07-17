!BOP
!
! !ROUTINE: task_coul
!
! !INTERFACE:
      subroutine task_coul

! !DESCRIPTION:
!
! This subroutine performs the test of the bare coulomb potential,
! calculates it with the mixed basis and planewaves and compares both.
!
! !USES:
 
      use mixbasis,   only: matsiz,mpwmix, wi0, mpwipw,  &
     &                      init_mixbasis,end_mixbasis
      use recipvec,   only: ngqbarc, ngq, ig0
      use barcoul
      use liboct_parser
      
      
! !LOCAL VARIABLES:
      
      implicit none
      
      integer(4) :: iq,ierr,iq1,iq2
            
      character(len=15) :: fnint

! !REVISION HISTORY:
!
! Created 16.09.2005 by RGA
!
!EOP
!BOC

      ierr=loct_parse_isdef('coul')
      if(ierr.eq.1) then 
        call loct_parse_block_int('coul',0,0,iq1)
        call loct_parse_block_int('coul',0,1,iq2)
      else 
        iq1=1
        iq2=1
      endif 

      do iq=iq1,iq2
        call init_mixbasis(iq) 
        call coul_mpwipw(iq)
        allocate(mpwmix(matsiz,ngqbarc(iq)))
        if(iq.eq.1) then 
          call coul_wmix0
        endif 
        call coul_mpwmix(iq)

        call testbarc(iq)
        deallocate(mpwmix)
        call end_mixbasis
      enddo  

      return
      end subroutine task_coul

!EOC      
