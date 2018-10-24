!BOP
!
! !ROUTINE: task_chkbz
!
! !INTERFACE:
      subroutine task_chkbz
      
! !DESCRIPTION:
!
! This subroutine calculate frequency-dependent U from the constrained
! RPA approach 
!
! !USES:
      use bands,    only: nspin
      use freq,     only: nomeg
      use kpoints,  only: nirkp,idikp,nkp,kqid,kpirind
      
! !LOCAL VARIABLES:
      implicit none
      integer(4) :: iq       ! (Counter) Runs over q-points
      integer:: iom 
      integer(4) :: ierr 

! !INTRINSIC ROUTINES: 
      intrinsic cpu_time      

! !REVISION HISTORY:
!      
!EOP
!BOC     

      do iq=1,nkp   !! Loop over q-points.
        write(6,*) "q-dependent weights iq=",iq
        call bz_calcqdepw(iq)  !! Calc the q-dependent integration weights
      enddo ! iqa

  100 format(10f12.6) 

      end subroutine task_chkbz
!EOC      
