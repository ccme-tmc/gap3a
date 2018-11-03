!BOP
!
! !ROUTINE: intw 
!
! !INTERFACE:
       subroutine intw(efer,kwt)
!
! !DESCRIPTION:
!
!   This subroutine calculates the integration weight of each k-point for
!   each band

! !USES:
       
       use order
       use tetra_internal
       use bzint, only: iop_bcor
       
       implicit none     

! !INPUT PARAMETERS:
       real(8), intent(in)  :: efer                 ! fermi energy
       
! !OUTPUT PARAMETERS:
       real(8), intent(out) :: kwt(nband,nirkp) ! the weight of each k-point for each band

!
! !LOCAL VARIABLES:
 
       integer(4) :: it,i,ib,kin
       integer(4), dimension(4) :: ik
       real(8) :: term
       real(8), dimension(4) :: ee
       real(8), dimension(4) :: w1t
       real(8), dimension(4) :: wcor
       external  intweight1t

! !REVISION HISTORY:
!
!   Created:  3th. March 2004. by RGA 
! 
!EOP
!BOC
      kwt=0.d0 
      do it=1,ntet
        do ib=1,nband 
          do i=1,4
            ee(i)=eband(ib,tetcorn(i,it))
          enddo
          call sort(4,ee,ik)
          call intweight1t(ee,efer,vt,w1t)

          if(iop_bcor.eq.1) then 
            call bloechlcor(ee,efer,vt,wcor)
          else
            wcor = 0.d0
          endif 

          do i=1,4
            term=(w1t(i)+wcor(i))*tetweig(it)
            kin= tetcorn(ik(i),it)
            kwt(ib,kin)=kwt(ib,kin)+term
          enddo
        enddo
      enddo
      return
      
      end subroutine intw
      
!EOC
