!BOP
!
! !ROUTINE: tetiw 
!
! !INTERFACE:
       subroutine tetiw(nik,nb,ebd,efer,iw)

!     
! !DESCRIPTION:
!
!   This subroutine gives the weight of on one k-point for a certain band
! in an operator integration.
!  

! !USES:
       use bzint,  only: nirtet, tndi, wirtet, tvol 
       use tetra_internal
       
       implicit none      
       
! !INPUT PARAMETERS:
 
       integer, intent(in) :: nik        ! Number of irreducible k-points
       integer, intent(in) :: nb         ! Number of bands
       real(8), intent(in) :: ebd(nb,nik)  ! Band energies
       real(8), intent(in) :: efer       ! fermi energy
       real(8), intent(out)   :: iw(nb,nik)    ! the value of the integral
       target ebd
       
! !INTRINSIC ROUTINES:
       
       intrinsic size
       
! !EXTERNAL ROUTINES:

       external intw
       external average_degen_weights
       
! !REVISION HISTORY:
!
!   Created: 4th. March 2004 by RGA
!
!EOP
!BOC
 
      nirkp = nik
      ntet  = nirtet
      nband = nb
      vt = tvol 
      tetcorn => tndi(1:4,1:ntet)
      tetweig => wirtet(1:ntet)
      eband   => ebd(1:nband,1:nik)

      call intw(efer,iw)
    
!      call average_degen_weights(iw)
    
      end subroutine tetiw
      
!EOC
