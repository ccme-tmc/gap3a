!BOP
!
! !ROUTINE: tetiwsurf
!
! !INTERFACE:
       subroutine tetiwsurf(nik,nb,ebd,omeg,iwsurf)
!     
! !DESCRIPTION:
!
!   This subroutine calculates the weight on one k-point of a certain band for 
! a normal surface integration which is not q-dependent. 

! !USES:
       use bzint,  only: nirtet, tndi, wirtet, tvol 
       use tetra_internal
       
       implicit none      
       
! !INPUT PARAMETERS:
 
       integer, intent(in) :: nik        ! Number of irreducible k-points
       integer, intent(in) :: nb         ! Number of bands
       real(8), target, intent(in) :: ebd(nb,nik)  ! Band energies
       real(8), intent(in) :: omeg       ! fermi energy
       
! !OUTPUT PARAMETERS:
       
       real(8), intent(out)   :: iwsurf(nb,nik)    ! the value of the weight
       
       
! !INTRINSIC ROUTINES:
       
       intrinsic size
       
! !EXTERNAL ROUTINES:
       
       external intwsurf
       
! !REVISION HISTORY:
!
!   Created: 10th. Jan 2005 by XZL
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

!
!     intwsurf is analogous to intw, difference is that it is for surface integration
!
      call intwsurf(omeg,iwsurf)

      end subroutine tetiwsurf
      
!EOC
