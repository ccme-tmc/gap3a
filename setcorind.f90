!BOP
!
! !ROUTINE: setcorind
!
! !INTERFACE:
      subroutine setcorind
!
! !DESCRIPTION:
!
! This subroutine sets a unique index for all the core states of all atoms
!
! !USES:

      use core,  only: corind,clmind, lcore, ncg, nclmmax, ncore,nclm
      use struk, only: nat, ndf, mult

! !LOCAL VARIABLES:

      implicit none
      
      integer(4) :: iat     ! (Counter) Runs over inequivalent atoms
      integer(4) :: ieq     ! (Counter) runs over equivalent atoms.
      integer(4) :: ic      ! (Counter) runs over core states of an atom
      integer(4) :: icg     ! (Counter) runs over all core states
      integer(4) :: iclm    ! index for core states at one atom 
      integer(4) :: idf  
      integer(4) :: lc      ! Angular momentum of the core state
      integer(4) :: mc
      
! !REVISION HISTORY:
!      
! Created 23.08.05 by RGA
!
!EOP
!BOC
      allocate(corind(5,nclmmax*ndf),clmind(nclmmax,ndf),nclm(ndf))
      idf=0
      icg=0
      do iat=1,nat
        do ieq=1,mult(iat)
          idf=idf+1
          iclm = 0
          do ic = 1, ncore(iat)
            lc = lcore(ic,iat) 
            do mc=-lc,lc
              iclm=iclm+1
              icg=icg+1
              clmind(iclm,idf) = icg 
              
              corind(1,icg)=iat
              corind(2,icg)=idf
              corind(3,icg)=ic
              corind(4,icg)=lc
              corind(5,icg)=mc
            enddo
          enddo
          nclm(idf) = iclm 
        enddo
      enddo  
      ncg=icg
      return
      end subroutine setcorind
!EOC            
              
 
      
