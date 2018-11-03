!BOP
!
! !ROUTINE: setlocmixind
!
! !INTERFACE:
      subroutine setlocmixind
!
! !DESCRIPTION:
!
! This subroutine sets an array that stores the general index of the mixed
!function for a given mixed function of a given atom
!
! !USES:

      use mixbasis,  only: bigl, locmixind, lmixmax,nmix,nmixlm
      use struk,     only: nat,ndf, mult
      use task,      only: fid_outmb

! !LOCAL VARIABLES:

      implicit none
      
      integer(4) :: iat     ! (Counter) Runs over inequivalent atoms
      integer(4) :: ieq     ! (Counter) runs over equivalent atoms.
      integer(4) :: im
      integer(4) :: imix     ! (Counter) runs over all core states
      integer(4) :: irm
      integer(4) :: idf  
      integer(4) :: l      ! Angular momentum of the core state
      integer(4) :: m
      
! !REVISION HISTORY:
!      
! Created 23.08.05 by RGA
!
!EOP
!BOC
!     call boxmsg('-','Indexes of MT-sphere mixed basis functions')
      write(fid_outmb,*)'   chi_i=v_(aNL)Y_(LM)   (a = atom)'
      write(fid_outmb,101)
      allocate(locmixind(lmixmax,ndf),nmixlm(ndf))
      locmixind(:,:)=0
      nmixlm=0
      idf=0
      imix=0
      do iat=1,nat
        do ieq=1,mult(iat)
          idf=idf+1
          im=0
          do irm = 1, nmix(iat)
            l = bigl(irm,iat)
            do m=-l,l
              im=im+1
              imix=imix+1
              locmixind(im,idf)=imix
              write(fid_outmb,102)imix,idf,irm,l,m
            enddo
          enddo
          nmixlm(idf)=im 
        enddo
      enddo  
      write(fid_outmb,*)
  101 format(5x,'i',5x,'a',5x,'N',5x,'L',5x,'M')
  102 format(5i6) 
      return
      end subroutine setlocmixind
!EOC            
              
 
      
