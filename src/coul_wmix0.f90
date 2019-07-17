!BOP
!
! !ROUTINE: coul_wmix0
!
! !INTERFACE:
      subroutine coul_wmix0 

! !DESCRIPTION:
!
!This subroutine calculate the matrix elements $\mathcal{W}^i_0$
!

! !USES:

      use constants,  only: czero,cone,pi
      use mixbasis,   only: nmix,bigl,locmatsiz,mpwipw,wi0,rtl,  &
     &                      mbsiz
      use recipvec,   only: ngq
      use struk,      only: mult, vi,nat

! !LOCAL VARIABLES:

      implicit none
      
      integer :: mixind  ! Counter, runs over mixed basis functions
      integer :: idf     ! Counter, runs over all atoms
      integer :: iat     ! Counter, runs over inequivalent atoms
      integer :: ieq     ! Counter, runs over equivalent atoms
      integer :: irm    ! Counter, runs over radial mixed functions
      integer :: l1      ! Angular momentum quantum number of the mixed function
      integer :: ippw    ! Counter, runs over interstitial mixed basis functions
      complex(8),allocatable::work(:)


! !INTRINSIC ROUTINES: 
      intrinsic sqrt

! !REVISION HISTORY:
!
! Created 11.02.05 by RGA
!
!EOP
!BOC
      wi0(1:mbsiz)=czero
      mixind=0
      idf=0
      do iat=1,nat
        do ieq=1,mult(iat)
          idf=idf+1
          do irm=1,nmix(iat)
            l1=bigl(irm,iat)
            if(l1.eq.0)then
              mixind=mixind+1
              wi0(mixind)=sqrt(4.0d0*pi*vi)*rtl(irm,iat)
            else
              mixind=mixind+2*l1+1
            endif
          enddo ! irm
        enddo ! ieq
      enddo ! iat
      
      do ippw=1,ngq(1)
        mixind=locmatsiz+ippw
        wi0(mixind)=mpwipw(ippw,1)
      enddo

      end subroutine coul_wmix0
!EOC                
