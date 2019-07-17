!BOP
!
! !ROUTINE: intstipw
!
! !INTERFACE: 
      subroutine intstipw
      
! !DESCRIPTION:
!
! This subroutine calculates the integral over the intestitial region
! between a star and a planewave
!
! !USES:

      use constants, only: czero, imag, pi
      use xcpot,     only: ksxc, nksxc, istpw
      use struk,     only: nsym, imat, tau
      use recipvec,  only: gindex, npw2apw
      
! !LOCAL VARIABLES:

      implicit none
      
      integer(4) :: i         !(Counter) Runs over cartesian coordinates.
      integer(4) :: iks       !(Counter) Runs over stars
      integer(4) :: ipw       !(Counter) Runs over G-vectors.
      integer(4) :: isym      !(Counter) Runs over symmetries
      integer(4) :: j         !(Counter) Runs over cartesian coordinates.
      
      integer(4), dimension(3) :: ig ! Integer coords. of G-G'
      integer(4), dimension(3) :: ktemp ! Temporary storage of the
!                                          G-vector
      
      
      real(8)    :: tpi       ! 2*\pi
      real(8)    :: phaset    ! Phase of each term summing up to the
!                               coefficient phi_m

      
      complex(8) :: phim
      complex(8) :: integ
      
!
! !REVISION HISTORY:
!
! Created 18.08.2005 by RGA
!
!EOP
!BOC
!
!     Set constants
!
      tpi=2.0d0*pi
!      
!     Loop over stars
!
      do iks=1,nksxc 
        istpw(1:npw2apw,iks)=czero
!         
!       loop over symmetry operations
!        
        do isym=1,nsym
!        
!         Calculate the phase G.tau and the rotated vector R.G
!          
          phaset=0.0d0
          do i=1,3
            phaset=phaset + tau(i,isym)*ksxc(i,iks)*tpi
            ktemp(i)=0
            do j=1,3
              ktemp(i)=ktemp(i)+imat(i,j,isym)*ksxc(j,iks)
            enddo ! j
          enddo ! i
          phim=exp(-imag*phaset)
!
!         Check if the rotated vector is already included
! 
          do ipw=1,npw2apw
            ig(1:3)=ktemp(1:3)-gindex(:,ipw)
            call int1ipw(ig,integ)
            istpw(ipw,iks)=istpw(ipw,iks)+phim*integ
          enddo ! ipw
        enddo ! isym

        do ipw=1,npw2apw
          istpw(ipw,iks) = istpw(ipw,iks)/nsym
        enddo ! ipw
      enddo ! iks  
      
      return
      
      end subroutine intstipw      
!EOC            
