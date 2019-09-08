      subroutine coul_barcq0_cutoff(icutoff)
!This subroutine calculates the matrix of the bare coulomb potential for
! q=0 and atomic functions with L=0

      use barcoul 
      use constants,  only: czero, pi
      use mixbasis,   only: bigl, maxbigl, mbl, nmix, nmixmax
      use q0barc,     only: nglen, glen0, sing, phase
      use struk,      only: mult, ndf, vi, pos, rmt,nat,rotij,rotloc


      implicit none

      integer, intent(in) :: icutoff ! the cut-off option. See barcoul

! !LOCAL VARIABLES:
      integer(4) :: iat       ! Indexes the atoms (inequivalent) 
      integer(4) :: idf       ! Indexes all the atoms 
      integer(4) :: ieq       ! Indexes the equivalent atoms 
      integer(4) :: igl       ! im(ipw) = locmatsiz + iipw)
      integer(4) :: im      ! Indexes the mixed wave function 
      integer(4) :: irm       ! Indexes the mixed basis functions  of atom iat
      integer(4) :: jat       ! Indexes the atoms (inequivalent) 
      integer(4) :: jdf       ! Indexes all the atoms
      integer(4) :: jeq       ! Indexes the equivalent atoms 
      integer(4) :: jm      ! Indexes the mixed wave function 
      integer(4) :: jrm       ! Indexes the mixed basis functions 
      integer(4) :: l1        ! Angular momentum quantum number of the atomic function irm
      integer(4) :: l2        ! Angular momentum quantum number of the  atomic function jrm
      integer(4) :: m1        ! z-component angular momentum quantum  number of the atomic function irm
      integer(4) :: m2        ! z-component angular momentum quantum number of the atomic function jrm
      
      real(8) :: qg1len       ! |qg1|
      real(8) :: invglen2
      real(8) :: len4
      real(8) :: t1, t2, t3
 
      complex(8) ::  sumipw
      integer(4) :: iprt 
 
! !EXTERNAL ROUTINES: 

      
      real(8), external :: getcgcoef

      complex(8), external :: getdjmm
      
      external calcjlam
      external k2cart
      external rotate
      external ylm
      external zgemm
!
! !REVISION HISTORY:
! 
!  Created 16th. March 2004 by RGA
! modified 31. March 2005 by RGA
! modified 08. Sept  2019 by MYZ
!
!EOP
!BOC
      iprt=0

      call coul_calcsing_cutoff
      idf=0
      im=0
      do iat = 1, nat
        do ieq = 1, mult(iat)
          idf=idf+1
          do irm = 1, nmix(iat)
            l1=bigl(irm,iat)
            if(l1.eq.0)then
              im = im + 1
              jm=0  
              jdf=0    
              do jat=1,nat
                do jeq=1,mult(jat)
                  jdf=jdf+1
                  do jrm = 1, nmix(jat)
                    l2=bigl(jrm,jat)
                    if(l2.eq.0)then
                      jm=jm+1
                      sumipw=czero
                      do igl = 2, nglen
                        qg1len = glen0(igl)
                        len4=qg1len*qg1len*qg1len*qg1len
                        invglen2=1.0d0/len4
                        sumipw=sumipw+phase(idf,jdf,igl)*               &
     &                         cmplx(sing(iat,irm,igl)*sing(jat,jrm,igl)&
     &                        *invglen2,0.0d0,8)
                      enddo ! igl  
                      vmat(im,jm)=16.0d0*pi*pi*vi*sumipw
                    else  
                      !do m2=-l2,l2
                      !  jm=jm+1
                      !enddo ! m2
                      jm=jm+2*l2+1
                    endif  
                  enddo ! jrm
                enddo ! jeq
              enddo ! jat
            else  
              !do m1=-l1,l1
              !  im = im + 1
              !enddo ! m1
              im = im +2*l1+1
            endif    
          enddo ! irm
        enddo ! ieq
      enddo ! iat                
      deallocate(sing)
      deallocate(phase)
      deallocate(glen0)
      return
  
      end subroutine coul_barcq0_cutoff
!EOC      
