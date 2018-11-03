!BOP
!
! !ROUTINE: calcmmatcc
!
! !INTERFACE:
      subroutine calcmmatcc(isp)
      
! !DESCRIPTION:
!
! This subroutine calculates the momentum matrix elements between pairs of
! core states.      
!
! !USES:

      use struk,      only: nat
      use constants, only: imag
      use core,      only: ncore, lcore
      use mommat,   only: iucl1ucl, iuclucl1, mmatcc
      implicit none
! !Input variables 
      integer(4) :: isp

      
! !LOCAL VARIABLES:
    
      
      integer(4) :: iat        ! Counter: runs over equivalent atmos
      
      integer(4) :: ic, jc     ! Counter: run over core states
      
      integer(4) :: l1, l2     ! angular momentum quantum numbers
      
      integer(4) :: iclm,jclm  ! run over core states including different
!                                m's

      integer(4) :: m1,m2      ! azimutal quantum number
      
      integer(4) :: l12        ! difference between l1 and l2
      
      integer(4) :: m12        ! difference between m1 and m2
      
      real(8) :: mxmy          ! matrix element of the (dx-idy) operator
      
      real(8) :: mxpy          ! matrix element of the (dx+idy) operator
      
      real(8) :: mz            ! matrix element of the (dz) operator
      
      real(8) :: num,den,flm   ! numerator, denominator and final value of
!                                the corresponding F^(i)_{lm} coefficient
      

! !INTRINSIC ROUTINES: 

 
      intrinsic iabs
      intrinsic isign
      intrinsic sqrt

! !REVISION HISTORY:
!
! Created 25. Oct. 2004 by RGA
! Last modified 27. Oct. 2004 by RGA
!
!EOP
!BOC
!
!     Loop over equivalent atoms
!
      do iat=1,nat
!
!       initialize the core states index
!
        iclm=0
!        
!       Loop over core states
!
        do ic=1, ncore(iat)
!
!         Calculate the angular momentum quantum number of the core state
!        
          l1= lcore(ic,iat) 
!          
!         Loop over z-component of the angular momentum
!
          do m1=-l1,l1
            iclm=iclm+1
            jclm=0
!        
!           Loop over core states
!
            do jc=1, ncore(iat)
!
!             Calculate the angular momentum quantum number of the core state
!        
              l2= lcore(jc,iat) 
!          
!             Loop over z-component of the angular momentum
!
              do m2=-l2,l2
                jclm=jclm+1
                l12=l1-l2
                m12=m1-m2
!
!               Initialize the operators matrix elements to zero
!                
                mxmy=0.0d0  
                mz=0.0d0  
                mxpy=0.0d0  

                select case (l12)

                case (1)   ! l1 = l2 + 1

                  den=dble((2*l2+1)*(2*l2+3))

                  select case (m12)

                  case (1)  ! m1 = m2 + 1

                    num=dble((l2+m2+1)*(l2+m2+2))
                    flm=-sqrt(num/den)
                    mxpy=flm*iucl1ucl(iat,ic,jc,isp)
                    
                  case (0)  ! m1 = m2 

                    num=dble((l2-m2+1)*(l2+m2+1))
                    flm=sqrt(num/den)
                    mz=flm*iucl1ucl(iat,ic,jc,isp)
                  
                  case (-1) ! m1 = m2 - 1

                    num=dble((l2-m2+1)*(l2-m2+2))
                    flm=sqrt(num/den)
                    mxmy=flm*iucl1ucl(iat,ic,jc,isp)
                  
                  end select ! m12
                  
                case (-1)  ! l1 = l2 - 1    

                  den=dble((2*l2-1)*(2*l2+1))
                  
                  select case (m12)  

                  case (1)  ! m1 = m2 + 1

                    num=dble((l2-m2)*(l2-m2-1))
                    flm=sqrt(num/den)
                    mxpy=flm*iuclucl1(iat,ic,jc,isp)
                    
                  case (0)  ! m1 = m2

                    num=dble((l2-m2)*(l2+m2))
                    flm=sqrt(num/den)
                    mz=flm*iuclucl1(iat,ic,jc,isp)
                  
                  case (-1) ! m1 = m2 - 1

                    num=dble((l2+m2)*(l2+m2-1))
                    flm=-sqrt(num/den)
                    mxmy=flm*iuclucl1(iat,ic,jc,isp)
                  
                  end select ! m12

                end select ! l12
!
!               determine the momentum matrix elements
!
                mmatcc(1,iclm,jclm,iat,isp)=-5.0d-1*imag*(mxpy+mxmy)
                mmatcc(2,iclm,jclm,iat,isp)=5.0d-1*(mxmy-mxpy)
                mmatcc(3,iclm,jclm,iat,isp)=-imag*mz

              enddo ! m2
            enddo ! jc
          enddo ! m1
        enddo ! ic
      enddo ! iat
!      
!     output formats
!
!   10 format(6i4,3g18.10)
!   12 format(6i4,2g18.10)
   11 format(6i4,3('(',g18.10,',',g18.10,') '))
!   20 format(4g18.10)
!   21 format(8('(',g18.10,',',g18.10,') '))
!   30 format(2i4,3g18.10)                
 
      end subroutine calcmmatcc
!EOC                
                
