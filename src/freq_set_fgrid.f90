!BOP
!
! !ROUTINE: freq_set_fgrid
!
! !INTERFACE:
      subroutine freq_set_fgrid
      
! !DESCRIPTION:
!
! This subroutine generates the frequecy mesh for integration
!
! !USES:
      use constants, only: pi 
      use freq

! !LOCAL VARIABLES:

      implicit none
      
      integer :: i,n
      real(8) :: w0
      real(8), allocatable :: u(:),wu(:)
      character(len=80)::fginfo(4)
      logical::lprt=.false.

! !INTRINSIC ROUTINES: 
      intrinsic exp
! !EXTERNAL ROUTINES: 
      external gaulag
      external gauleg

! !REVISION HISTORY:
!
! Created 20.06.05 by RGA.
!
!EOP
!BOC
      call linmsg(6,'-',"Setup freq mesh")

      fginfo(1)='Equally spaced mesh'
      fginfo(2)='Gauss-Laguerre quadrature'
      fginfo(3)='Double Gauss-Legendre quadrature '
      fginfo(4)='Fermion Matsubara frequencies '

      allocate(omega(1:nomeg),womeg(1:nomeg))

      if(nomeg.eq.1) then 
        omega(1)=omegmax
        womeg(1)=1.d0
        return 
      endif  
      
      write(6,*)'Frequency grid'
      write(6,101) iop_fgrid
      write(6,102) fginfo(iop_fgrid)
      select case(iop_fgrid) 

      case(1)    ! Equaly spaced mesh c calcultions, not for integration)

        do i=1,nomeg
          omega(i)= omegmin + dble(i-1)*(omegmax-omegmin) /dble(nomeg-1)
          womeg(i)= (omegmax-omegmin) /dble(nomeg-1)
        enddo  

      case(2)    ! Grid for Gauss-Laguerre quadrature

        allocate(wu(nomeg))
        call gaulag(omega(1:nomeg),wu(1:nomeg),nomeg,-1.0d0)
        do i = 1, nomeg
          womeg(i) = wu(i)*exp(omega(i))*omega(i)
        enddo  
        deallocate(wu)

      case(3)    ! Double Gauss-Legendre quadrature from 0 to omegmax and from omegmax to infinity 
        n=nomeg/2
        nomeg=n*2
        allocate(u(n),wu(n))
        call gauleg(-1.0d0,1.0d0,u,wu,n)
        do i = 1, n
          omega(i)=omegmax*(u(i)+1.0d0)/2.0d0
          omega(2*n-i+1)=2.0d0*omegmax/(u(i)+1.0d0)
          womeg(i)=omegmax*wu(i)/2.0d0
          womeg(2*n-i+1)=2.0d0*omegmax*wu(i)/((u(i)+1.0d0)*(u(i)+1.0d0))
        enddo ! i  
        deallocate(u,wu)

      case(4) ! Fermion Matsubara frequencies 
        do i=1,nomeg 
          omega(i) = omegmin*(2*i-1) 
          womeg(i) = 1.0 
        enddo 
      end select

      do i=1,nomeg
        write(6,104) i,omega(i),womeg(i) 
      enddo 

  101 format(' iop_fgrid=',i2,a)
  102 format(2x,a)
  103 format('number of frequencies:',i4,' omegmax=',f18.10)
  104 format(i4,2f12.4)
      return

      end subroutine freq_set_fgrid
!EOC      
            
            
