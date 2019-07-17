!BOP
!
! !ROUTINE: edifwt
!
! !INTERFACE:

      subroutine stweight_numeric(iop_omg,de_vert,omg,wt_vert)
!
! !DESCRIPTION:
!
! This subroutine calculates by numerical integration the weight 
! on the whole small tetrahedron 
! in which the $k$ states are fully occupied and $k-q$ states are fully 
! unoccupied. This is for the $sigfreq=3$ case when we consider the 
! imaginary frequency. 
!                                             
! !USES:
      use tetra_internal, only: weighttol, weightwarn,fout,eta_refreq,& 
     &                          vol_small_tetra,n_gauq,x_gauq,w_gauq  

! !INPUT PARAMETERS:
      implicit none
     
      integer, intent(in) :: iop_omg       ! 0 - real frequency 
                                           ! 1 - imaginary frequency  
      real(8), intent(in) :: de_vert(4)    ! difference of the energy 
                                               ! in k-mesh tetrahedron vertices 
                                               ! and k-q mesh tetrahedron vertices.
     

      real(8), intent(in) :: omg              ! the frequency omega to be calculated

! !OUTPUT PARAMETERS:            
      real(8), intent(out) :: wt_vert(4)   ! the weight on the whole tetrahedron.

! !LOCAL VARIABLES:

      integer(4) :: j,k
      integer(4) :: ivert

      real(8)    :: d0,d10,d20,d30
      real(8)    :: x,y,z         !! common variables used by different internal subroutines

      character(20):: sname="stweight_numeric"

      integer:: nmax = 10, nmin = 3 
      real(8):: eps = 1.e-4
      logical:: l_conv

! 
! !SYSTEM ROUTINES:

      intrinsic datan
      intrinsic dlog
      intrinsic dsign
 
! !REVISION HISTORY:
!
! Created 05.11.2004 by XZL.
!

!EOP
!BOC
      wt_vert(1:4) = 0.0d0
      d0 = de_vert(1) 
      d10 = de_vert(2)-de_vert(1)
      d20 = de_vert(3)-de_vert(1)
      d30 = de_vert(4)-de_vert(1)

      do ivert=1,4
        l_conv = .true.
        call gauq_z(1.d0, wt_vert(ivert)) 
        if(.not.l_conv) then 
          write(fout,1) sname,ivert,omg,vol_small_tetra,     &
     &              de_vert,wt_vert(ivert)
        endif 
      enddo  

 1    format(a,'- WARNING: numerical integration not converged !',/ &
     &' vertix:',i4,' omeg =',g16.6,' vol_tetra=',g16.6,/,&
     &' de_vert =',4g16.6,/,&
     &' weight=',g16.6)

      contains

!
!     This defines the integrand 
!
      real(8) function func(x,y,z)
      real(8):: x,y,z

      real(8):: de,comm,fv(4) 
      
      de = d0 + x*d10 + y*d20 + z*d30
 
      if(iop_omg.eq.0) then 
        comm = (omg-de)/((omg-de)**2+eta_refreq**2)
      else
        comm = -2.d0*de/(omg*omg + de*de)
      endif 

      fv(1) = (1.d0-x-y-z)*comm
      fv(2) = x*comm
      fv(3) = y*comm
      fv(4) = z*comm
      
      func = fv(ivert) 
      end function 

      real(8) function f(xx) 
      real(8):: xx
      x = xx
      f = func(x,y,z)  
      end function

      real(8) function g(yy) 
      real(8):: yy
      y = yy
      call gauq_x(1.0-y-z,g)
      end function  

      real(8) function h(zz) 
      real(8):: zz
      z = zz 
      call gauq_y(1.0-z,h)
      end function  

      subroutine gauq_x(a,s)
      real(8),intent(in) :: a
      real(8),intent(out):: s
      integer :: i
      s=0.d0
      do i=1,n_gauq
        s = s + w_gauq(i)*a*f(a*x_gauq(i))
      enddo 
      end subroutine
 
      subroutine gauq_y(a,s)
      real(8),intent(in) :: a
      real(8),intent(out):: s
      integer :: i
      s=0.d0
      do i=1,n_gauq
        s = s + w_gauq(i)*a*g(a*x_gauq(i))
      enddo
      end subroutine

      subroutine gauq_z(a,s)
      real(8),intent(in) :: a
      real(8),intent(out):: s
      integer :: i
      s=0.d0
      do i=1,n_gauq
        s = s + w_gauq(i)*a*h(a*x_gauq(i))
      enddo
      end subroutine

      subroutine qint1d_x(a,b,s)
      real(8),intent(in) :: a,b
      real(8),intent(out):: s

      integer :: n,it,tnm,j
      real(8) :: olds,del,x,ss

      olds = -1.e30
      do n=1,nmax
        if(n.eq.1) then 
          s = 0.5*(b-a)*(f(a)+f(b))
        else 
          it = 2**(n-2) 
          tnm = it 
          del = (b-a)/tnm 
          x = a+0.5*del 
          ss = 0.0 
          do j=1,it
            ss = ss + f(x) 
            x = x + del 
          enddo 
          s = 0.5*(s+(b-a)*ss/tnm)
        endif 
        if(n.gt.nmin) then 
           if(abs(s-olds).lt.eps*abs(olds).or.&
     &        (s.eq.0.0.and.olds.eq.0.0)) return 
        endif  
        olds = s
      enddo 
      l_conv = .false.
      end subroutine



      subroutine qint1d_y(a,b,s)
      real(8),intent(in) :: a,b
      real(8),intent(out):: s

      integer :: n,it,tnm,j
      real(8) :: olds,del,x,ss

      olds = -1.e30
      do n=1,nmax
        if(n.eq.1) then
          s = 0.5*(b-a)*(g(a)+g(b))
        else
          it = 2**(n-2)
          tnm = it
          del = (b-a)/tnm
          x = a+0.5*del
          ss = 0.0
          do j=1,it
            ss = ss + g(x)
            x = x + del
          enddo
          s = 0.5*(s+(b-a)*ss/tnm)
        endif
        if(n.gt.nmin) then
           if(abs(s-olds).lt.eps*abs(olds).or.&
     &        (s.eq.0.0.and.olds.eq.0.0)) return
        endif
        olds = s
      enddo
      l_conv = .false.
      end subroutine

      subroutine qint1d_z(a,b,s)
      real(8),intent(in) :: a, b
      real(8),intent(out):: s

      integer :: n,it,tnm,j
      real(8) :: olds,del,x,ss

      olds = -1.e30
      do n=1,nmax
        if(n.eq.1) then
          s = 0.5*(b-a)*(h(a)+h(b))
        else
          it = 2**(n-2)
          tnm = it
          del = (b-a)/tnm
          x = a+0.5*del
          ss = 0.0
          do j=1,it
            ss = ss + h(x)
            x = x + del
          enddo
          s = 0.5*(s+(b-a)*ss/tnm)
        endif
        if(n.gt.nmin) then
           if(abs(s-olds).lt.eps*abs(olds).or.&
     &        (s.eq.0.0.and.olds.eq.0.0)) return
        endif
        olds = s
      enddo
      l_conv = .false.
      end subroutine

      end subroutine stweight_numeric
!EOC      
