!BOP
!
! !ROUTINE: bessj
!
! !INTERFACE:
      subroutine bessj(x,nl,rj,rjp)

! !DESCRIPTION:

!Returns the regular Bessel functions of first kind
!$\texttt{rj(l)}=J_{l-\frac{1}{2}}$ and their derivatives
!$\texttt{rjp(l)}=J'_{l-\frac{1}{2}}$, for positive \texttt{x} and for
!$1\le l \le nl+1$, using Steed's method.
!
! !INPUT PARAMETERS:

      implicit none
      
      real(8), intent(in) :: x
      
      integer(4), intent(in) :: nl
      
! !OUTPUT PARAMETERS:
      
      real(8), intent(out) :: rj(*)            

      real(8), intent(out) :: rjp(*)            
      
! !LOCAL VARIABLES:
      
      integer(4) :: i 
      integer(4) :: isign 
      integer(4) :: l 

      real(8) :: xnu
      real(8) :: b 
      real(8) :: c 
      real(8) :: d 
      real(8) :: del 
      real(8) :: f 
      real(8) :: fact 
      real(8) :: gam 
      real(8) :: h 
      real(8) :: p 
      real(8) :: q 
      real(8) :: rjl 
      real(8) :: rjl1 
      real(8) :: rjmu 
      real(8) :: rjp1 
      real(8) :: rjpl 
      real(8) :: rjtemp 
      real(8) :: w 
      real(8) :: xi 
      real(8) :: xi2 
      real(8) :: xmu 
      real(8) :: xmu2

! !DEFINED PARAMETERS:

      integer(4), parameter :: maxit=10000
      
      real(8), parameter :: eps=1.e-16

      real(8), parameter :: fpmin=1.e-30

      real(8), parameter :: pi=3.141592653589793d+0

! !REVISION HISTORY:
!
! Original subroutine: bessjy.for (c) copr. 1986-92 numerical recipes software &124i..
! Last modified: 30th. Nov. 2004 by RGA
!
!EOP
!BOC
      if(x.le.0..or.nl.lt.0) stop 'bad arguments in bessj'
      
      xnu=dble(nl)+0.5d0
      xmu=0.5d0
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
!
!     The Wronskian
!      
      w=xi2/pi
!
!     Evaluate the continued fraction expansion for J'_(nl-1/2)/J_(nl-1/2)
!     by the modified Lentz's method. isign keeps track of sign changes in
!     the denominator
! 
      isign=1
      h=xnu*xi
      if(h.lt.fpmin)h=fpmin
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,maxit
        b=b+xi2
        d=b-d
        if(abs(d).lt.fpmin)d=fpmin
        c=b-1.d0/c
        if(abs(c).lt.fpmin)c=fpmin
        d=1.d0/d
        del=c*d
        h=del*h
        if(d.lt.0.d0)isign=-isign
        if(abs(del-1.d0).lt.eps)goto 1
11    continue
      stop 'x too large in bessjy; try asymptotic expansion'
1     continue
!
!     Initialize J and J' for downward recurrence
!
      rjl=isign*fpmin
      rjpl=h*rjl
!      
!     Store values for later rescaling      
!      
      rjl1=rjl
      rjp1=rjpl
!
!     Downward recurrence (unnormalized)
!      
      fact=xnu*xi
      do 12 l=nl,0,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
12    continue
      if(rjl.eq.0.d0)rjl=eps
      f=rjpl/rjl
!
!     Equation 6.7.3 Numerical Recipies in Fortran
!      
      p=-.5d0*xi
      q=1.d0
!
!     Equation 6.7.6 Numerical Recipies in Fortran
!      
      gam=(p-f)/q
!
!     Equation 6.7.7 Numerical Recipies in Fortran
!      
      rjmu=sqrt(w/((p-f)*gam+q))
      rjmu=sign(rjmu,rjl)
!
!     Scale original J and J'
!      
      fact=rjmu/rjl
      rj(nl+1)=rjl1*fact
      rjp(nl+1)=rjp1*fact
      fact=xnu*xi
! 
!     Downward recurence
!      
      do l=nl,1,-1
        rj(l)=fact*rj(l+1)+rjp(l+1)
        fact=fact-xi
        rjp(l)=fact*rj(l)-rj(l+1)
      enddo ! l
      end subroutine bessj
      
!EOC
