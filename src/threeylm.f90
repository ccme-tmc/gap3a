!BOP
! !ROUTINE: threeylm
!
! !INTERFACE:
      function threeylm(l1,l2,l3,m1,m2,m3)

! !DESCRIPTION:
!
!This function calculate the Gaunt Coefficients (three spherical harmonics
!integral):
!
!\begin{equation}
!\int{Y_{l_1m_1}(\hat{r})Y_{l_2m_2}(\hat{r})Y^*_{l_3m_3}
!(\hat{r})d\hat{r}}
!\end{equation}
!
!\noindent
!numerically, using the grid points and weights according to the
!quadrature developed by Levedev and Laikov. 
!
! !USES:

      use angintegrals, only: n_angrid, sphar, w_angrid
      use constants,    only: czero, pi

! !INPUT PARAMETERS:     
      
      implicit none
      
      integer(4) :: l1 ! l index of the first spherical harmonic.
!
      integer(4) :: l2 ! l index of the second spherical harmonic.
!
      integer(4) :: l3 ! l index of the third spherical harmonic.
!
      integer(4) :: m1 ! m index of the first spherical harmonic.
!
      integer(4) :: m2 ! m index of the second spherical harmonic.
!
      integer(4) :: m3 ! m index of the third spherical harmonic.
!
!
! !LOCAL VARIABLES:

      integer(4) :: ll1 ! index of the value of Y_l1m1 in the matrix 
!                         shpar: ll1=l1^2+l1+m1+1
!
      integer(4) :: ll2 ! index of the value of Y_l2m2 in the matrix 
!                         shpar: ll2=l2^2+l2+m2+1
!
      integer(4) :: ll3 ! index of the value of Y_l3m3 in the matrix 
!                         shpar: ll3=l3^2+l3+m3+1
!
      integer(4) :: ileb !  index of the Lebedev-Laikov grid point
!      
      complex(8) :: yprod12 ! the product Y_l1m1*Y_l3m3
      
      real(8) :: yp12r,yp12i
      
      real(8) :: y3r,y3i
!
      real(8) :: fyint   ! temporary storage of the four 
!                                      spherical harmonics integral
!
      real(8) :: threeylm ! final value of the four 
!                                      spherical harmonics integral
!
!
! !REVISION HISTORY:
! 
! Created Jan. 2004
! Last Modified: Feb. 2nd. 2004
!
!EOP
!BOC

      fyint = 0.0d+0

      ll1=l1*l1+l1+m1+1
      ll2=l2*l2+l2+m2+1
      ll3=l3*l3+l3+m3+1
      
      do ileb=1,n_angrid
        yprod12=sphar(ileb,ll1)*sphar(ileb,ll2)
        yp12r=real(yprod12)
        yp12i=aimag(yprod12)
        y3r=real(sphar(ileb,ll3))
        y3i=aimag(sphar(ileb,ll3))
        fyint=fyint+4.0d+0*pi*w_angrid(ileb)*(yp12r*y3r+yp12i*y3i)
      enddo
      threeylm=fyint
      return
      end function threeylm
!EOC

