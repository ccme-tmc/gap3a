!BOP
!
! !ROUTINE: getcgcoef
!
! !INTERFACE:      
      real(8) function getcgcoef(l1,l2,l3,m1,m2)
      
! !DESCRIPTION:
!
!This subroutine gets the gaunt coefficient $G^{l3,m1+m2}_{l1,l2,m1,m2}$
!from the vector \verb"cgcoef" by:
!\begin{equation}
!G^{l3,m1+m2}_{l1,l2,m1,m2}=\alpha \verb"cgcoef"(i)
!\end{equation}
!calculating i by:
!\begin{equation}
!\begin{aligned}
!i=&\tfrac{1}{60}(16l^2-3l-3)(l+2)(l+1)l+\tfrac{1}{3}ll'(l'+1)(4l'-1)+%
!\tfrac{1}{6}l'(l'-1)(4l'+7)+\\
!&+(2l+1)(l'+1)(L-l+l')+(l'+1)(m+l)+m'+l'+1 
!\end{aligned}
!\end{equation}
!where
!
!\begin{subequations}
!\begin{align}
!l&=l1 & l'&=l2 & m&=m1 & m'&=m2 & \alpha&=1 &\text{ if $l1 \ge l2$} \\
!l&=l2 & l'&=l1 & m&=m2 & m'&=m1 & \alpha&=1 &\text{ if $l1 < l2$} \\
! &    &   &    & m&=m1 & m'&=m2 & \alpha&=1 &\text{ if $m2 \ge 0$} \\
! &    &   &    & m&=-m1 & m'&=-m2 & \alpha&=(-1)^{l+l'-L} &\text{ if $m2 < 0$}
!\end{align}
!\end{subequations}
!
! !USES:      

      use mixbasis, only: cgcoef

! !INPUT PARAMETERS:
      
      implicit none

      integer(4), intent(in) :: l1
      integer(4), intent(in) :: l2
      integer(4), intent(in) :: l3
      integer(4), intent(in) :: m1
      integer(4), intent(in) :: m2
      
! !LOCAL VARIABLES:

      integer(4) :: j1,j2,mj1,mj2
      integer(4) :: par,ing
      integer(4) :: ind1,ind2,ind3,ind4

      real(8) :: fact
      
      logical :: trcond

 
!
! !INTRINSIC ROUTINES: 
!


      intrinsic mod
      intrinsic abs

!
! !REVISION HISTORY:
!
! Created: Apr. 2004 by RGA
! Last modified  May 21st. 2004 by RGA
!
!EOP
!BOC    
     
      par=mod(abs(l1+l2-l3),2)    
      fact=1.0d0
      trcond=(abs(m1+m2).le.l3).and.(abs(l1-l2).le.l3).and.(l1+l2.ge.l3)
      if(trcond)then
        if(l1.lt.l2)then
          j1=l2
          mj1=m2
          j2=l1
          mj2=m1
!          fact=-2.0d0*par+1.0d0
        else
          j1=l1
          mj1=m1
          j2=l2
          mj2=m2
        endif
        if(mj2.lt.0)then
          mj2=-mj2
          mj1=-mj1
          fact=(-2.0d0*par+1.0d0)
        endif
        ind1=(16*j1*j1-3*j1-3)*(j1+2)*(j1+1)*j1/60
        ind2=j1*j2*(j2+1)*(4*j2-1)/3
        ind3=j2*(j2-1)*(4*j2+7)/6
        ind4=(2*j1+1)*(j2+1)*(l3-j1+j2)
        ing=ind1+ind2+ind3+ind4+(j2+1)*(mj1+j1)+mj2+j2+1 
        getcgcoef=fact*cgcoef(ing)   
      else
        getcgcoef=0.0d0
      endif  
!      write(*,10)l1,l2,m1,m2,j1,j2,l3,mj1,mj2,ing,getcgcoef
     
!   10 format('gaunt',9i4,i8,d15.7) 
      end function getcgcoef
!EOC      
