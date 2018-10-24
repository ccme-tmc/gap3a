!BOP
!
! !ROUTINE: calcgcoef
!
! !INTERFACE:
      subroutine calcgcoef(maxj)

! !DESCRIPTION:
!      
! This subroutine calculates the gaunt coefficients:
! 
! \begin{equation}
! \mathcal{G}^{LM}_{ll',mm'}=\int\limits_0^{2\pi}{\int\limits_0^{\pi}{%
! Y_{lm}(\theta,\phi) Y_{l'm'}(\theta,\phi)
! Y^*_{LM}(\theta,\phi)\sin(\theta)d\theta}d\phi}
! \end{equation}
! 
! for $l$ and $l' = 0,1,$... \verb"maxj". The integral is done numerically
! on a special grid (see Lebedev-Laicov.f90).
! The values are calculated only for $l\ge l'$ and $m' \ge 0$. 
! 
! The storage is optimized by saving the values in a vector
! (\verb"cgcoef") only for those coefficients that are different from zero.
! The size of the vector is:
!\begin{equation}
!n=\tfrac{1}{60}(l_{max}+1)(l_{max}+2)(l_{max}+3)(16l_{max}^2+29l_{max}+10)
!\end{equation}
!
!The gaunt coefficient $\mathcal{G}^{LM}_{ll',mm'}$ can be accesed directly at
!\verb"cgcoef(i)" by applying the function:
!
!\begin{equation}
!\begin{aligned}
!i=&\tfrac{1}{60}(16l^2-3l-3)(l+2)(l+1)l+\tfrac{1}{3}ll'(l'+1)(4l'-1)+%
!\tfrac{1}{6}l'(l'-1)(4l'+7)+\\
!&(2l+1)(l'+1)(L-l+l')+(l'+1)(m+l)+m'+l'+1 
!\end{aligned}
!\end{equation}
!
!of course, $M=m+m'$ is already taken into account. For the cases $l'>l$
!and $m'<0$ see the \verb"getcgcoef" subroutine.
!
 
! !USES:
             
      use mixbasis,     only: cgcoef, init_gaunt
      use angintegrals, only: end_angint

! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: maxj
      
! !LOCAL VARIABLES:
      
      integer(4) :: i
      integer(4) :: l1 
      integer(4) :: l2 
      integer(4) :: l3
      integer(4) :: m1 
      integer(4) :: m2 
      integer(4) :: m3 
      integer(4) :: ntot
      integer(4) :: ngrid
      
 
!
! !EXTERNAL ROUTINES: 
!


      real(8), external :: threeylm

      external prep_ant_int
      

 
! 
! !INTRINSIC ROUTINES: 
!


      intrinsic iabs

!
! !REVISION HISTORY:
!
! Created: Apr. 2004 by RGA
! Last modified: May 21st. 2004 by RGA
!
!EOP
!BOC

      ntot=(maxj+1)*(maxj+2)*(maxj+3)*(16*maxj*maxj+29*maxj+10)/60
      ngrid=(4*maxj+1)*(4*maxj+1)/3

      call prep_ang_int(2*maxj,ngrid)
      
      call init_gaunt(ntot)
      i=0      
      do l1=0,maxj
        do l2=0,l1
          do l3=l1-l2,l1+l2
            do m1=-l1,l1
              do m2=0,l2
                m3=m1+m2
                i=i+1
                if(mod(l1+l2+l3,2).eq.0)then
                  if(iabs(m3).le.l3)then
                    cgcoef(i)=threeylm(l1,l2,l3,m1,m2,m3)
                  else
                    cgcoef(i)=0.0d0
                  endif
                else
                  cgcoef(i)=0.0d0
                endif    
              enddo
            enddo 
          enddo
        enddo
      enddo

      call end_angint

      end subroutine calcgcoef
!EOC
