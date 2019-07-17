!BOP
!
! !ROUTINE: intipw
!
! !INTERFACE:
      subroutine intipw
      
! !DESCRIPTION:

!This function calculates the integral of a plane wave with wave vector
!belonging to the reciprocal Bravais lattice in the
!interstitial region by the difference between the integral over the whole
!unit cell and the Muffin Tin spheres:
!
!\begin{equation}
!\frac{1}{\Omega}\int\limits_{I}{e^{i\vec{G}\cdot\vec{r}}d^3r}=%
!\frac{1}{\Omega}\int\limits_{\Omega}{e^{i\vec{G}\cdot\vec{r}}d^3r}-%
!\frac{1}{\Omega}\sum\limits_a{\int\limits_{MT_a}{e^{i\vec{G}\cdot\vec{r}}d^3r}}
!\end{equation}
!being:
!
!\begin{equation}
!\frac{1}{\Omega}\int\limits_{\Omega}{e^{i\vec{G}\cdot\vec{r}}d^3r}=\delta_{\vec{G},0}
!\end{equation}
!
!and
!
!\begin{equation}
!\frac{1}{\Omega}\int\limits_{MT_a}{e^{i\vec{G}\cdot\vec{r}}d^3r}=\frac{V_{MT_a}}{\Omega}\left[\delta_{\vec{G},0}+%
!3(1-\delta_{\vec{G},0})\frac{\sin(|\vec{G}|R_a)-(|\vec{G}|R_a)%
!\cos(|\vec{G}|R_a)}{(|\vec{G}|R_a)^3} e^{i\vec{G}\cdot\vec{r}_a}\right]
!\end{equation}
!
! !USES:
      
      use constants, only: czero,pi
      use struk,     only: nat,mult, pos, rmt
      use mixbasis,  only: ipwint, igm
      use recipvec,  only: npw, gindex, ig0
      
! !INPUT PARAMETERS:

      implicit none
      

! !LOCAL VARIABLES:

      integer(4) :: idxg
      
      integer(4) :: ig(3)
      

! !EXTERNAL ROUTINES:      
      
      external k2cart

! !REVISION HISTORY:
!
! Created 30. April 2004 by RGA
!
!EOP
!BOC
      allocate(ipwint(1:npw))
      do idxg=1,npw
        ig(1:3)=gindex(:,idxg)
        call int1ipw(ig,ipwint(idxg))
      enddo      
!      do idxg=1,npw
!        write(*,*)idxg,ipwint(idxg)
!      enddo  
!      ig(1:3)=0
!      sum = 0.0d0
!      do ig1=-igm(1),igm(1)
!        ig(1)=ig1
!        do ig2=-igm(2),igm(2)
!          ig(2)=ig2
!          do ig3=-igm(3),igm(3)
!            ig(3)=ig3
!            sum=sum+ipwint(ig(1),ig(2),ig(3))*conjg(ipwint(ig(1),ig(2),ig(3)))
!          enddo
!        enddo
!      enddo      
!      write(*,*)'sumipwint =', sum
!      write(*,10)ig,ipwint(0,0,0)
!   10 format(i4,'(',g18.10,',',g18.10,') ')
      
      end subroutine intipw
      
!EOC      
              
            

