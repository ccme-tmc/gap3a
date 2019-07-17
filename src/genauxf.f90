!BOP
!
! !ROUTINE: genauxf
!
! !INTERFACE: 
      subroutine genauxf(iq,f1,f2,alfa)

! !DESCRIPTION:
!
! Given the index \texttt{iq} of $\vec{q}$, this
! subroutine generates the auxiliary functions $F_1(\vec{q})$ and $F_2(\vec{q})$
! according to the formulas:
!
! \begin{subequations}\label{genauxf-01}
! \begin{align}
! F_1(\vec{q})=&\frac{1}{\Omega}\sum\limits_i{\frac{e^{-\alpha|\vec{q}+\vec{G}_i|^2}}{|\vec{q}+\vec{G}_i|}}\\
! F_2(\vec{q})=&\frac{1}{\Omega}\sum\limits_i{\frac{e^{-\alpha|\vec{q}+\vec{G}_i|^2}}{|\vec{q}+\vec{G}_i|^2}}
! \end{align}
! \end{subequations}
! 
! where $\alpha=\left(\frac{\Omega}{6\pi^2}\right)^{\frac{1}{3}}$
!
! !USES:

      use constants, only: pi
      use kpoints,   only: qlist, idvq, nkp
      use recipvec,  only: ngqbarc, indgqlen, gqleng, indgq
      use struk,     only: vi
      
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: iq    ! Index of the q-vector
      
! !OUTPUT PARAMETERS:

      integer(4) :: igl

      real(8),    intent(out) :: f1   ! Vaule of the auxiliary function F_1
!                                       at iq       
      real(8),    intent(out) :: f2   ! Vaule of the auxiliary function F_2
!                                       at iq       
      real(8),    intent(out) :: alfa ! Parameter of the function
      
! !LOCAL VARIABLES:

      integer(4) :: ipw              ! (Counter) runs over plane waves.
      integer(4) :: ipwin            ! Initial value of ipw

      real(8) :: modgpq              ! Length of q+G squared
      real(8) :: expagpq             ! exp(alfa*|q+g|^2

! !INTRINSIC ROUTINES: 

      intrinsic exp
      intrinsic sqrt

! !REVISION HISTORY:
!
! Created 23.06.05 by RGA.
! Last modified Nov 2006 by RGA
!
!EOP
!BOC
!
!     Initializations
!
      f1=0.0d0
      f2=0.0d0
      ipwin=1
      alfa=(1.0d0/(6.0d0*pi*pi*vi))**(1.0d0/3.0d0)
      if(iq.eq.1) ipwin=2 
!
!     Loop over G-vectors
!
      do ipw = ipwin, ngqbarc(iq)
!
!       Calculate the cartessian coordinates of G
!
        igl=indgqlen(ipw,iq)
!
!       Calculate the length of G+q squared 
!
        modgpq=gqleng(igl,iq)*gqleng(igl,iq)
        expagpq=exp(-alfa*modgpq)
!
!       Accumulate the terms into f1 and f2
!
        f1=f1+expagpq/sqrt(modgpq)
        f2=f2+expagpq/modgpq

      enddo ! ipw
      f1 = f1 / dble(nkp)
      f2 = f2 / dble(nkp)
      
      return
      
      end subroutine genauxf
!EOC                
