!BOP
!
! !ROUTINE: genauxf
!
! !INTERFACE: 
      subroutine genauxf_2d(iq,iaxis,alfa,f1,f2)

! !DESCRIPTION:
!
! Given the index \texttt{iq} of $\vec{q}$, this
! subroutine generates the auxiliary functions $F_1(\vec{q})$ and $F_2(\vec{q})$
! according to the formulas:
!
! \begin{subequations}\label{genauxf-01}
! \begin{align}
! F_1(\vec{q})=&\frac{1}{\Omega}\sum\limits_i{\frac{e^{-\alpha^2|\vec{q}+\vec{G}_i|^2}}{\sqrt{|\vec{q}+\vec{G}_i|}}}\\
! F_2(\vec{q})=&\frac{1}{\Omega}\sum\limits_i{\frac{e^{-\alpha|\vec{q}+\vec{G}_i|^2}}{|\vec{q}+\vec{G}_i|}}
! \end{align}
! \end{subequations}
! 
! where $\alpha=\left(\frac{\Omega}{6\pi^2}\right)^{\frac{1}{3}}$
!
! !USES:

      use kpoints,   only: qlist, idvq, nkp
      use recipvec,  only: ngqbarc, indgqlen, gqleng, indgq,gindex
      use struk,     only: vi
      
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: iq    ! Index of the q-vector
      integer(4), intent(in) :: iaxis ! index of axis along which v is cutoff
      real(8),    intent(in) :: alfa  ! Parameter of the function
      
! !OUTPUT PARAMETERS:

      integer(4) :: igl

      real(8),    intent(out) :: f1   ! Vaule of the auxiliary function F_1
!                                       at iq       
      real(8),    intent(out) :: f2   ! Vaule of the auxiliary function F_2
!                                       at iq       
      
! !LOCAL VARIABLES:

      integer(4) :: ipw              ! (Counter) runs over plane waves.
      integer(4) :: ipwin            ! Initial value of ipw

      real(8) :: modgpq              ! Length of q+G squared
      real(8) :: expagpq             ! exp(alfa*|q+g|^2
      integer(4),dimension(3) :: igvec

! !EXTERNAL ROUTINES: 
      !real(8),external ::

! !INTRINSIC ROUTINES: 

      intrinsic exp
      intrinsic sqrt

! !REVISION HISTORY:
!
! Created 23.06.05 by RGA.
!         Nov 2006 by RGA
!      24 Oct 2018 by ZMY:
!               change alfa as input instead of output
!      09 Sep 2019 by ZMY:
!
!EOP
!BOC
!
!     Initializations
!
      if((iaxis>3).or.(iaxis<1)) then
        call errmsg(.true.,"genauxf_2d","invalid cutoff axis number")
      endif
      f1=0.0d0
      f2=0.0d0
      ipwin=1
      ! skip G=0 for q=0
      if(iq.eq.1) ipwin=2 
!
!     Loop over G-vectors
!
      do ipw = ipwin, ngqbarc(iq)
!
!       Calculate the cartessian coordinates of G
!
        igvec(1:3) = gindex(1:3,indgq(ipw,iq))
        ! skip ipw that has finite Gz component
        if(igvec(iaxis).ne.0) continue
        igl=indgqlen(ipw,iq)
!
!       Calculate the length of G+q squared 
!
        modgpq=gqleng(igl,iq)
        !expagpq=exp(-alfa*modgpq)
        expagpq=exp(-alfa*alfa*modgpq*modgpq)
!
!       Accumulate the terms into f1 and f2
        f1=f1+expagpq/sqrt(modgpq)
        f2=f2+expagpq/modgpq

      enddo ! ipw

      ! nkp = N_c
      f1 = f1 / dble(nkp)
      f2 = f2 / dble(nkp)
      
      return
      
      end subroutine genauxf_2d
!EOC                
