!BOP
!
! !ROUTINE: coul_strcnst
!
! !INTERFACE:
      subroutine coul_strcnst(iq,lammax)

! !DESCRIPTION:
!
!This subroutine calculates the lattice sums $\Sigma^{a,a'}_{\lambda,\mu}(\vec{q})=\sum_{\vec{R}}{\frac{e^{i\vec{q}\cdot\left(\vec{R}+\vec{r}_{aa'}\right)}}%
!{|\vec{R}+\vec{r}_{aa'}|^{(\lambda+1)}}Y_{\lambda\mu}\left(\hat{R}_{aa'}\right)}$ using the method
!described in appendix \ref{ewaldmeth}, in particular equation \ref{strconstdef}.
!
! !USES:

      use barcoul,   only: eta, gcf,rcf, rstr, sgm, &
     &                     rbas, genrstr
      use constants, only: czero, imag, pi
      use kpoints,   only: idvq, qlist
      use struk,     only: alat, ndf, pos, vi, br2, ortho 
      use task,      only: fid_outdbg

! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: iq        ! index of the q-point for which sigma is calculated
      integer(4), intent(in) :: lammax ! Maximum value of lambda
     
! !LOCAL VARIABLES:

      logical :: ldbg=.true.
      integer(4) :: np ! Number of points for the real space summation
      integer(4) :: ng ! Number of points for the reciprocal  space summation
      integer(4) :: idf  !!  indexes all atoms
      integer(4) :: jdf  !!  indexes all atoms
      integer(4) :: ijdf !!  index of the matrix element (idf,jdf) in 
                         !!  packed form sigma(idf,jdf)=sgm(ijdf) with 
                         !!  ijdf= idf+jdf*(jdf-1)/2
      integer(4) :: i
      integer(4) :: i1
      integer(4) :: iqvec(3)
      integer(4) :: lam
      integer(4) :: mu
      integer(4) :: ilmu
      integer(4) :: fid

      real(8) :: rleng               ! the length of rpaa
      real(8) :: qtraa               ! the scalar product qvec.rpaa
      real(8) :: gausr               ! value of the gaussian function e^(-(rleng/eta)^2)
      real(8) :: pref                ! prefactor for the reciprocal lattice sum = 4 pi^(3/2)/v
      real(8) :: gleng               ! the length of gqv
      real(8) :: gqtra               ! the scalar product gqv.raa
      real(8) :: gausg               ! value of the gaussian function e^(-(gleng*eta/2)^2)
      real(8) :: gtolam              ! value of gleng^(lambda-2)/2^(lambda-0.5)
      real(8) :: gamlam              ! gamlam = Gamma[lambda+1/2]
      real(8) :: erfr
      real(8) :: gammaor             

      real(8), dimension(3) :: qvec  ! cartesian coords. of the q point
      real(8), dimension(3) :: raa   ! Vector going from atom 1 to atom 2.
      real(8), dimension(3) :: rdif
      real(8), dimension(3,3) :: rbs
      real(8), dimension(3) :: rpaa  ! the corresponding sum R+raa
      real(8), dimension(3) :: qtemp
      real(8), dimension(3) :: g     ! vector belonging to the reciprocal space lattice
      real(8), dimension(3) :: gqv   ! the corresponding sum G+qvec

      complex(8) :: ylam((lammax+1)*(lammax+1)) ! the values of the spherical harmonics
      complex(8) :: stmp1((lammax+1)*(lammax+1)) ! temporary allocation of the values of sigma
      complex(8) :: stmp2((lammax+1)*(lammax+1)) ! temporary allocation of  the values of sigma

      complex(8) :: expqdr                            ! e^(qvec.rpaa)
      complex(8) :: itolam               !i^lambda

            

! !EXTERNAL ROUTINES: 
      real(8), external :: derfc
      external k2cart

! !REVISION HISTORY:
!
! Created 21. Jan. 2004 by RGA
! Last Modified 3. Apr 2010 by JH
!
!EOP
!BOC      
    
      if(ldbg)  call boxmsg(fid_outdbg,'+',"structure constant (Sigma)")

      call k2cart(qlist(1:3,iq),idvq,qtemp)
      qvec(1:3)=-1.0d0*qtemp(1:3)

      do idf=1,ndf
        do jdf=idf,ndf

          rdif(:)=pos(:,idf)-pos(:,jdf)
          if(ortho)then
            raa(1)=rdif(1)*alat(1)
            raa(2)=rdif(2)*alat(2)
            raa(3)=rdif(3)*alat(3)
          else
            do i=1,3
              raa(i)=rbas(1,i)*rdif(1)+rbas(2,i)*rdif(2)+rbas(3,i)*rdif(3)
            enddo
          endif
          ijdf=idf+jdf*(jdf-1)/2

          !! Calculate all the R's such that R+r_aa < rcf
          rbs=transpose(rbas)
          call genrstr(rcf,raa,rbs,np)

          !! Initialize the temporal storage of the lattice sums
          stmp1=czero
          stmp2=czero

          !! Calculate the sum over the real space lattice
          do i1=1,np
            rpaa(1:3)=rstr(1:3,i1)
            rleng=rstr(4,i1)

            !! calculate the values of the spherical harmonics at rpaa
            call ylm(rpaa,lammax,ylam)
            qtraa=qvec(1)*rpaa(1)+qvec(2)*rpaa(2)+qvec(3)*rpaa(3)
            expqdr= cmplx(cos(qtraa),sin(qtraa),8)
            gausr = exp(-1.0d0*rleng*rleng/(eta*eta))
            erfr =  erfc(rleng/eta)
            gammaor = sqrt(pi)*erfr/rleng
            stmp1(1)=stmp1(1)+expqdr*cmplx(gammaor,0.0d0,8)*ylam(1)
            do lam=1,lammax
              gammaor=(dble(lam)-0.5d0)*gammaor/rleng+&
     &                rleng**(lam-2)*gausr/(eta**(2*lam-1))
              do mu=-lam,lam
                ilmu=lam*lam+lam+mu+1
                stmp1(ilmu)=stmp1(ilmu)+cmplx(gammaor,0.0d0,8)*     &
     &                       expqdr*ylam(ilmu)
              enddo ! mu
            enddo ! lam
          enddo ! i1

          !! calculate the vectors for the sum in reciprocal space
          qtemp(1:3)=-1.0d0*qvec(1:3)
          call genrstr(gcf,qtemp,br2,ng)

          !! Calculate the reciprocal lattice sum
          pref=4.0d0*pi*dsqrt(pi)*vi
          do i1=1,ng
            gqv(1:3)=rstr(1:3,i1)
            g(1:3)=gqv(1:3)-qtemp(1:3)
            gleng=rstr(4,i1)

            !!calculate the values of the spherical harmonics at rpaa
            call ylm(gqv,lammax,ylam)
            gqtra=g(1)*raa(1)+g(2)*raa(2)+g(3)*raa(3)
            expqdr=cmplx(dcos(gqtra),dsin(gqtra),8)
            gausg = dexp(-2.5d-1*eta*gleng*eta*gleng)
            gtolam = 1.0d0/(gleng*gleng)
            stmp2(1)=stmp2(1)+cmplx(pref*gtolam*gausg,0.0d0,8)*expqdr*ylam(1)
            itolam=cmplx(1.0d0,0.0d0,8)
            do lam=1,lammax
              gtolam=-1.0d0*gleng*gtolam/2.0d0
              itolam=itolam*imag
              do mu=-lam,lam
                ilmu=lam*lam+lam+mu+1
                stmp2(ilmu)=stmp2(ilmu)+cmplx(pref*gtolam*gausg,    &
     &                       0.0d0,8)*expqdr*ylam(ilmu)*itolam
              enddo ! mu
            enddo ! lam
          enddo !i1

          gamlam=dsqrt(pi)

          stmp1(1)=stmp1(1)*cmplx(1.0d0/gamlam,0.0d0,8)
          stmp2(1)=stmp2(1)*cmplx(1.0d0/gamlam,0.0d0,8)
          sgm(1,ijdf)=stmp1(1)+stmp2(1)

          do lam=1,lammax
            gamlam= 5.0d-1*dble(2*lam-1)*gamlam
            do mu=-lam,lam
              ilmu=lam*lam+lam+mu+1
              stmp1(ilmu)=stmp1(ilmu)*cmplx(1.0d0/gamlam,0.0d0,8)
              stmp2(ilmu)=stmp2(ilmu)*cmplx(1.0d0/gamlam,0.0d0,8)
              sgm(ilmu,ijdf)=stmp1(ilmu)+stmp2(ilmu)
            enddo
          enddo

          if(idf.eq.jdf) sgm(1,ijdf)=sgm(1,ijdf)-1.0d0/(eta*pi)

        enddo ! jdf
      enddo ! idf

      if(ldbg) then 
        do idf=1,ndf
          do jdf=idf,ndf
            ijdf=idf+jdf*(jdf-1)/2
            do lam=1,lammax
              do mu=-lam,lam
                ilmu=lam*lam+lam+mu+1
                write(fid_outdbg,'(4i5,2e12.4)') mu,lam,idf,jdf,sgm(ilmu,ijdf)
              enddo
            enddo
          enddo
        enddo
      endif 
   
   10 format('at1',1x,'at2',2x,'l',2x,'m',15x,'Sigma',27x,'real sum',  &
     &         21x,'reciprocal sum')
   11 format(4(i3,1x),e16.8,e16.8,1x,e16.8,e16.8,1x,e16.8,e16.8)
      return
      end subroutine coul_strcnst
