!BOP
!
! !ROUTINE: calcmmatcv
!
! !INTERFACE:
      subroutine calcmmatcv(nv,ik,isp,mcv)
      
! !DESCRIPTION:
!
! This subroutine calculates the momentum matrix element between the core
! state nc and the band nv
!      
! !USES:
      
      use core,        only: corind, ncg
      use constants,   only: czero, imag
      use eigenvec,    only: alfa, beta, gama
      use lapwlo,      only: nt,nlo_at
      use mommat,      only: iucl1ul, iucl1udl, iucl1ulol, iuclul1,   &
     &                       iucludl1, iuclulol1, iul1ucl, iudl1ucl,  &
     &                       iulol1ucl, iulucl1, iudlucl1, iulolucl1, &
     &                       lhermit
      use struk,       only: mult,nat

! !INPUT PARAMETERS:

      
      integer, intent(in) :: nv  ! Second band index
      integer, intent(in) :: ik  ! Index of the irreducible k-point
      integer, intent(in) :: isp ! Index of spin
      complex(8), intent(out):: mcv(3,ncg)

! !OUTPUT PARAMETERS:
      
      complex(8) :: mxcv,mxvc ! MT contribution to the x component      
      complex(8) :: mycv,myvc ! MT contribution to the x component      
      complex(8) :: mzcv,mzvc ! MT contribution to the x component 

! !LOCAL VARIABLES:           

      integer :: nc,iclm  ! core state index
      integer :: iat      ! run over inequivalent atoms
      integer :: icg      ! over all core states 
      integer :: idf      ! run over equivalent atoms
      integer :: idf0     ! run over all atoms
      integer :: l        ! angular momentum quantum number
      integer :: m        ! quantum number of the z component of the angular momentum
      integer :: lm, l1m1, l1mm1, l1m,lm1m,lm1mm1,lm1m1
      integer :: ilo      ! index for LO 
      
      real(8) :: denom1,denom2,renorm
      real(8) :: flm1,flm2,flm3,flm4,flm5,flm6 ! F_lm coefficients of the mmat expansion in 

      complex(8) :: xpy1,xpy2,xmy1,xmy2,z1,z2
      complex(8) :: pxpy,pxmy,pz
      complex(8) :: mxpy,mxmy,mz
      complex(8) :: p(3)
 
! !INTRINSIC ROUTINES: 
      
      intrinsic conjg
      intrinsic sqrt


! !REVISION HISTORY:
! 
! Created 2. Sep. 2004 by RGA
!
!EOP  
!BOC  
      mxpy = czero
      mxmy = czero
      mz = czero
      xpy1= czero
      xpy2=czero
      xmy1 = czero
      xmy2 = czero
      z1 = czero
      z2 = czero
      
!     -------------------------------------------------------
!                   contribution from MT spheres
!                               BEGIN
!     -------------------------------------------------------
!
!     Loop over core states
!
      idf0=0
      do icg = 1, ncg
        iat=corind(1,icg)
        idf=corind(2,icg)

        !!reset iclm if it is a new atom 
        if(idf.ne.idf0)then
          iclm=0
          idf0=idf
        endif  
        nc=corind(3,icg)
        l=corind(4,icg)
        m=corind(5,icg)
        iclm=iclm+1
        denom1=dble((2*l+1)*(2*l+3))
        denom2=dble((2*l-1)*(2*l+1))

        lm = l*l+l+m+1      !* -> l,m
        l1mm1=lm+2*l+1      !* -> l+1,m-1
        l1m = l1mm1 + 1     !* -> l+1,m
        l1m1 = l1m + 1     !* -> l+1,m+1
        lm1mm1=lm-2*l-1    !* -> l-1,m-1
        lm1m = lm1mm1 + 1  !*  -> l-1,m
        lm1m1 = lm1m + 1   !* -> l-1,m+1


        !!calculate the F_lm coefficients

        flm1=-sqrt(dble((l+m+1)*(l+m+2))/denom1)
        flm2=sqrt(dble((l-m-1)*(l-m))/denom2)
        flm3=sqrt(dble((l-m+1)*(l-m+2))/denom1)
        flm4=-sqrt(dble((l+m-1)*(l+m))/denom2)
        flm5=sqrt(dble((l-m+1)*(l+m+1))/denom1)
        flm6=sqrt(dble((l-m)*(l+m))/denom2)
!
!       calculate the matrix elements <core/p/valence>
!
        if (l.gt.0.and.l+m-2.ge.0) then 
          xpy1=alfa(nv,lm1mm1,idf)*iucl1ul(nc,iat,isp)     &
     &        +beta(nv,lm1mm1,idf)*iucl1udl(nc,iat,isp)
          do ilo=1,nLO_at(l-1,iat) 
            xpy1=xpy1+gama(nv,ilo,lm1mm1,idf)*iucl1ulol(ilo,nc,iat,isp)
          enddo  
        endif
            
        if (l.gt.0.and.l-m-2.ge.0)then
          xmy1=alfa(nv,lm1m1,idf)*iucl1ul(nc,iat,isp)     &
     &        +beta(nv,lm1m1,idf)*iucl1udl(nc,iat,isp)
          do ilo=1,nLO_at(l-1,iat)
            xmy1= xmy1+gama(nv,ilo,lm1m1,idf)*iucl1ulol(ilo,nc,iat,isp) 
          enddo 
        endif
            
        if (l.gt.0.and.(l+m-1.ge.0).and.(l-m-1.ge.0))then
          z1=alfa(nv,lm1m,idf)*iucl1ul(nc,iat,isp)       &
     &      +beta(nv,lm1m,idf)*iucl1udl(nc,iat,isp)
          do ilo=1,nLO_at(l-1,iat) 
            z1=z1+gama(nv,ilo,lm1m,idf)*iucl1ulol(ilo,nc,iat,isp) 
          enddo 
        endif

        xpy2=alfa(nv,l1mm1,idf)*iuclul1(nc,iat,isp)         &
     &      +beta(nv,l1mm1,idf)*iucludl1(nc,iat,isp)

        xmy2=alfa(nv,l1m1,idf)*iuclul1(nc,iat,isp)          &
     &      +beta(nv,l1m1,idf)*iucludl1(nc,iat,isp)

        z2=alfa(nv,l1m,idf)*iuclul1(nc,iat,isp)             &
     &    +beta(nv,l1m,idf)*iucludl1(nc,iat,isp)

        do ilo=1,nLO_at(l+1,iat) 
          xpy2 = xpy2 + gama(nv,ilo,l1mm1,idf)*iuclulol1(ilo,nc,iat,isp)
          xmy2 = xmy2 + gama(nv,ilo,l1m1, idf)*iuclulol1(ilo,nc,iat,isp)
          z2   = z2   + gama(nv,ilo,l1m,  idf)*iuclulol1(ilo,nc,iat,isp)
        enddo 

!
!       we make use of the fact that:
!       F^1-(l-1,m-1)=f^4_{lm}
!       F^2-(l-1,m-1)=f^3_{lm}
!       F^3-(l-1,m-1)=f^2_{lm}
!       F^4-(l-1,m-1)=f^1_{lm}
!       F^5-(l-1,m-1)=f^6_{lm}
!       F^6-(l-1,m-1)=f^5_{lm}
!
        pxpy = flm4 * xpy1 + flm3 * xpy2
        pxmy = flm2 * xmy1 + flm1 * xmy2
        pz   = flm6 * z1   + flm5 * z2 
           
        mxcv=-5.0d-1*imag*(pxpy+pxmy)
        mycv=5.0d-1*(pxmy-pxpy)
        mzcv=-imag*pz

!       write(64,10)nc,nv,m,mxcv,mycv,mzcv
            
!
!       calculate the matrix elements <valence/p/core>
!         
        xpy1=conjg(alfa(nv,l1m1,idf))*iul1ucl(nc,iat,isp)+   &
     &       conjg(beta(nv,l1m1,idf))*iudl1ucl(nc,iat,isp)

        xmy1=conjg(alfa(nv,l1mm1,idf))*iul1ucl(nc,iat,isp)+  &
     &       conjg(beta(nv,l1mm1,idf))*iudl1ucl(nc,iat,isp)

        z1=conjg(alfa(nv,l1m,idf))*iul1ucl(nc,iat,isp)+      &
     &     conjg(beta(nv,l1m,idf))*iudl1ucl(nc,iat,isp)

        do ilo=1,nLO_at(l+1,iat)  
          xpy1=xpy1+conjg(gama(nv,ilo,l1m1, idf))*iulol1ucl(ilo,nc,iat,isp) 
          xmy1=xmy1+conjg(gama(nv,ilo,l1mm1,idf))*iulol1ucl(ilo,nc,iat,isp)
          z1  =z1 + conjg(gama(nv,ilo,l1m,  idf))*iulol1ucl(ilo,nc,iat,isp) 
        enddo 

        if (l.gt.0.and.l-m-2.ge.0) then 
          xpy2=conjg(alfa(nv,lm1m1,idf))*iulucl1(nc,iat,isp)            &
     &        +conjg(beta(nv,lm1m1,idf))*iudlucl1(nc,iat,isp)
          do ilo=1,nLO_at(l-1,iat) 
            xpy2=xpy2+conjg(gama(nv,ilo,lm1m1,idf))*iulolucl1(ilo,nc,iat,isp)
          enddo 
        endif 

        if (l.gt.0.and.l+m-2.ge.0) then 
          xmy2=conjg(alfa(nv,lm1mm1,idf))*iulucl1(nc,iat,isp) &
     &        +conjg(beta(nv,lm1mm1,idf))*iudlucl1(nc,iat,isp)

          do ilo=1,nLO_at(l-1,iat) 
            xmy2=xmy2+conjg(gama(nv,ilo,lm1mm1,idf))*iudlucl1(nc,iat,isp)
          enddo 
        endif 
            
        if (l.gt.0.and.(l+m-1.ge.0).and.(l-m-1.ge.0))then
          z2=conjg(alfa(nv,lm1m,idf))*iulucl1(nc,iat,isp)               &
     &      +conjg(beta(nv,lm1m,idf))*iudlucl1(nc,iat,isp)
          do ilo=1,nLO_at(l-1,iat) 
            z2=z2+conjg(gama(nv,ilo,lm1m,idf))*iulolucl1(ilo,nc,iat,isp)
          enddo 
        endif

        pxpy = flm1 * xpy1 + flm2 * xpy2
        pxmy = flm3 * xmy1 + flm4 * xmy2
        pz   = flm5 * z1   + flm6 * z2 
      
        mxvc=-5.0d-1*imag*(pxpy+pxmy)
        myvc=5.0d-1*(pxmy-pxpy)
        mzvc=-imag*pz

        if(lhermit) then             
          p(1)=0.5d0*(mxcv+conjg(mxvc))           
          p(2)=0.5d0*(mycv+conjg(myvc))           
          p(3)=0.5d0*(mzcv+conjg(mzvc)) 
        else 
          p(1) = mxcv 
          p(2) = mycv
          p(3) = mzcv
        endif 

        mcv(1:3,icg)=p(1:3)
      enddo ! icg
      
      return

      end subroutine calcmmatcv
!EOC      
