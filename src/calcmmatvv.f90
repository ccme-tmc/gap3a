!BOP
!
! !ROUTINE: calcmmatvv
!
! !INTERFACE:
      subroutine calcmmatvv(n1,n2,ik,isp,p12)
      
! !DESCRIPTION:
!
!This subroutine calculates the contribution of the MT-Spheres to the
!momentum matrix element of the bands n1 and n2
!      
! !USES:
     
      use bands,       only: nv
      use constants,   only: czero, imag
      use eigenvec,    only: alfa, beta, gama,zzk
      use kpoints,     only: klist, idvk, kpirind
      use lapwlo,      only: nt,nLO_at 
      use mixbasis,    only: ipwint
      use mommat,     only: iul1ul,iul1udl,iudl1ul,iudl1udl,iulul1,  &
     &                       iuludl1,iudlul1,iudludl1,iul1ulol,       &
     &                       iulol1ul,iulol1udl,iudl1ulol,iulol1ulol, &
     &                       iululol1,iulolul1,iuloludl1,iulolulol1,  &
     &                       iudlulol1 
      use recipvec,    only: ig0, gindex, indgk
      use struk,       only: nat,mult
      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: n1  ! First band index
      integer, intent(in) :: n2  ! Second band index
      integer, intent(in) :: ik  ! index for the k-point
      integer, intent(in)  :: isp ! index for the spin 
! !OUTPUT PARAMETERS:
      complex(8), intent(out) :: p12(3) ! MT contribution to the x/y/z component      

! !LOCAL VARIABLES:           
      
      integer :: iat    ! run over inequivalent atoms
      integer :: ieq    ! run over equivalent atoms
      integer :: idf    ! run over all atoms
      integer :: l      ! angular momentum quantum number
      integer :: m      ! quantum number of the z component of the  angular momentum
      integer :: lm, l1m1, l1mm1, l1m  !! collective index for (l,m) etc. 
      integer :: ilo,jlo ! index for LO 

      integer :: ig1,ig2, irk,igv1(1:3),igv2(1:3),ikv(1:3)
      integer :: iv1, iv2, iv3

      real(8), dimension(3) :: kv
      real(8), dimension(3) :: gv1,gv2
      real(8), dimension(3) :: kpgv1,kpgv2
      real(8) :: denom
      real(8) :: flm1,flm3,flm5 ! F_lm coefficients of the mommat expansion in spherical harmonics

      complex(8) :: xpy1,xpy2,xmy1,xmy2,z1,z2
      complex(8) :: pxpy,pxmy,pz
      complex(8) :: mxpy,mxmy,mz 
      complex(8) :: a1(6),a2(6),b1(6),b2(6),c1(6),c2(6)
      complex(8) :: mt12(3),mpw(3) ! MT contribution to the x/y/z component      
      complex(8) :: intmom, modmom
 
! !INTRINSIC ROUTINES: 

      
      intrinsic conjg
      intrinsic sqrt


! !REVISION HISTORY:
! 
! Created 2. Sep. 2004 by RGA
!
!EOP  
!BOC  
!     -------------------------------------------------------
!                   contribution from MT spheres
!                               BEGIN
!     -------------------------------------------------------a
! Notations:

!    lm -> l,m                
!    l1mm1 -> l+1,m-1                
!    l1m -> l+1,m                
!    l1m1 -> l+1,m+1                
      mxpy = czero
      mxmy = czero
      mz = czero
      idf=0
      do iat=1,nat
        do ieq=1,mult(iat)
          idf=idf+1

          pxpy = czero
          pxmy = czero
          pz   = czero
          lm   = 0

          do l=0,nt-2 
            denom=dble((2*l+1)*(2*l+3))

            do m = -l, l

              lm = lm + 1              !! collective index for (l,m)
              l1mm1 = lm+2*l+1         !! (l+1,m-1)
              l1m = l1mm1 + 1          !! (l+1,m )
              l1m1 = l1m + 1           !! (l+1,m+1) 
!
!             calculate the F_lm coefficients in Eqs.(N2.6a,c,e)
!
              flm1= -sqrt(dble((l+m+1)*(l+m+2))/denom)
              flm3=  sqrt(dble((l-m+1)*(l-m+2))/denom)
              flm5=  sqrt(dble((l-m+1)*(l+m+1))/denom)

             xpy1= conjg(alfa(n1,l1m1,idf))*alfa(n2,lm,idf)*        &
     &              cmplx(iul1ul(l,iat,isp),0.0d0,8)+                    &
     &             conjg(alfa(n1,l1m1,idf))*beta(n2,lm,idf)*        &
     &              cmplx(iul1udl(l,iat,isp),0.0d0,8)+                   &
     &             conjg(beta(n1,l1m1,idf))*alfa(n2,lm,idf)*        &
     &              cmplx(iudl1ul(l,iat,isp),0.0d0,8)+                   &
     &             conjg(beta(n1,l1m1,idf))*beta(n2,lm,idf)*        &
     &              cmplx(iudl1udl(l,iat,isp),0.0d0,8)

              xpy2=conjg(alfa(n1,lm,idf))*alfa(n2,l1mm1,idf)*       &
     &             cmplx(iulul1(l,iat,isp),0.0d0,8)+                    &
     &             conjg(alfa(n1,lm,idf))*beta(n2,l1mm1,idf)*       &
     &             cmplx(iuludl1(l,iat,isp),0.0d0,8)+                   &
     &             conjg(beta(n1,lm,idf))*alfa(n2,l1mm1,idf)*       &
     &             cmplx(iudlul1(l,iat,isp),0.0d0,8)+                   &
     &             conjg(beta(n1,lm,idf))*beta(n2,l1mm1,idf)*       &
     &             cmplx(iudludl1(l,iat,isp),0.0d0,8)

              xmy1=conjg(alfa(n1,l1mm1,idf))*alfa(n2,lm,idf)*       &
     &             cmplx(iul1ul(l,iat,isp),0.0d0,8)+                    &
     &             conjg(alfa(n1,l1mm1,idf))*beta(n2,lm,idf)*       &
     &             cmplx(iul1udl(l,iat,isp),0.0d0,8)+                   &
     &             conjg(beta(n1,l1mm1,idf))*alfa(n2,lm,idf)*       &
     &             cmplx(iudl1ul(l,iat,isp),0.0d0,8)+                   &
     &             conjg(beta(n1,l1mm1,idf))*beta(n2,lm,idf)*       &
     &             cmplx(iudl1udl(l,iat,isp),0.0d0,8)

              xmy2=conjg(alfa(n1,lm,idf))*alfa(n2,l1m1,idf)*        &
     &             cmplx(iulul1(l,iat,isp),0.0d0,8)+                    &
     &             conjg(alfa(n1,lm,idf))*beta(n2,l1m1,idf)*        &
     &             cmplx(iuludl1(l,iat,isp),0.0d0,8)+                   &
     &             conjg(beta(n1,lm,idf))*alfa(n2,l1m1,idf)*        &
     &             cmplx(iudlul1(l,iat,isp),0.0d0,8)+                   &
     &             conjg(beta(n1,lm,idf))*beta(n2,l1m1,idf)*        &
     &             cmplx(iudludl1(l,iat,isp),0.0d0,8)

              z1 = conjg(alfa(n1,l1m,idf))*alfa(n2,lm,idf)*         &
     &             cmplx(iul1ul(l,iat,isp),0.0d0,8)+                    &
     &             conjg(alfa(n1,l1m,idf))*beta(n2,lm,idf)*         &
     &             cmplx(iul1udl(l,iat,isp),0.0d0,8)+                   &
     &             conjg(beta(n1,l1m,idf))*alfa(n2,lm,idf)*         &
     &             cmplx(iudl1ul(l,iat,isp),0.0d0,8)+                   &
     &             conjg(beta(n1,l1m,idf))*beta(n2,lm,idf)*         &
     &             cmplx(iudl1udl(l,iat,isp),0.0d0,8)

              z2 = conjg(alfa(n1,lm,idf))*alfa(n2,l1m,idf)*        &
     &             cmplx(iulul1(l,iat,isp),0.0d0,8)+                    &
     &             conjg(alfa(n1,lm,idf))*beta(n2,l1m,idf)*        &
     &             cmplx(iuludl1(l,iat,isp),0.0d0,8)+                   &
     &             conjg(beta(n1,lm,idf))*alfa(n2,l1m,idf)*        &
     &             cmplx(iudlul1(l,iat,isp),0.0d0,8)+                   &
     &             conjg(beta(n1,lm,idf))*beta(n2,l1m,idf)*        &
     &             cmplx(iudludl1(l,iat,isp),0.0d0,8)

!
!             the loop for the Local orbitals on the right side 
!
              do jlo = 1,nLO_at(l,iat) 
                xpy1 = xpy1 + &
     &               conjg(alfa(n1,l1m1,idf))*gama(n2,jlo,lm,idf)* &
     &                cmplx(iul1ulol(jlo,l,iat,isp),0.0d0,8)+         &
     &               conjg(beta(n1,l1m1,idf))*gama(n2,jlo,lm,idf)*      &
     &                cmplx(iudl1ulol(jlo,l,iat,isp),0.0d0,8)

                xmy1 = xmy1 + &
     &               conjg(alfa(n1,l1mm1,idf))*gama(n2,jlo,lm,idf)*&
     &                cmplx(iul1ulol(jlo,l,iat,isp),0.0d0,8)+                &
     &               conjg(beta(n1,l1mm1,idf))*gama(n2,jlo,lm,idf)*     &
     &                cmplx(iudl1ulol(jlo,l,iat,isp),0.0d0,8)

                z1 = z1 + &
     &               conjg(alfa(n1,l1m,idf))*gama(n2,jlo,lm,idf)*  &
     &                cmplx(iul1ulol(jlo,l,iat,isp),0.0d0,8)+                &
     &               conjg(beta(n1,l1m,idf))*gama(n2,jlo,lm,idf)*       &
                      cmplx(iudl1ulol(jlo,l,iat,isp),0.0d0,8)
              enddo
              
              do jlo = 1,nLO_at(l+1,iat) 
                xpy2 = xpy2 + & 
     &               conjg(alfa(n1,lm,idf))*gama(n2,jlo,l1mm1,idf)*&
     &                cmplx(iululol1(jlo,l,iat,isp),0.0d0,8)+               &
     &               conjg(beta(n1,lm,idf))*gama(n2,jlo,l1mm1,idf)*     &
     &                cmplx(iudlulol1(jlo,l,iat,isp),0.0d0,8)
               
                xmy2 = xmy1 + &
     &               conjg(alfa(n1,lm,idf))*gama(n2,jlo,l1m1,idf)* &
     &                cmplx(iululol1(jlo,l,iat,isp),0.0d0,8)+                &
     &               conjg(beta(n1,lm,idf))*gama(n2,jlo,l1m1,idf)*      &
     &                cmplx(iudlulol1(jlo,l,iat,isp),0.0d0,8)

                z2 = z2 + &
     &               conjg(alfa(n1,lm,idf))*gama(n2,jlo,l1m,idf)*  &
     &                cmplx(iululol1(jlo,l,iat,isp),0.0d0,8)+                &
     &               conjg(beta(n1,lm,idf))*gama(n2,jlo,l1m,idf)*       &
     &                cmplx(iudlulol1(jlo,l,iat,isp),0.0d0,8)
              enddo
!
!             the loop for the Local orbitals on the left side 
!
              do ilo = 1,nLO_at(l+1,iat)
                xpy1 = xpy1 +  &
     &               conjg(gama(n1,ilo,l1m1,idf))*alfa(n2,lm,idf)*    &
     &                cmplx(iulol1ul(ilo,l,iat,isp),0.0d0,8)+         &
     &               conjg(gama(n1,ilo,l1m1,idf))*beta(n2,lm,idf)*      &
     &                cmplx(iulol1udl(ilo,l,iat,isp),0.0d0,8)

                xmy1 = xmy1 + &
     &               conjg(gama(n1,ilo,l1mm1,idf))*alfa(n2,lm,idf)*     &
     &                cmplx(iulol1ul(ilo,l,iat,isp),0.0d0,8)+                &
     &               conjg(gama(n1,ilo,l1mm1,idf))*beta(n2,lm,idf)*     &
     &                cmplx(iulol1udl(ilo,l,iat,isp),0.0d0,8) 

                z1 = z1 + &
     &               conjg(gama(n1,ilo,l1m,idf))*alfa(n2,lm,idf)*     &
     &                cmplx(iulol1ul(ilo,l,iat,isp),0.0d0,8)+         &
     &               conjg(gama(n1,ilo,l1m,idf))*beta(n2,lm,idf)*     &
     &                cmplx(iulol1udl(ilo,l,iat,isp),0.0d0,8) 
              enddo

              do ilo = 1,nLO_at(l,iat)
                xpy2 = xpy2 +  &
     &               conjg(gama(n1,ilo,lm,idf))*alfa(n2,l1mm1,idf)*   &
     &                cmplx(iulolul1(ilo,l,iat,isp),0.0d0,8) +        &
     &               conjg(gama(n1,ilo,lm,idf))*beta(n2,l1mm1,idf)*   &
     &                cmplx(iuloludl1(ilo,l,iat,isp),0.0d0,8)

                xmy2 = xmy2 + &
     &               conjg(gama(n1,ilo,lm,idf))*alfa(n2,l1m1,idf)*    &
     &                cmplx(iulolul1(ilo,l,iat,isp),0.0d0,8)+         &
     &               conjg(gama(n1,ilo,lm,idf))*beta(n2,l1m1,idf)*    &
     &                cmplx(iuloludl1(ilo,l,iat,isp),0.0d0,8)

                z2 = z2 + &
     &               conjg(gama(n1,ilo,lm,idf))*alfa(n2,l1m,idf)*     &
     &                cmplx(iulolul1(ilo,l,iat,isp),0.0d0,8)+         &
     &               conjg(gama(n1,ilo,lm,idf))*beta(n2,l1m,idf)*     &
     &                cmplx(iuloludl1(ilo,l,iat,isp),0.0d0,8)
              enddo
    
              do jlo = 1,nLO_at(l,iat)
                do ilo = 1, nLO_at(l+1,iat) 
                  xpy1 = xpy1 +  &
     &             conjg(gama(n1,ilo,l1m1,idf))*gama(n2,jlo,lm,idf)*  &
     &              cmplx(iulol1ulol(ilo,jlo,l,iat,isp),0.0d0,8)

                  xmy1 = xmy1 + &
     &             conjg(gama(n1,ilo,l1mm1,idf))*gama(n2,jlo,lm,idf)* &
     &                cmplx(iulol1ulol(ilo,jlo,l,iat,isp),0.0d0,8)

                  z1 = z1 + &
     &             conjg(gama(n1,ilo,l1m,idf))*gama(n2,jlo,lm,idf)*   &
     &                cmplx(iulol1ulol(ilo,jlo,l,iat,isp),0.0d0,8)

                  xpy2 = xpy2 +  &
     &             conjg(gama(n1,jlo,lm,idf))*gama(n2,ilo,l1mm1,idf)* &
     &                cmplx(iulolulol1(jlo,ilo,l,iat,isp),0.0d0,8)

                  xmy2 = xmy2 + &
     &             conjg(gama(n1,jlo,lm,idf))*gama(n2,ilo,l1m1,idf)*  &
     &                cmplx(iulolulol1(jlo,ilo,l,iat,isp),0.0d0,8)

                  z2 = z2 + &
     &               conjg(gama(n1,jlo,lm,idf))*gama(n2,ilo,l1m,idf)* &
     &                cmplx(iulolulol1(jlo,ilo,l,iat,isp),0.0d0,8)

                enddo
              enddo 

              pxpy = pxpy + flm1 * xpy1 + flm3 * xpy2
              pxmy = pxmy + flm3 * xmy1 + flm1 * xmy2
              pz   = pz   + flm5 * ( z1 + z2 )
            
            enddo ! m
          enddo ! l
          mxpy = mxpy + pxpy
          mxmy = mxmy + pxmy
          mz = mz + pz
         
        enddo ! ieq
      enddo ! iat
      mt12(1)=-5.0d-1*imag*(mxpy+mxmy)
      mt12(2)=5.0d-1*(mxmy-mxpy)
      mt12(3)=-imag*mz

!
! Interstitial contributions 
!

      ikv(1:3)=klist(1:3,ik)
      irk=kpirind(ik)
      call k2cart(ikv,idvk,kv)
      mpw=czero
      do ig2=1,nv(irk)
        igv2(1:3)=gindex(:,indgk(ig2,ik))
        call k2cart(igv2,1,gv2)
        kpgv2(1:3)=gv2(1:3)+kv(1:3)
        intmom = czero
        do ig1=1,nv(irk)
          igv1(1:3)=gindex(1:3,indgk(ig1,ik))
          iv1=igv2(1)-igv1(1)
          iv2=igv2(2)-igv1(2)
          iv3=igv2(3)-igv1(3)
          intmom=intmom+conjg(zzk(ig1,n1))*ipwint(ig0(iv1,iv2,iv3))
        enddo ! ig1
        modmom=zzk(ig2,n2)*intmom
        mpw=mpw+kpgv2*modmom
      enddo ! ig2

!
! Sum of MT and IS contributions 
!

      p12=mt12+mpw

      end subroutine calcmmatvv
!EOC      
