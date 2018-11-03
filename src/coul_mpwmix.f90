!BOP
!
! !ROUTINE: coul_mpwmix
!
! !INTERFACE:
      subroutine coul_mpwmix(iq)
      
! !DESCRIPTION:
!
! This subroutine calculates the matrix elements between mixed basis
! functions and plane waves.
!
! !USES:

      use constants,  only: czero,imag, pi
      use kpoints,    only: qlist, idvq
      use mixbasis,   only: bigl, maxbigl, mbl, nmix, nmixmax,locmatsiz,&
     &                      mbsiz,mpwipw, mpwmix, wi0
      use recipvec,   only: gindex, indgq, indgqlen, ngqbarc,maxngqlen
      use barcoul,    only: jlam
      use struk,      only: mult, pos, vi,nat,rotij, rotloc,dh, nrpt
      
! !INPUT PARAMETERS:      
      implicit none
      
      integer, intent(in) :: iq
!
!
! !LOCAL VARIABLES:

      integer :: ipin
      integer :: iat  ! (Counter): runs over inequivalent atoms
      integer :: i1   ! (Counter): runs over local mixed basis functions
      integer :: irm  ! (Counter): runs over radial mixed basis functions
      integer :: ieq  ! (Counter): runs over equivalent atoms 
      integer :: idf  ! (Counter): runs over equivalent and inequivalent atoms
      integer :: igl
      integer :: ippw ! (Counter): runs over plane waves
      integer :: l1   ! angular momentum quantum number of the mixed  basis function
      integer :: m1   ! z component of the angular momentum of the mbf
      integer :: lm1
      integer :: ylsize 
      integer, dimension(3) :: iqvec,igvec  ! integer coodintates of the G-vector
      
      real(8) :: dd ! exponential step of the radial mesh
      real(8) :: gpr ! Scalar product G.pos(atom)
      
      real(8), dimension(3) :: qvec ! the q-vector which is zero
      real(8), dimension(3) :: gvec ! the G-vector 
      real(8), dimension(3) :: qg1 ! the G-vector half rotated
      real(8), dimension(3) :: qg2 ! the G-vector half rotated
      real(8), dimension(3) :: qg3 ! the G-vector rotated
      
      complex(8) :: prefac
      complex(8) :: y((maxbigl+1)*(maxbigl+1))! spherical harmonics

      complex(8), allocatable :: expg(:) ! Exp(gpr)
      complex(8), allocatable :: sph(:,:)

 
! !EXTERNAL ROUTINES: 


      external calcjlam
      external k2cart
      external rotate
      external ylm


! !INTRINSIC ROUTINES: 


      intrinsic sqrt
      intrinsic exp
      intrinsic cos
      intrinsic sin

!
! !REVISION HISTORY:
! 
! Created: 23. Nov. 2004 by RGA

!EOP
!BOC

!
!      Local functions:
!
      ylsize=(maxbigl+1)*(maxbigl+1)
      allocate(expg(ngqbarc(iq)))
      allocate(sph(ylsize,ngqbarc(iq)))
      allocate(jlam(nmixmax,maxngqlen)) 
      iqvec(1:3)=qlist(1:3,iq)
      call k2cart(iqvec,idvq,qvec)
      
      if(iq.eq.1)then
        mpwmix(1:mbsiz,1)=wi0(1:mbsiz)
        ipin=2
      else
        ipin=1
      endif    
!
!     Loop over inequivalent atoms:
!
      idf=0
      i1=0
      do iat = 1, nat
        dd=exp(dh(iat))
        call calcjlam(iat,mbl(iat),nrpt(iat),iq)

        do ieq = 1, mult(iat)
          idf=idf+1
!
!     Calculate Y_lm(q+G) for all G
!
          if(ipin.eq.2) then 
            sph(:,1) = 0.d0
            expg(1) = 0.d0
            mpwmix(:,1) = 0.d0
          endif 

          do ippw=ipin,ngqbarc(iq)
            igvec(1:3)=gindex(:,indgq(ippw,iq))
            call k2cart(igvec,1,gvec)
            qg1(1:3)=gvec(1:3)+qvec(1:3)
            call rotate(qg1,rotij(1:3,1:3,idf),qg2)
            call rotate(qg2,rotloc(1:3,1:3,iat),qg3)
            call ylm(qg3,maxbigl,y)
            sph(1:ylsize,ippw)=y(1:ylsize)
            gpr=dble(igvec(1))*pos(1,idf)+ &
     &          dble(igvec(2))*pos(2,idf)+ &
     &          dble(igvec(3))*pos(3,idf)
            expg(ippw)=cmplx(cos(2.0d0*pi*gpr),sin(2.0d0*pi*gpr),8)
          enddo
          
          do irm = 1, nmix(iat)
            l1=bigl(irm,iat)
            prefac=cmplx(4.0d0*pi*sqrt(vi),0.0d0,8)*imag**l1
            do m1=-l1,l1
              i1 = i1 + 1
              lm1=l1*l1+l1+m1+1
              do ippw=ipin,ngqbarc(iq)
                mpwmix(i1,ippw)=czero
                igl=indgqlen(ippw,iq)
                mpwmix(i1,ippw)= prefac*expg(ippw)*jlam(irm,igl) &
     &                          *conjg(sph(lm1,ippw))
              enddo ! ippw
            enddo ! m1  
          enddo ! irm
        
        enddo ! ieq
      enddo ! iat    
            
      do ippw=1,ngqbarc(iq)
        do i1=locmatsiz+1,mbsiz
          mpwmix(i1,ippw)=mpwipw(i1-locmatsiz,ippw)
        enddo
      enddo  
      deallocate(sph,expg,jlam)

      end subroutine coul_mpwmix
!EOC
