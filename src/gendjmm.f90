!BOP
!
! !ROUTINE: gendjmm
!
! !INTERFACE: 
      subroutine gendjmm
      
! !DESCRIPTION:
!
! Generates the rotation matrices for spherical harmonics $D^j_{mm'}$ up
!to \texttt{maxbigl} once
! $D^1_{mm'}$ is known using the inverse Clebsch-Gordan series:
!
!\begin{equation}
!D^j_{\mu m}=\sum\limits_{\mu_1=-1}^1{\sum\limits_{m_1=1}^1{C(1,j-1,j;%
!\mu_1,\mu-\mu_1)C(1,j-1,j;m_1,m-m_1)D^1_{\mu_1m_1}D^{j-1}_{\mu-\mu_1,m-m_1}}}
!\end{equation}.
!
! Since the program already calculates the Gaunt coefficients
!$G^{LM}_{l_1m_1,l_2m_2}$ we replace the product of Clebsch-Gordan
!coefficients using:
!
! \begin{equation}
!C(1,j-1,j;\mu_1,\mu-\mu_1)C(1,j-1,j;m_1,m-m_1)=\sqrt{\frac{4\pi(2j+1)}%
!{3(2j-1)}}\frac{G^{j\mu}_{1\mu_1,j-1\mu-\mu_1}G^{jm}_{1m_1,j-1m-m_1}}%
!{G^{j0}_{10,j-10}}
!\end{equation}
!
! !USES:

      use constants, only: czero,pi
      use mixbasis,  only: maxbigl
      use rotylm,    only: djmm
      use struk,     only: nat,ndf,mult,rotij, rotloc


! !LOCAL VARIABLES:
      
      implicit none
      
      integer(4) :: dimdj ! second dimension of the array djmm
      integer(4) :: iat   ! indexes inequivalent atoms
      integer(4) :: idf   ! indexes all atoms
      integer(4) :: idx   ! index of 2nd. dimension of array djmm
      integer(4) :: ieq   ! indexes equivalent atoms
      integer(4) :: i,j,k ! couters
      integer(4) :: l     ! angular momentum quantum number
      integer(4) :: mu,mu1,m,m1,mu2,m2 ! azimutal quantum numbers
      
      real(8) :: cg1,cg2  ! Gaunt coefficients
      real(8) :: prefac   ! prefactor for djmm =sql/gaunt(1,l-1,l,0,0,0)
      real(8) :: sql      ! sqrt(4pi(2l+1)/3(2l-1)
      real(8), dimension(3,3) :: rotcart ! local rotation matrix in the
!                                          cartesian basis
      
      complex(8) :: dj1   ! D^1_mu1m1
      complex(8) :: dj2   ! D^(l-1)_(mu-mu1,m-m1) 
      complex(8), dimension(3,3) :: rotsph ! local rotation matrix in the
!                                          spherical basis
      
      
 
!
! !EXTERNAL ROUTINES: 

      real(8), external :: getcgcoef
      
      complex(8), external :: getdjmm
      
      external rot1tosph

! !INTRINSIC ROUTINES: 

      intrinsic dsqrt


! !REVISION HISTORY:
!
! Created 10th. August 2004 by RGA
!
!EOP
!BOC
!
!     Calculate the dimension of djmm
!
      dimdj=(maxbigl+1)*(maxbigl+2)*(4*maxbigl+3)/6
!
!     Allocate djmm and initialize it to zero
!
      allocate(djmm(ndf,dimdj))
      djmm=czero
!
!     Loop over all atoms
!     
      idf=0
      do iat=1,nat
        do ieq=1,mult(iat)
          idf=idf+1
!
!         D^0_00 =0
!          
          djmm(idf,1)=1.0d0
!
!         Calculate the rotation matrix in the cartesian basis
!         rotcart=transpose(rotloc x rotij)
!         
          rotcart(1:3,1:3)=0.0d0
          do i=1,3
            do j=1,3
              do k=1,3
                rotcart(i,j)=rotcart(i,j)+rotloc(i,k,iat)*rotij(k,j,idf)
              enddo
            enddo
          enddo


!
!         Transform the rotation matrix to the spherical basis
!          
          call rot1tosph(rotcart,rotsph)
!
!         Obtain the rotation matrix for j=1 D^1_{m,m')=transpose(rotsph)
!           
          idx=1
          do i=-1,1
            do j=0,1
              idx=idx+1
              djmm(idf,idx)=rotsph(j+2,i+2)
            enddo
          enddo
!
!         Obtain the rotation matrix for j> 1 by recurrence          
      
          do l=2,maxbigl    
            sql=dsqrt(4.0d0*pi*dble(2*l+1)/dble(2*l-1)/3.0d0)
            prefac=sql/getcgcoef(1,l-1,l,0,0)
            do mu=-l,l
              do mu1=-1,1
                mu2=mu-mu1
                if(abs(mu2).le.l-1)then
                  cg1=getcgcoef(1,l-1,l,mu1,mu2)
                  do m=0,l
                    idx=l*(l+1)*(4*l-1)/6+(l+1)*(mu+l)+m+1
                    do m1=-1,1
                      m2=m-m1
                      if(abs(m2).le.l-1)then
                        dj1=getdjmm(idf,1,mu1,m1)
                        dj2=getdjmm(idf,l-1,mu2,m2)
                        cg2=getcgcoef(1,l-1,l,m1,m2)
                        djmm(idf,idx)=djmm(idf,idx)+cg1*cg2*dj1*dj2
                      endif
                    enddo ! m1
                  enddo ! m
                endif  
              enddo ! mu1
            enddo ! mu
            do mu=-l,l
              do m=0,l
                idx=l*(l+1)*(4*l-1)/6+(l+1)*(mu+l)+m+1
                djmm(idf,idx)=prefac*djmm(idf,idx)
              enddo
            enddo    
          enddo ! l
        enddo ! ieq
      enddo ! iat    
      end subroutine gendjmm
!EOC                
      
