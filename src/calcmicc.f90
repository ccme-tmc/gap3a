!BOP
!
! !ROUTINE: calcmicc
!
! !INTERFACE:
      subroutine calcmicc(isp,ccdim,micc)

! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $M^i_{nm}(\vec{k},\vec{q})$ 
!of equation \ref{Mdef} for n => core and m => LAPW vector
!
! !USES:

      use bands,       only: nspin 
      use barcoul,     only: barcvm
      use constants,   only: czero, cone, pi
      use core,        only: ncore,ncg,corind,clmind,nclmmax,nclm
      use kpoints,     only: klist,idvk,kqid
      use mixbasis,    only: nmix,bigl,lmixmax,locmatsiz,locmixind,&
     &                       mbsiz,matsiz,s3r
      use struk,       only: mult,vi,pos,nat,ndf,inddf 
      use task,        only: time_lapack

! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: isp
      integer(4), intent(in) :: ccdim
      complex(8), intent(out):: micc(1:matsiz,ccdim)

! !LOCAL VARIABLES:

      integer(4) :: bl         ! Angular momentum quantum number of the mixed atomic functions (big l)
      integer(4) :: bm         ! z-component angular momentum quantum number of the mixed atomic function (big m)
      integer(4) :: iat        ! index for inequivalent atoms.
      integer(4) :: idf        ! index for atoms, including equivalent ones
      integer(4) :: icg1, icg2 ! index for core states for all atoms 
      integer(4) :: iclm1,iclm2 ! index for core states at each atom including (lm) 
      integer(4) :: ic1,ic2     ! index for core states at each atom 
      integer(4) :: ic12

      integer(4) :: im      ! index for MT mixed basis functions (of all the atoms)
      integer(4) :: imix    ! index for mixed basis functions for each atom 
      integer(4) :: irm     ! index the radial mixed basis functions of each atom.
      integer(4) :: l1      ! l-quantum number of the (L)APW eigenfunction at k
      integer(4) :: il2,l2  ! l-quantum number of the (L)APW eigenfunction at k' = k - q
      integer(4) :: m1      ! m-quantum number of the (L)APW eigenfunction at k 
      integer(4) :: m2      ! m-quantum number of the (L)APW eigenfunction at k' = k - q
      integer(4) :: l2min   ! Lowest allowed value of lambda = | bl - l2 |
      integer(4) :: l2max   ! Highest allowed value of lambda = bl + l2
      integer(4) :: l2m2    ! Counter: Runs over MTS (L)APW basis   functions l2,m2)
      integer(4) :: la1,lb1,lc1,la2,lb2,lc2,i
      integer(4) :: ierr

      complex(8) :: angint ! Angular integral calculated by fourylm
      complex(8) :: aaterm  ! the term A^{a*}_{l1,m1}(k-q+K')*A^a_{\l2\m2}(k+G)<bl,l2,l1|l2> in equation \ref{dbracketorg} 
 
      complex(8), allocatable:: mcc(:,:)

      real(8) :: time1,time2
      logical :: lprt = .false.

!
! !EXTERNAL ROUTINES: 


      real(8), external :: getcgcoef


! !INTRINSIC ROUTINES: 


      intrinsic abs
      intrinsic min
      intrinsic nint

! !REVISION HISTORY:
! 
! Created  March 21, 2007 by JH

!EOP
!
!BOC
      if(lprt) call linmsg(6,'-','calcmicc') 

      allocate(mcc(locmatsiz,ccdim))
      mcc = czero
      
      ic12=0
      do idf=1,ndf 
        iat=inddf(idf)
        do iclm2 = 1, nclm(idf)               !! loop over core states 
          icg2 = clmind(iclm2,idf) 
          ic2  = corind(3,icg2) 
          l2   = corind(4,icg2) 
          m2   = corind(5,icg2)
 
          do iclm1 = 1, nclm(idf) 
            icg1= clmind(iclm1,idf) 
            ic1 = corind(3,icg1) 
            l1  = corind(4,icg1) 
            m1  = corind(5,icg1) 
            ic12=ic12+1

            imix=0
            do irm = 1, nmix(iat)      !! Loop over mixed functions
              bl=bigl(irm,iat)
              do bm=-bl,bl            !! Loop over bm
                imix = imix + 1
                im=locmixind(imix,idf)
                angint=getcgcoef(l1,bl,l2,m1,bm)    ! ! Calculate the angular integral
                aaterm=s3r(irm,ic1,ic2,iat,isp)*angint
                mcc(im,ic12)=aaterm
              enddo ! bm
            enddo   ! irm
          enddo  ! iclm1 
        enddo  ! iclm2
      enddo ! idf

!
! Transform to the eigenvectors of the coulomb matrix
!
      call cpu_time(time1)
      call zgemm('c','n',matsiz,ccdim,locmatsiz,cone,barcvm,mbsiz,mcc, &
     &           locmatsiz,czero,micc,matsiz)
      call cpu_time(time2)
      time_lapack=time_lapack+time2-time1

      end subroutine calcmicc
!EOC      


