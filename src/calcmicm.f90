!BOP
!
! !ROUTINE: calcmicm
!
! !INTERFACE:
      subroutine calcmicm(ik,iq,cstart,cend,mstart,mend,isp,micm)

! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $M^i_{nm}(\vec{k},\vec{q})$ 
!of equation \ref{Mdef} for n => core and m => LAPW vector
!
! !USES:

      use bands,       only: nspin 
      use barcoul,     only: barcvm,iop_coul,im_g0
      use constants,   only: czero, cone, pi
      use core,        only: ncore,ncg,corind,core_ortho
      use eigenvec,    only: alfa, beta, gama, alfp, betp, gamp
      use kpoints,     only: klist,idvk,kqid
      use lapwlo,      only: nt,lomax,nlo_at 
      use mixbasis,    only: nmix,bigl,locmatsiz,matsiz,locmixind,mbsiz,&
     &                       s3r,lo_index_s3r
      use struk,       only: mult,vi,pos,nat
      use task,        only: time_lapack,time_minm

! !INPUT PARAMETERS:

      implicit none

      integer, intent(in) :: ik  ! the order index of kvec
      integer, intent(in) :: iq
      integer, intent(in) :: cstart,cend,mstart,mend
      integer, intent(in) :: isp  ! index for spin
      complex(8), intent(out):: micm(matsiz,mstart:mend,ncg)


! !LOCAL VARIABLES:

      integer :: bl      ! Angular momentum quantum number of the mixed atomic functions (big l)
      integer :: bm      ! z-component angular momentum quantum number of the mixed atomic function (big m)
      integer :: i 
      integer :: iat     ! Counter: runs over inequivalent atoms.
      integer :: icg,idf,idf0     ! Counter: runs over atoms, including equivalent ones
      integer :: ic1 
      integer :: ie12,ie1,ie2 ! Counter: run over the eigenvalues in  the corresponding k-point.
      integer :: ieq     ! Counter: runs over the equivalent atoms of each type.
      integer :: im,imix    ! Counter: runs over all MT-sphere mixed basis functions (of all the atoms)
      integer :: irm     ! Counter: runs over the radial mixed basis functions of each atom.
      integer :: ilo,irlo
      integer :: l1      ! l-quantum number of the (L)APW eigenfunction at k
      integer :: il2,l2  ! l-quantum number of the (L)APW eigenfunction at k' = k - q
      integer :: m1      ! m-quantum number of the (L)APW eigenfunction at k 
      integer :: m2      ! m-quantum number of the (L)APW eigenfunction at k' = k - q
      integer :: l2min   ! Lowest allowed value of lambda = | bl - l2 |
      integer :: l2max   ! Highest allowed value of lambda = bl + l2
      integer :: l1m1    ! Counter: Runs over MTS (L)APW basis  functions (l1,m1) 
      integer :: l2m2    ! Counter: Runs over MTS (L)APW basis   functions l2,m2)
      integer :: cmdim

      integer, dimension(3) :: igc ! Indexes of G_1+G'-G
      integer ::ierr

      real(8) :: sqvi,x, arg
      real(8) :: kqvec1(3),kqvec2(3)
      real(8) :: time1,time2,tstart,tend

      complex(8) :: angint ! Angular integral calculated by fourylm
      complex(8) :: suml12,suml2m2
      complex(8) :: phs
      complex(8) :: sumterms
      complex(8) :: aa  ! the term A^{a*}_{l1,m1}(k-q+K')*A^a_{\l2\m2}(k+G)<bl,l2,l1|l2> in equation \ref{dbracketorg} 
      complex(8) :: ab  ! the term A^{a*}_{l1,m1}(k-q+K')*B^a_{l2,m2}(k+G)<bl,l2,l1|l2> in equation \ref{dbracketorg}
      complex(8) :: ac  ! the term A^{a*}_{l1,m1}(k-q+K')*C^a_{l2,m2}(k+G)<bl,l2,l1|l2>in equation \ref{dbracketorg}
      character(8), parameter :: sname = 'calcmicm'
      logical :: lprt =  .false.

      complex(8),allocatable::mmat(:,:,:)
 
!
! !EXTERNAL ROUTINES: 


      real(8), external :: getcgcoef


! !INTRINSIC ROUTINES: 


      intrinsic abs
      intrinsic min
      intrinsic nint

! !REVISION HISTORY:
! 
! Created  23th. Feb. 2004 by RGA
! Last modified on Oct. 2, 2013 by JH

!EOP
!
!BOC

      call cpu_time(tstart)
      if(lprt) call linmsg(6,'-','calcmicm')

      cmdim=ncg*(mend-mstart+1)
      allocate(mmat(locmatsiz,mstart:mend,ncg),stat=ierr)
      call errmsg(ierr.ne.0,sname,'Fail to allocate mmat')
      mmat=czero

      do i=1,3
        kqvec1(i)=-dble(klist(i,ik))/dble(idvk)
        kqvec2(i)= dble(klist(i,kqid(ik,iq)))/dble(idvk)
      enddo

      idf0=0
      do icg = 1, ncg             !! loop over core states 
        iat=corind(1,icg)
        idf=corind(2,icg)
        arg=sum(pos(:,idf)*kqvec2(:))
        phs=cmplx(cos(2.0d0*pi*arg),sin(2.0d0*pi*arg),8)
        if(idf.ne.idf0)then        !! reset l1m1 if it is a new atom
          l1m1=0
          idf0=idf
        endif  
        ic1=corind(3,icg)
        l1=corind(4,icg)
        m1=corind(5,icg)
        l1m1=l1m1+1

        if(icg.lt.cstart.or.icg.gt.cend) cycle 

        do ie2 = mstart, mend

          imix=0
          do irm = 1, nmix(iat)      !! Loop over mixed functions
            bl=bigl(irm,iat)

            do bm=-bl,bl            !! Loop over bm
              imix = imix + 1
              im=locmixind(imix,idf)
              l2min=iabs(bl-l1)
              l2max=min(bl+l1,nt-1)
              suml2m2 = czero

              do l2=l2min,l2max

                m2=bm+m1
                if(iabs(m2).gt.l2) cycle 

                l2m2=l2*l2+l2+m2+1
                angint=getcgcoef(l1,bl,l2,m1,bm)    ! ! Calculate the angular integral
                if(abs(angint).le.1.0d-8) cycle 

                il2= l2 + ncore(iat)+1
                aa=s3r(irm,ic1,il2,iat,isp)*alfp(ie2,l2m2,idf)
                ab=s3r(irm,ic1,il2+nt,iat,isp)*betp(ie2,l2m2,idf)

                ac = 0.d0 
                do ilo=1,nLO_at(l2,iat)
                  irlo = lo_index_s3r(ilo,l2,iat) 
                  ac = ac + s3r(irm,ic1,irlo,iat,isp)*gamp(ie2,ilo,l2m2,idf)
                enddo 
                suml2m2 = suml2m2 + angint * (aa + ab + ac) 
              enddo ! l2

              mmat(im,ie2,icg) = suml2m2*phs
            enddo   ! bm
          enddo  ! irm 
        enddo  ! icg
      enddo ! ie2

!
! Transform to the eigenvectors of the coulomb matrix
!
      call cpu_time(time1)
      call zgemm('c','n',matsiz,cmdim,locmatsiz,cone,barcvm,mbsiz,mmat, &
     &           locmatsiz,czero,micm,matsiz)
      call cpu_time(time2)
      time_lapack=time_lapack+time2-time1
      deallocate(mmat)
      call cpu_time(tend)

      time_minm = time_minm + tend - tstart 
      if(iop_coul.gt.0.and.iq.eq.1) then
        do icg=cstart,cend
          do ie2=mstart,mend
            micm(im_g0,ie2,icg) = 0.d0
          enddo
        enddo
      endif


      end subroutine calcmicm
!EOC      


