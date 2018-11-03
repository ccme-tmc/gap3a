!BOP
!
! !ROUTINE: calcminc
!
! !INTERFACE:
      subroutine calcminc(ik,iq,nstart,nend,cstart,cend,isp,minc)

! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $M^i_{nm}(\vec{k},\vec{q})$ 
!of equation \ref{Mdef} for 
!      n => lapw states, m  => core states
!
!
! !USES:
 
      use bands,       only: nspin 
      use barcoul,     only: barcvm,iop_coul,im_g0 
      use constants,   only: czero, cone, pi
      use core,        only: ncore,ncg,corind
      use eigenvec,    only: alfa, beta, gama
      use kpoints,     only: klist,idvk,kqid
      use lapwlo,      only: lomax,lmax,nlo_at,nt
      use mixbasis,    only: mbsiz,nmix,bigl,locmatsiz,matsiz,locmixind,&
     &                       s3r,lo_index_s3r
      use struk,       only: mult,vi,pos,nat
      use task,        only: time_lapack,time_minm

! !INPUT PARAMETERS:

      implicit none

      integer, intent(in) :: ik  ! the order index of kvec
      integer, intent(in) :: iq
      integer, intent(in) :: nstart,nend
      integer, intent(in) :: cstart,cend
      integer, intent(in) :: isp
      complex(8), intent(out):: minc(matsiz,cstart:cend,nstart:nend)


! !LOCAL VARIABLES:

      integer :: bl      ! Angular momentum quantum number of the mixed atomic functions (big l)
      integer :: bm      ! z-component angular momentum quantum number of the mixed atomic function (big m)
      integer :: i 
      integer :: iat     ! Counter: runs over inequivalent atoms.
      integer :: icg,idf,idf0     ! Counter: runs over atoms, including equivalent ones
      integer :: ic      ! index for the core state on a atom ( m not counted) 
      integer :: ie1     ! Counter: run over the eigenvalues in  the corresponding k-point.
      integer :: ieq     ! Counter: runs over the equivalent atoms of each type.
      integer :: im,imix    ! Counter: runs over all MT-sphere mixed basis functions (of all the atoms)
      integer :: irm     ! Counter: runs over the radial mixed basis functions of each atom.
      integer :: l1      ! l-quantum number of the (L)APW eigenfunction at k
      integer :: il2,l2  ! l-quantum number of the (L)APW eigenfunction at k' = k - q
      
      integer :: m1      ! m-quantum number of the (L)APW eigenfunction at k 
      integer :: m2      ! m-quantum number of the (L)APW eigenfunction at k' = k - q
      integer :: l2min   ! Lowest allowed value of lambda = | bl - l2 |
      integer :: l2max   ! Highest allowed value of lambda = bl + l2
      integer :: ilo, irlo ! index for LO orbitals 
      integer :: l1m1    ! Counter: Runs over MTS (L)APW basis  functions (l1,m1) 
      integer :: l2m2    ! Counter: Runs over MTS (L)APW basis   functions l2,m2)
      integer :: cdim,ndim,ncdim

      integer ::ierr
      
      real(8) :: sqvi, x, arg
      real(8) :: qvec(3),kqvec1(3)
      real(8) :: time1,time2,tstart,tend

      complex(8) :: angint ! Angular integral calculated by fourylm
      complex(8) :: suml12,suml2m2
      complex(8) :: phs
      complex(8) :: sm
      complex(8) :: aa  ! the term A^{a*}_{l1,m1}(k-q+K')*A^a_{\l2\m2}(k+G)<bl,l2,l1|l2> in equation \ref{dbracketorg} 
      complex(8) :: bb  ! the term B^{a*}_{l1,m1}(k-q+K')*B^a_{l2,m2}(k+G)<bl,l2,l1|l2> in equation \ref{dbracketorg}
      complex(8) :: cc  ! the term C^{a*}_{l1,m1}(k-q+K')*C^a_{l2,m2}(k+G)<bl,l2,l1|l2> in equation \ref{dbracketorg}

      character(8), parameter :: sname = 'calcminc'
      logical :: notf,lprt=.false.
      logical :: lcont  !! define whether minc to calculated in this subroutine is contiguous in memory

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
! Last modified Jan. 16, 2007 by JH

! !TO DO:
!
! - Check complex and real cases
!
!EOP
!
!BOC

      if(lprt) call linmsg(6,'-','calcminc')
      call cpu_time(tstart)

      ndim = nend-nstart+1
      cdim = cend-cstart+1
      ncdim = cdim*ndim
      allocate(mmat(locmatsiz,cstart:cend,nstart:nend),stat=ierr)
      call errmsg(ierr.ne.0,sname,'Fail to allocate mmat')
     
      do i=1,3
        kqvec1(i)=-dble(klist(i,ik))/dble(idvk)
      enddo  

      mmat=czero
      do ie1=nstart, nend

        idf0=0
        do icg = 1, ncg             !! loop over core states
          iat=corind(1,icg)
          idf=corind(2,icg)

          arg=sum(pos(:,idf)*kqvec1(:))
          phs=cmplx(cos(2.0d0*pi*arg),sin(2.0d0*pi*arg),8)

          if(idf.ne.idf0)then        !! reset l1m1 if it is a new atom
            l1m1=0
            idf0=idf
          endif

          ic = corind(3,icg)
          l1 = corind(4,icg)
          m1 = corind(5,icg)
          l1m1 = l1m1+1

          if(icg.lt.cstart.or.icg.gt.cend) cycle 

          imix=0
          do irm = 1, nmix(iat)      !! Loop over mixed functions
            bl=bigl(irm,iat)
            do bm=-bl,bl            !! Loop over bm
              imix = imix + 1
              im = locmixind(imix,idf)

              l2min = iabs(bl-l1)
              l2max = min(bl+l1,lmax)

              suml2m2 = czero
              do l2=l2min,l2max   ! loop over l2
                m2=m1-bm

                if(iabs(m2).gt.l2) cycle 

                l2m2=l2*l2+l2+m2+1

                angint=getcgcoef(l2,bl,l1,m2,bm)  !Calculate the angular integral:

                if(abs(angint).lt.1.0d-8) cycle 

                il2= l2 + ncore(iat) + 1
                aa = s3r(irm,ic,il2,iat,isp)*alfa(ie1,l2m2,idf)
                bb = s3r(irm,ic,il2+nt,iat,isp)*beta(ie1,l2m2,idf)
                cc = 0.d0
                do ilo=1,nLO_at(l2,iat)
                  irlo = lo_index_s3r(ilo,l2,iat) 
                  cc=cc+ s3r(irm,ic,irlo,iat,isp)*gama(ie1,ilo,l2m2,idf)
                enddo
                sm = aa + bb + cc
                suml2m2 = suml2m2 + angint * sm
              enddo ! l2
              mmat(im,icg,ie1) = suml2m2*phs
            enddo !bm  
          enddo ! irm  
        enddo ! ie1
      enddo ! icg

!
! Transform to the eigenvectors of the coulomb matrix
!
      call cpu_time(time1)     
      call zgemm('c','n',matsiz,ncdim,locmatsiz,cone,barcvm,mbsiz,    &
     &             mmat,locmatsiz,czero,minc,matsiz)
      call cpu_time(time2) 
      time_lapack=time_lapack+time2-time1    
      deallocate(mmat)

      if(iop_coul.gt.0.and.iq.eq.1) then
        do ie1=nstart,nend
          do icg=cstart,cend
            minc(im_g0,icg,ie1) = 0.d0
          enddo
        enddo
      endif

      call cpu_time(tend) 
      time_minm=time_minm+tend-tstart    
      end subroutine calcminc
!EOC      


