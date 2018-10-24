!BOP
! !ROUTINE: diagipw
! !INTERFACE:
      subroutine diagipw(iq)

! !DESCRIPTION:
!
! This subroutine generates the overlap matrix of the product functions in the
! interstitial region and
! diagonalizes it. The output is the matrix $S_{\vec{G}i}$.
!
! !USES:

      use constants,  only: pi
      use kpoints,    only: nkp
      use mixbasis,   only: ipwint, sgi
      use recipvec,   only: ngq, gindex, ig0, indgq
      use struk,      only: vi
      use task,       only: fid_outdbg
      
! !INPUT PARAMETERS:
       
      implicit none

      integer, intent(in) :: iq      

! !LOCAL VARIABLES:
      character*20:: sname="diagipw"
      integer :: i     ! Counter: Runs over cartesian coordinates
      integer :: iipw  ! Counter: Runs over iipw's
      integer :: jipw  ! Counter: Runs over iipw's
      integer :: info  ! Output status of the diagonalization routine.
      integer :: lwork ! Size of the workspace array work.
      integer, dimension(3) :: iik ! integer coordinates of G-G'
      integer :: ng

      real(8) :: pi2vi, alfa
      real(8), allocatable :: epsipw(:) ! the eigenvalues of sgi
      real(8), allocatable :: rwork(:)  ! Workspace array for diagonalization
      complex(8) :: cfact               ! Normalization factor
      complex(8), allocatable :: work(:)! Workspace array for diagonalization

      character(len=1), parameter :: jobz = 'v'
      character(len=1), parameter :: uplo = 'u'

 
! !EXTERNAL ROUTINES: 

      external zheev      ! Lapack diagonalization subroutine

! !INTRINSIC ROUTINES: 

      intrinsic dabs
      intrinsic dsqrt
!
! !REVISION HISTORY:
!
! Created Dec. 2003. by RGA
! Last Modification 10. Nov. 2005 by RGA
!
!EOP
! ----------------------------------------------------------------
!BOC

!
!     Allocate the working space needed by the diagonalizing subroutine
!
      lwork = 2*ngq(iq)-1
      allocate(work(lwork),rwork(3*ngq(iq)-2),epsipw(1:ngq(iq)))
!
!     Calculate the overlap matrix between product plane waves:
!
      sgi=0.d0
      do iipw=1,ngq(iq)
        sgi(iipw,iipw)=ipwint(1)
        do jipw=iipw+1,ngq(iq)
          iik(:)=gindex(:,indgq(iipw,iq))-gindex(:,indgq(jipw,iq))
          sgi(iipw,jipw)=ipwint(ig0(iik(1),iik(2),iik(3)))
          sgi(jipw,iipw)=conjg(sgi(iipw,jipw))
        enddo ! jipw
      enddo ! iipw
!
!     Diagonalize sgi:
!

#ifdef DEBUG
      write(fid_outdbg,*) "### sgi-0 ###"
      do jipw=1,ngq(iq),ngq(iq)/10
        do iipw=1,ngq(iq),ngq(iq)/10
          write(fid_outdbg,'(2i5,4f12.6)') iipw,jipw,sgi(iipw,jipw)
        enddo
      enddo
#endif 
      call zheev('V','U',ngq(iq),sgi,ngq(iq),epsipw,work,lwork,rwork,info)
      call errmsg0(info,sname,"Fail in calling zheev")

#ifdef DEBUG
      write(fid_outdbg,*) "### sgi-1 ###"
      do jipw=1,ngq(iq),ngq(iq)/10
        do iipw=1,ngq(iq),ngq(iq)/10
          write(fid_outdbg,'(2i5,4f12.6)') iipw,jipw,sgi(iipw,jipw)
        enddo
      enddo

      write(fid_outdbg,*)
      write(fid_outdbg,*) "### epsipw ###"
      do iipw=1,ngq(iq)
        write(fid_outdbg,'(i5,2f12.6)') iipw,epsipw(iipw)
      enddo 
#endif

!     Normalize sgi:
!
      do iipw=1,ngq(iq)
        cfact=cmplx(1.0d0/sqrt(dabs(epsipw(iipw))),0.0d0,8)
        sgi(:,iipw)=cfact*sgi(:,iipw)
      enddo
      
      deallocate(rwork,work,epsipw)
      end subroutine diagipw
!EOC

