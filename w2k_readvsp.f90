!BOP
!
! !ROUTINE: w2k_readvsp
!
! !INTERFACE:
      subroutine w2k_readvsp
      
! !USES:
      use bands,   only: nspin
      use struk,   only: nat,nrpt
      use radwf,   only: nrad,vr
      use task,    only: casename
      
      implicit none
      
! !LOCAL VARIABLES:

      integer(4) :: fid    ! file id for *.vsp file 
      integer(4) :: isp    ! index for spin 
      integer(4) :: idummy ! a very dummy variable
      integer(4) :: ir     ! index for radial mesh points
      integer(4) :: iat    ! index for inequivalent atom
      integer(4) :: npt    ! number of radial mesh points 
      integer:: ierr
      
      character(len=10), parameter :: sname = 'w2k_readvsp'
      
      
! !DESCRIPTION:
!
! {\bf WIEN2k interface:}
!
! This subroutine reads the spherical potential at each muffin-tin sphere
! from file case.vsp
!
! !REVISION HISTORY:
!
! Created Feb. 5th. 2004
! Modified Dec 08, 2008 by JH
!
!EOP
!
!BOC
      fid=999
      open(fid,file=trim(casename)//".vsp",action='read',iostat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to open vsp file")

      allocate(vr(nrad,nat,nspin))
      do isp=1,nspin 
        read(fid,5070,err=901)
        do iat = 1, nat
          npt=nrpt(iat)
          read(fid,5010,err=901)
          read(fid,5020,err=901) idummy
          read(fid,5060,err=901)
          read(fid,5040,err=901) (vr(ir,iat,isp),ir=1,npt)
          read(fid,5060,err=901)
          read(fid,5050,err=901)
          vr(1:npt,iat,isp) = vr(1:npt,iat,isp)/2.0d+0
        enddo 
      enddo  ! isp 
      close(999)
      return 

 901  call errmsg(.true.,sname,"ERROR in reading case.vsp")
!
! formats for wien2k.03
!
 5010 format (3x)
 5020 format (15x,i3,/,/)
 5040 format (3x,4e19.12)
 5050 format (/,/,/)
 5060 format (/)
 5070 format (/,/)
!
! formats for wien2k.02
!
! 5010 FORMAT (3X)
! 5020 FORMAT (16X,I2,/,/)
! 5040 FORMAT (3X,4E19.12)
! 5050 FORMAT (/,/,/)
! 5060 FORMAT (/)
! 5070 FORMAT (/,/)

      end subroutine w2k_readvsp
!EOC      
