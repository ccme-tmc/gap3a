!BOP
!
! !ROUTINE: kip_readenk
!
! !INTERFACE:
      subroutine kip_readenk 
      
! !DESCRIPTION:
!
! This subroutine  reads the band energy for a general k-mesh and store them in arrays defined
! module kmeshintp
!
! {\bf WIEN2k interface}
!      
! !USES:
      
      use kmeshintp,only: kvecs2,eks2,nbands2,nkp2,wk2,nsp_kip, &
     &                    init_kmeshintp
      use struk,    only: nat 
      use task,     only: casename,spflag

! !LOCAL VARIABLES:
      implicit none

      integer,parameter::nkp_max=10000  !! maximal number of k-points that is possible (?!
      integer,parameter::nb_max=10000  !! maximal number of k-points that is possible (?!
      integer(4) :: iat      ! Indexes inequivalent atoms
      integer(4) :: ik      ! Indexes irreducible k-points
      integer(4) :: i,ib       ! Indexes bandsa
      integer(4) :: isp      ! Index for spin 
      integer(4) :: fid
      integer(4) :: nwf,nbk,nkp,nbmax
      integer    :: ierr 
      real(8)    :: kvec(3),w
      real(8), allocatable:: enk(:)
      
      character(10) :: kname
      character(10), parameter :: sname = 'kip_readenk'
!EOP
!BOC
!

!----------------------------------------------------------------------!
      fid=999
      open(unit=fid,file=trim(casename)//".energy"//spflag(1),&
     &     action='read',iostat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to open case.energy")

      do iat = 1,nat
        read(fid,*)
        read(fid,*)
      enddo

      nkp2=0
      nbands2=nb_max
      nbmax=0
      do ik=1,nkp_max
        read(fid,'(67x,2i6,5x)',iostat=ierr) nwf,nbk
        if(ierr.ne.0) exit
        nkp2=nkp2+1
        if(nbk.gt.nbmax) nbmax=nbk
        if(nbk.lt.nbands2) nbands2=nbk
        do ib=1,nbk
          read(fid,*)
        enddo
      enddo
      close(fid)

      call init_kmeshintp(2)

      allocate(enk(nbmax))

      do isp=1,nsp_kip 
        open(unit=fid,file=trim(casename)//".energy"//spflag(isp),&
     &     action='read',iostat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to open case.energy")

        do iat = 1,nat
          read(fid,*) 
          read(fid,*) 
        enddo
      
        do ik = 1, nkp2
          read(fid,100,err=907) kvec(1:3),kname,nwf,nbk,w
          do ib=1,nbk
            read(fid,*,err=908) i,enk(ib)
          enddo

          kvecs2(1:3,ik)=kvec(1:3)
          wk2(ik)=w
          eks2(1:nbands2,ik,isp)=enk(1:nbands2)
        enddo
        close(fid)
      enddo ! isp

      !! Convert band energy values to the unit of Hartree unit
      eks2=eks2*0.5d+0
      return

  100 format(3e19.12,a10,2i6,f5.1)
  907 call outerr(sname,'error reading energy file, code 907')
  908 call outerr(sname,'error reading energy file, code 908')
      end subroutine kip_readenk
      
!EOC      
