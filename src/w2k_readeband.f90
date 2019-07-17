!BOP
!
! !ROUTINE: readeband
!
! !INTERFACE:
      subroutine readeband
      
! !DESCRIPTION:
!
! This subroutine reads the DFT eigenvalues for every k-point
!
! {\bf WIEN2k interface}
!      
! !USES:
      
      use bands,     only: nspin
      use kmeshintp, only: nbands2,nkp2,kvecs2,eks2,eqp2 
      use struk,     only: nat
      use task,      only: casename,spflag

! !LOCAL VARIABLES:

      implicit none
 
      
      integer(4) :: iat  ! Indexes inequivalent atoms
      integer(4) :: ik,idv   ! Indexes irreducible k-points
      integer(4) :: ib   ! Indexes bands
      integer(4) :: i,isp,itmp
      integer(4) :: io
      integer(4) :: iu
      integer(4) :: nwf,nbk,nbmax
      integer(4) :: fid
      integer(4) :: ierr
      integer(4) :: ikvec(3)
      real(8) :: kvec(3),w
      real(8) :: eryd
      
      character(10) :: kname
      character(13), parameter :: sname = 'readeband'
      character(120) :: fname
      character(8) :: aaa
!EOP
!BOC

!----------------------------------------------------------------------!
!                  Read k-list for band plot                           !
!----------------------------------------------------------------------!

!
!  check how many points to be drawn
!
      fid=999
      open(unit=fid,file=trim(casename)//".klist_band",action='read', &
     &     iostat=ierr) 
      call errmsg(ierr.ne.0,sname,"Fail to open case.klist_band")

      i=1
      do
        read(fid,98) aaa
        if(aaa.eq.'END') exit
        i=i+1
      enddo
      nkp2=i-1

!  read the k-points in case.klist\_band
!
      allocate(kvecs2(3,nkp2))

      rewind(fid)
      do ik=1,nkp2
        read(fid,99) ikvec(1:3),idv
        kvecs2(1:3,ik)=dble(ikvec)/dble(idv)
      enddo
      close(fid) 

!----------------------------------------------------------------------!
!                   Read band energy on the kvecs2                     !
!----------------------------------------------------------------------!


! Get the number of bands 
      fid=999 
      nbands2=10000
      open(unit=fid,file=trim(casename)//".energy_band"//spflag(1), &
     &     action='read',iostat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to open case.energy_band")

      do iat = 1, nat
        read(fid,*)
        read(fid,*)
      enddo
      do ik = 1, nkp2
        read(fid,100) kvec,kname,nwf,nbk,w
        if(nbk.le.nbands2) nbands2 = nbk
        do ib=1,nbk
          read(fid,*)i,eryd
        enddo
      enddo
      write(6,*) ' nkp2,nbands2=',nkp2,nbands2
      close(fid)
!
! Initialize array in kmeshintp
!
      allocate(eks2(nbands2,nkp2,nspin),eqp2(nbands2,nkp2,nspin))

      do isp=1,nspin
        open(unit=fid,file=trim(casename)//".energy_band"//spflag(isp), &
     &     action='read',iostat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to open case.energy_band")

        do iat = 1, nat
          read(fid,*)      
          read(fid,*)      
        enddo
        do ik = 1, nkp2
          read(fid,100) kvec(1:3),kname,nwf,nbk,w
          do ib=1,nbk
            read(fid,*) itmp,eryd
            if(ib.le.nbands2) then 
              eks2(ib,ik,isp)=0.5d+0*eryd     ! Transform the energies to Hartree units
            else 
              cycle
            endif  
          enddo
        enddo ! ik 

        close(fid) 
      enddo  !isp

 98   format(a3)
 99   format(10x,4i5)

  100 format(3e19.12,a10,2i6,f5.1)
  101 format(a,i4)
  
      return
      
      end subroutine readeband
      
!EOC      
