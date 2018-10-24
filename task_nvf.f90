!BOP
!
! !ROUTINE: task_nvf
!
! !INTERFACE:
      subroutine task_nvf
      
! !DESCRIPTION:
!
! This task subroutine 
!   (1) add extrapolated GW quasi-particle quasi-particle corrections to LDA eigen-energies 
!   (2) Write new eigen-energies to *.energy_gw and *.vector_gw, which is going to be used in 
!      WIEN2k programs   
! !USES:
     
      use eigenvec,    only: lcmplx
      use kmeshintp,   only: iop_kip,kvecs2,nbands2,nkp2,nsp_kip,eqp2,  &
     &                       eferqp2,eferqp1
      use lapwlo,      only: lmax,lomax,elapw,elo,nloat,l_newlo,nlmax
      use struk,       only: nat
      use task,        only: casename,spflag
      use liboct_parser
      
! !LOCAL VARIABLES:

      implicit none

      integer:: isp            ! index for spin
      integer:: iat            ! index for non-equivalent atoms 
      integer:: ib,itmp             ! index for energy band  
      integer:: ig             ! index for G-vector 
      integer:: ik         ! index for k-points 
      integer:: l,k
      integer:: ilo
      integer:: nwf,nbk,nwfmax
      integer:: ierr

      
      

      real(8):: kvec(3)
      real(8):: weight,eigval,rtmp

      integer:: fin,fout            ! id for output file  
      character(80)::fnm_nef,fnm_ef      ! name for new and old energy file 
      character(80)::fnm_nvf,fnm_vf      ! name for new and old vector file 
      character(900) :: str       ! string used to read and write the head of energy file a

      character(10) :: kname=""
      character(3)  :: ipgr 
      character(10) :: sname="task_nvf"
      character(10) :: blk_nvf="nvf"
      character(500):: sline


      real(8),    allocatable::  rzzk(:)
      complex(8), allocatable::  zzk(:)
      integer, allocatable :: kzz(:,:)
      
! !REVISION HISTORY:
!
! Created 21.06.2007 by Hong Jiang
!
!EOP
!BOC            

!
!  Obtain QP energies on the fine k-mesh by interpolating the QP energies on 
!  the sparse k-mesh 
!
      call kip_readeqp 
      call kip_readenk          !* Read the eigenenergies from the Wien2k "case.energy" file
      call kip_qpeintp
!
! extrapolate Quasi-particle energy correction 
!

      do isp=1,nsp_kip 
!
!  set names of the new vector and energy file 
!
        fnm_nef=trim(casename)//".energy"//trim(spflag(isp))//"_nvf"
        fnm_nvf=trim(casename)//".vector"//trim(spflag(isp))//"_nvf"
        fnm_ef=trim(casename)//".energy"//trim(spflag(isp))
        fnm_vf=trim(casename)//".vector"//trim(spflag(isp))

        fin=200
        fout=300
!
!  write energy file
!
        open(unit=fout,file=trim(fnm_nef),action='write',iostat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fnm_nef))

        open(unit=fin,file=trim(fnm_ef),action='read',iostat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fnm_ef))

        write(6,*) "task_nvf: write new energy file"
        do iat=1,nat
          read(fin,*)  
          read(fin,*)
          if (l_newlo) then
            write(fout,'(1000f12.5)') elapw(0:lmax,iat,isp)
            write(fout,'(1000f12.5)') ((elo(l,ilo,iat,isp),l=0,lomax),ilo=0,nloat-1)
          else 
            write(fout,'(100f9.5)') elapw(0:lmax,iat,isp)
            write(fout,'(100f9.5)') elo(0:lomax,0:nloat-1,iat,isp)
          endif 
        enddo 

        nwfmax=0
        do ik = 1, nkp2
          read(fin,100)  kvec,kname,nwf,nbk,weight
          if(nwf.gt.nwfmax) nwfmax=nwf
          write(fout,100) kvec,kname,nwf,nbands2,weight
          do ib=1,nbk
            if(ib.le.nbands2) write(fout,*) ib,eqp2(ib,ik,isp)*2.d0
            read(fin,*) 
          enddo
        enddo
        close(fin)
        close(fout)

        allocate(kzz(1:3,nwfmax))
        if(lcmplx)  then
          allocate(zzk(1:nwfmax))
          zzk = 0.d0
        else 
          allocate(rzzk(nwfmax))
          rzzk = 0.d0
        endif
!
! Write vector file 
!
        write(6,*) "task_nvf: write new vector file(s)"
        open(unit=fout,file=trim(fnm_nvf),action='write',iostat=ierr,&
     &       form='unformatted')
        call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fnm_nvf))

        open(unit=fin,file=trim(fnm_vf),action='read',iostat=ierr,&
     &       form='unformatted')
        call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fnm_vf))
        do iat=1,nat
          read(fin) rtmp
          read(fin) rtmp
        enddo

        do iat=1,nat
          write(fout)  elapw(0:nlmax-1,iat,isp)
          write(fout)  elo(0:lomax,0:nloat-1,iat,isp)
        enddo 

        do ik=1,nkp2 
!
! Read from old-vector file
!
          read(fin,iostat=ierr) kvec(1:3),kname,nwf,nbk,weight
          call errmsg(ierr.ne.0,sname,"Fail to read vector file - 1 ")

          read(fin,iostat=ierr) (kzz(1:3,ig),ig=1,nwf)
          call errmsg(ierr.ne.0,sname,"Fail to read vector file - 2 ")

! Write to new vector file
          write(fout) kvec(1:3),kname,nwf,nbands2,weight,ipgr
          write(fout) (kzz(1:3,ig),ig=1,nwf)

          do ib=1,nbk

! Read from old-vector file
            read(fin,iostat=ierr)  itmp,eigval
            call errmsg(ierr.ne.0,sname,"Fail to read vector file - 3 ")
            if(lcmplx) then
              read(fin,iostat=ierr) (zzk(ig),ig=1,nwf)
            else
              read(fin,iostat=ierr) (rzzk(ig),ig=1,nwf)
            endif
            call errmsg(ierr.ne.0,sname,"Fail to read vector file - 4 ")

! Write to new vector file
            if(ib.le.nbands2) then 
              write(fout) ib,eqp2(ib,ik,isp)*2.d0 
              if(lcmplx) then
                write(fout) (zzk(ig),ig=1,nwf)
              else
                write(fout) (rzzk(ig),ig=1,nwf)  
              endif 
            endif 

          enddo   
        enddo  
        close(fout) 
        close(fin)
        deallocate(kzz)
        if(lcmplx) then 
          deallocate(zzk) 
        else 
          deallocate(rzzk)
        endif 
      enddo ! isp
      write(6,'(a,f10.4,a)')  ":E_FERMI_QP1 ",eferqp1*2.d0," Ry."
      write(6,'(a,f10.4,a)')  ":E_FERMI_QP2 ",eferqp2*2.d0," Ry." 
      return

  100 format(3e19.12,a10,2i6,f5.1)
  101 format(a900)
  102 format(a)
  99  format("ERROR in task_nvf -- ", a) 

      end subroutine task_nvf
!EOC
