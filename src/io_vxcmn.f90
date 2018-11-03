!BOP
!
! !ROUTINE: io_vxcmn
!
! !INTERFACE:
      subroutine io_vxcmn(io,flag,isym,istat)

! !DESCRIPTION:
!
!This subroutine handls the input/output of the full vxc matrix (vxcmn) 
!
! !USES:

      use selfenergy,  only: vxcmn,vxcnn
      use xcpot,       only: vorbnn,lvorb
      use bands,       only: ibgw,nbgw,nspin
      use kpoints,     only: nkp,nirkp,get_kvec
      use task,        only: casename, fid_outgw
      implicit none

      character,intent(in):: io   ! 'r'/'w' (read/write)
      character,intent(in):: flag ! 'f'/'d' (full/diagonal)
      integer, intent(in):: isym  ! 0/1 - full/irreducible BZ 
      integer,intent(out):: istat ! return status
      

! !LOCAL VARIABLES:

      integer :: fid
      integer:: ibgw_in,nbgw_in,nsp_in,nkp_in  ! the index range and the number of spin channels 
                                               ! that can be read from the external file, 
                                               ! used only when io_read == 1
      integer :: ierr
      integer :: iik, ik, irk  ! the order index of the irred. k-point
      integer :: isp     ! index for spin 
      integer :: ien,iem ! Counter: run over the eigenvalues in the corresponding k-point.
      integer :: itmp
      integer :: nktot
      real(8) :: rtmp
      logical :: lerr
      real(8) :: kvec(3)  

      character(len=80)::sname='io_vxcmn',fn_vxc
     
! !REVISION HISTORY:
! 
! Created  Jan 22, 2010 by H. Jiang
!
!
!EOP
!
!BOC

      if(isym.eq.0) then 
        nktot = nkp
      else 
        nktot = nirkp 
      endif 

      fid=200
      istat = 0
      if(flag.eq.'f'.or.flag.eq.'F') then 
        fn_vxc=trim(casename)//".vxcmn"
        if(io.eq.'r'.or.io.eq.'R') then 
          open(unit=fid,file=fn_vxc,action='read',status='old',iostat=ierr)
          if(ierr.ne.0) then 
            !write(6,*) "WARNING in io_vxc: Fail to open vxcmn file"
            write(fid_outgw,*) "WARNING in io_vxc: Fail to open vxcmn file"
            istat=1
            return 
          endif 

          read(fid,*) nsp_in,nkp_in,ibgw_in,nbgw_in 
          lerr=nsp_in.ne.nspin.or.nkp_in.ne.nktot.or.ibgw_in.gt.ibgw &
     &       .or.nbgw_in.lt.nbgw
          if(lerr) then 
            !write(6,*) "WARNING in io_vxc: inconsistency occurs"
            write(fid_outgw,*) "WARNING in io_vxc: inconsistency occurs"
            istat = 1
            close(fid) 
            return 
          endif 
       
          do isp=1,nspin
            do iik=1,nktot
              read(fid,*) 
              do ien=ibgw_in,nbgw_in
                do iem=ibgw_in,nbgw_in
                  if( (ien.ge.ibgw.and.ien.le.nbgw).and.&
     &                (iem.ge.ibgw.and.iem.le.nbgw)) then  
                   read(fid,2) itmp,itmp,vxcmn(iem,ien,iik,isp)
                 else 
                   read(fid,2) itmp,itmp,rtmp,rtmp
                 endif 
                enddo
              enddo
            enddo
          enddo

        else  !! write 
          open(unit=fid,file=fn_vxc,action='write',iostat=ierr)
          call errmsg(ierr.ne.0,sname,"Fail to open > "//trim(fn_vxc))
          write(fid,*) nspin,nktot,ibgw,nbgw
          do isp=1,nspin
            do iik=1,nktot
              call get_kvec(isym,iik,ik,irk,kvec) 
              write(fid,1) iik,kvec(1:3) 
              do ien=ibgw,nbgw
                do iem=ibgw,nbgw
                  write(fid,2) iem,ien,vxcmn(iem,ien,iik,isp)
                enddo 
              enddo
            enddo
          enddo 
        endif 
!
!     I/O of Diagonal vxc matrix 
!
      else 
        fn_vxc=trim(casename)//".vxcnn"
        if(io.eq.'r'.or.io.eq.'R') then
          open(unit=fid,file=fn_vxc,action='read',status='old',iostat=ierr)
          
          if(ierr.ne.0) then 
            !write(6,*) "WARNING in io_vxc: Fail to open to vxcnn file"
            write(fid_outgw,*) "WARNING in io_vxc: Fail to open to vxcnn file"
            istat=1
            return 
          endif 

          read(fid,*) nsp_in,nkp_in,ibgw_in,nbgw_in
          lerr=nsp_in.ne.nspin.or.nkp_in.ne.nktot.or.ibgw_in.gt.ibgw &
     &     .or.nbgw_in.lt.nbgw

          if(lerr) then 
            !write(6,*) "WARNING in io_vxc: inconsistency occurs"
            write(fid_outgw,*) "WARNING in io_vxc: inconsistency occurs"
            istat = 1
            close(fid)
            return 
          endif 

          do isp=1,nspin
            do iik=1,nktot
              read(fid,*)
              do ien=ibgw_in,nbgw_in
                if((ien.ge.ibgw.and.ien.le.nbgw)) then
                  if(lvorb) then
                    read(fid,12) itmp,vxcnn(ien,iik,isp),vorbnn(ien,iik,isp)
                  else
                    read(fid,13) itmp,vxcnn(ien,iik,isp)
                  endif
                  !write(6,*) itmp, vxcnn(ien,iik,isp)
                else
                  read(fid,*)
                endif
              enddo
            enddo
          enddo

        else
          open(unit=fid,file=fn_vxc,action='write',iostat=ierr)
          call errmsg(ierr.ne.0,sname,"Fail to open-w "//trim(fn_vxc))
          write(fid,*) nspin,nktot,ibgw,nbgw
          do isp=1,nspin
            do iik=1,nktot
              call get_kvec(isym,iik,ik,irk,kvec) 
              write(fid,1) iik, kvec(1:3)
              do ien=ibgw,nbgw
                if(lvorb) then
                  write(fid,12) ien,vxcnn(ien,iik,isp),vorbnn(ien,iik,isp)
                else
                  write(fid,13) ien,vxcnn(ien,iik,isp)
                endif
              enddo
            enddo
          enddo
        endif
      endif 

      close(fid)
          
    1 format(i6,3e16.6)
    2 format(2i6,2e24.12)     

   11 format(3e16.6,2i6)
   12 format(i6,2e24.12)
   13 format(i6,e24.12)

      end subroutine io_vxcmn

