!BOP
! !ROUTINE: io_initfiles
!
! !INTERFACE: 
      subroutine io_initfiles
!
!
! !DESCRIPTION:
!
! This subroutine set up some global I/O files including scratch files 
!
! !USES:

      use bands,      only: nspin
      use task,       only: casename,fid_outdbg,fid_outmb,iop_scratch,&
     &                   fid_outmom,fid_outkpt,scrdir,spflag,scrfn,&
     &                   savdir,save_prefix,fid_outqp
      use eigenvec,   only: vfunit,vfname
      use modmpi

! !LOCAL VARIABLES:
!
      implicit none

      integer(4) :: i,nlen
      integer(4) :: ierr
      real(8) :: r
      character(len=9), parameter :: sname = 'io_initfiles'
      character(len=120)::pbsjobid
      character(len=10),external::Int2Str
      integer, external:: system 
!EOP
!BOC

!
! In the case of parallel calculation, attach a process rank to some output files 
!
      write(6,*) " setup scratch (direct access) files"

      call getenv("SCRATCH",scrdir)
      nlen=len_trim(scrdir)
      if(nlen.eq.0) then 
        write(6,*) "WARNING: SCRATCH is not defined!"
        write(6,*) "  -- using the current directory as scratch"
        scrdir='scr/'
      endif 

!     check whether scrdir exisits 
      open(999,file=trim(scrdir)//'.scr'//trim(procflag),iostat=ierr)
      if(ierr.ne.0) then 
        write(6,*) "WARNING in io_initfiles: scratch unaccessible! ",& 
     &            trim(scrdir)
        write(6,*) "  -- using the current directory as scratch"
        scrdir='scr/'
        if(myrank.eq.0) ierr = system("mkdir -p "//trim(scrdir))
      endif 
      close(999,status='delete')

      nlen=len_trim(scrdir)
      if(scrdir(nlen:nlen).ne.'/') then
        scrdir=trim(scrdir)//'/'
      endif

      write(6,*) " Scratch directory: ", trim(scrdir) 

      !! To avoid possible confilction between different calculations, the names of scratch files are 
      !! appended by a random number  
      call init_random_seed
      call random_number(r)
      scrfn=trim(scrdir)//trim(casename)//trim(int2str(nint(1000.0*r)))

!
!     check save directory
!
      if(myrank.eq.0) then 
        open(999,file=trim(savdir)//'/.sav'//trim(procflag),iostat=ierr)
        if(ierr.ne.0) then
          write(6,*) "WARNING in io_initfiles: savdir does not exist! ",&
     &            trim(savdir)
          savdir='tmp/'
          ierr = system("mkdir -p "//trim(savdir))
        else 
          close(999,status='delete')
        endif 
      endif 

#ifdef MPI
      call mpi_barrier(MPI_COMM_WORLD,ierr) 
#endif 

      nlen=len_trim(savdir)
      if(savdir(nlen:nlen).ne.'/') then
        savdir=trim(savdir)//'/'
      endif
      write(6,*) " Save directory: ", trim(savdir)

      save_prefix = trim(savdir)//trim(casename)

! vector file 

      if(iop_scratch.gt.0) then 
        vfunit=101
        if(iop_scratch.eq.1) then 
          vfname=trim(scrfn)//".vectord"//trim(procflag)
        else
          vfname=trim(save_prefix)//".vectord"
        endif 
      endif 
        
      write(6,*) "  Direct access vector file: ",trim(vfname)

! setup some other I/O files
      open(fid_outmb, file=trim(save_prefix)//".outmb", action='write')
      open(fid_outdbg,file=trim(save_prefix)//".outdbg",action='write')
      open(fid_outmom,file=trim(save_prefix)//".outmom",action='write')
      open(fid_outkpt,file=trim(save_prefix)//".outkpt",action='write')
      open(fid_outqp, file=trim(save_prefix)//".outqp", action='write')

      return
      end subroutine io_initfiles
!EOC
