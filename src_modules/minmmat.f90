!-----------------------------------------------------------------------
!BOI
!
! !MODULE: minmmat
!  This module defines variables related to Minm, the matrix elements between
!  the product of two eigenvectors and a mixed basis function
! 
      module minmmat
!EOE
!BOP
      use bands,   only: nbmax,nspin
      use core,    only: ncg
      use kpoints, only: nkp
      use mixbasis,only: matsiz 
      use task,    only: save_prefix 
      implicit none 
      
! !PUBLIC VARIABLES:
      integer :: mblksiz=10    ! to reduce the memory usage, the m-index is divided 
                               ! into blocks, each with a size of mblksiz   
      integer :: fid_minm 
      integer :: recl_minm     ! record length for minm direct access file  
      integer :: irec_last

      integer,allocatable ::idx_minm(:,:,:,:)  ! this array is used to store the information about 
                                               ! the position of minm for given n,m, 
                                               ! zero indicating the corresponding minm
                                               ! is not calculated yet 


      character(len=120) :: fn_idx_minm,fn_minm 
     
!EOP
      
!BOC 
! !LOCAL VARIABLES
      character(10),external:: int2str
      integer,private::ierr

      contains

      subroutine init_minm(iq,iop)
!
!     Initialize idx_minm array
!       iop < 0 -> Minm will be always calculated on the fly 
!       iop ==0 -> Minm is calculated and saved 
!       iop > 0 -> Minm is read from the file 
!
      implicit none 
      integer,intent(in) :: iq, iop
      integer:: isp,ik

      if(iop .lt. 0) return 

      fid_minm = 1000+iq
      fn_idx_minm=trim(save_prefix)//".idx_minm-q"//trim(int2str(iq))
      fn_minm=trim(save_prefix)//".minm-q"//trim(int2str(iq))
      recl_minm = 16*matsiz  
      
      allocate(idx_minm(nbmax+ncg,nbmax+ncg,nkp,nspin))

      !! check whether idx_minm exists
      idx_minm = 0
      irec_last = 0
      if(iop.gt.0) then 
        open(999,file=fn_idx_minm,action='read',iostat=ierr)
        if(ierr.ne.0) then 
          write(6,*) "WARNING in init_minm: idx_minm does not exist"
          write(6,*) " -- calc. Minm from the scratch "
        else
          do isp=1,nspin 
            do ik=1,nkp
              read(999,*) idx_minm(:,:,ik,isp)
            enddo
          enddo
          irec_last = maxval(idx_minm)
        endif 
        close(999)
      endif 

      open(fid_minm,file=fn_minm,action='readwrite',recl=recl_minm,  &
     &      form='unformatted',access='direct',iostat=ierr)
      call errmsg(ierr.ne.0,"init_minm","Fail to open minm d.a. file "  &
     &         //trim(fn_minm))
      end subroutine 
!
!     Clean up Minm 
!
      subroutine end_minm
      implicit none

      if(allocated(idx_minm)) then
        call save_idx_minm
        deallocate(idx_minm) 
        close(fid_minm)
      endif 
      end subroutine 
!
!     Save idx_minm array
!
      subroutine save_idx_minm
        implicit none
        integer:: ik,iq,isp

        open(999,file=fn_idx_minm,action='write',iostat=ierr)
        if(ierr.ne.0) then
          write(6,*)"WARNING in save_idx_minm: fail to open the idx &
     &file -- idx_minm is not saved !!!"
        else  
          do isp=1,nspin
            do ik=1,nkp
              write(999,*) idx_minm(:,:,ik,isp)
            enddo
          enddo
          close(999)
        endif
      end subroutine

!
!     Check the status of Minm in the given range in the file 
!
      subroutine check_minm(nmflag,nst,nend,mst,mend,ik,isp,stat) 
        implicit none
        character,intent(in):: nmflag*2
        integer,intent(in):: nst,nend,mst,mend,ik, isp
        integer,intent(out) :: stat !!  0 -- all absent, 1 -- all present, 2 -- partially present 

        integer:: irec
        integer:: n, m,n0,m0
        character(20):: sname="minmmat%io_minm"
        logical:: l_allnew, l_allold

        l_allnew = .true.
        l_allold = .false. 
        call set_n0m0(nmflag,n0,m0)
        do m = mst,mend
          do n= nst, nend
            irec = idx_minm(m0+m,n0+n,ik,isp)  !! record index for writing
            if(irec.eq.0) then
              l_allold = .false.
            else 
              l_allnew = .false.
            endif 
          enddo
        enddo
        if(l_allnew) then 
          stat = 0
        elseif(l_allold) then 
          stat = 1
        else 
          stat = 2
        endif 
      end subroutine 

      subroutine set_n0m0(nmflag,n0,m0) 
        implicit none 
        character,intent(in):: nmflag*2
        integer, intent(out):: n0, m0
        if(nmflag.eq.'nm') then
          n0=ncg
          m0=ncg
        elseif(nmflag.eq.'nc') then
          n0=ncg
          m0=0
        elseif(nmflag.eq.'cm') then
          n0=0
          m0=ncg
        else
          n0=0
          m0=0
        endif
      endsubroutine 

      subroutine io_minm(io,nmflag,minm,nst,nend,mst,mend,ik,isp)
        implicit none 
        character,intent(in):: io,nmflag*2
        complex(8),intent(inout):: minm(matsiz,mst:mend,nst:nend)
        integer,intent(in):: nst,nend,mst,mend,ik, isp

        integer:: irec
        integer:: n, m,n0,m0
        logical:: lprt=.true.
        character(20):: sname="minmmat%io_minm"
        call set_n0m0(nmflag,n0,m0) 
        do n= nst, nend
          do m = mst,mend
            if (io.eq.'w'.or. io.eq.'W') then 
              irec = idx_minm(m0+m,n0+n,ik,isp)  !! record index for writing
              if(irec.eq.0) then 
                irec_last = irec_last + 1
                irec = irec_last 
                idx_minm(m0+m,n0+n,ik,isp) = irec 
              endif 

              write(fid_minm,rec = irec,iostat=ierr) minm(:,m,n)
              call errmsg(ierr.ne.0,sname,"Fail to write minm")  
            else 
              irec = idx_minm(m0+m,n0+n,ik,isp) 
              if(irec.eq.0) then 
                write(6,100) "ERROR in io_minm: mi"//nmflag//" for n=", &
     &                     n," m=",m," not available"
                stop "ERROR in io_minm"
              endif 
              read(fid_minm,rec=irec,iostat=ierr) minm(:,m,n)
              call errmsg(ierr.ne.0,sname,"Fail to read minm")
            endif 
          enddo
        enddo
  100   format(a,i4,a,i4,a)
      end subroutine 

      end module 
!EOC
