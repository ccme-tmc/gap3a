      subroutine io_eps(io,iq,iomfirst,iomlast,ierr)
      use dielmat,    only: head,epsw1,epsw2,eps
      use constants,  only: size_cmplx,size_int
      use mixbasis,   only: matsiz
      use task,       only: savdir,casename
      use modmpi,     only: myrank 
      implicit none 
      character,intent(in):: io !! 'r/w/m'
      integer,intent(in) :: iq
      integer,intent(in) :: iomfirst,iomlast    !! the range of frequency points
      integer,intent(out):: ierr

      integer   :: fid 
      integer   :: iom,iq_in,itmp
      integer(8):: recl_eps            !! record length for eps files
      character(10):: sname="io_eps"
      character(120):: msg,msg_r,msg_w
      character(len=120)::fn        !! file name 
      character(len=10),external:: int2str

      fid = 900
      fn=trim(savdir)//trim(casename)//".eps-q"//trim(int2str(iq))

!      if(io.eq.'W') then 
!        fn=trim(fn)//'-'//trim(int2str(myrank))
!      endif 

      if(io.eq.'w' .or. io.eq.'W') then 
        write(6,*) 'Write the eps matrix to ',fn
      elseif(io.eq.'r' .or. io.eq.'R') then 
        write(6,*) 'Read the eps matrix from ',fn
      endif 

      msg="Fail to open "//trim(fn)//" for direct access: "//io  
      msg_r="Fail to read from "//trim(fn)  
      msg_w="Fail to wrtie to  "//trim(fn)  

      if(io.eq.'m') then !! just read the matsiz 
        open(fid,file=fn,action='read',form='unformatted',&
     &       access='direct',recl=size_int,iostat=ierr)
        if(ierr.ne.0) then 
          call wrnmsg(.true.,sname,msg) 
          return 
        endif  

        read(fid,rec=1) matsiz
        close(fid) 
        write(6,*) "matsiz from "//trim(fn),matsiz 
        return 
      endif 

      recl_eps = size_int+size_cmplx*matsiz*matsiz 
      if(iq.eq.1) then 
        recl_eps = recl_eps + size_cmplx*(matsiz*2+1) 
      endif 

      if(io.eq.'w' .or. io.eq.'W') then !! write binary  
        open(fid,file=fn,action='write',form='unformatted',&
     &       access='direct',recl=recl_eps,iostat=ierr)
        if(ierr.ne.0) then 
          call wrnmsg(.true.,sname,msg) 
          return 
        endif 
        do iom=iomfirst,iomlast
          if(iq.eq.1) then 
            write(fid,rec=iom,iostat=ierr) matsiz,head(iom),epsw1(:,iom),&
     &                         epsw2(:,iom),eps(:,:,iom)
          else 
            write(fid,rec=iom,iostat=ierr) matsiz,eps(:,:,iom)
          endif 
          if(ierr.ne.0) then 
            call wrnmsg(.true.,sname,msg_w) 
            write(6,*) "iom, ierr=",iom, ierr 
            close(fid) 
            return 
          endif 
        enddo 
        close(fid)

      elseif(io.eq.'r') then !! read binary 
        open(fid,file=fn,action='read',form='unformatted',&
     &       access='direct',recl=recl_eps,iostat=ierr)

        if(ierr.ne.0) then 
          call wrnmsg(.true.,sname,msg) 
          return 
        endif 
        do iom=iomfirst,iomlast
          if(iq.eq.1) then 
            read(fid,rec=iom,iostat=ierr) itmp,head(iom),epsw1(:,iom),  &
     &       epsw2(:,iom),eps(:,:,iom)
          else 
            read(fid,rec=iom,iostat=ierr) itmp,eps(:,:,iom)
          endif 
          if(ierr.ne.0) then
            call wrnmsg(.true.,sname,msg_r)
            write(6,*) "ierr=",ierr 
            close(fid)
            return
          endif
        enddo
        close(fid)

      endif
 100  format(4g24.12)  


      end subroutine
