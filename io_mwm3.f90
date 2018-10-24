      subroutine io_mwm3(io,mwm,m1,m2,iom1,iom2,isp,iq,irk,ierr)
!
!     I/O of M*W*M 
!
      use bands,      only: nspin,nbandsgw,nbands_c
      use constants,  only: size_cmplx 
      use task,       only: save_prefix
      implicit none
      character,intent(in) :: io   ! '[rR]/[wW]'
      integer,intent(in):: m1,m2,isp,iq,irk,iom1,iom2
      complex(8),intent(inout)::mwm(nbandsgw,nbandsgw,m1:m2,iom1:iom2)  
      integer,intent(out):: ierr

      character(len=2):: spflag
      integer :: m,iom
      integer :: irec
      integer :: fid=999
      integer(8):: recl_mwm 
      complex(8) :: xtmp(nbandsgw,nbandsgw) 

      character(len=120)::fn_mwm
      character(len=10),external::int2str
      character(len=10):: sname="io_mwm"

      !! set the mwm file name
      if(nspin.eq.1) then
        spflag=''
      else
        if(isp.eq.1) then
          spflag='up'
        else
          spflag='dn'
        endif
      endif

      recl_mwm = size_cmplx*nbandsgw*nbandsgw 
      fn_mwm=trim(save_prefix)//".mwm"//trim(spflag)//"-q" &
     &       //trim(int2str(iq))//"-k"//trim(int2str(irk))

      if(io.eq.'w'.or.io.eq.'W') then  !! write 
        open(unit=fid,file=fn_mwm,action='write',form='unformatted',&
     &       recl=recl_mwm,access='direct',iostat=ierr)

        if(ierr.ne.0) return 
           
        do iom=iom1,iom2
          do m=m1,m2 
            irec = (iom-1)*nbands_c + m
            write(fid,rec=irec,iostat=ierr) mwm(:,:,m,iom)
            if(ierr.ne.0) then 
              close(fid)
              return 
            endif 
          enddo
        enddo 
        close(fid)
      else                !! read 
        open(unit=fid,file=fn_mwm,action='read',form='unformatted',&
     &       recl=recl_mwm,access='direct',iostat=ierr)

        if(ierr.ne.0) return 

        do iom=iom1,iom2 
          do m=m1,m2
            irec = (iom-1)*nbands_c + m
            read(fid,rec=irec,iostat=ierr) mwm(:,:,m,iom)

            if(ierr.ne.0) then 
              close(fid)
              return 
            endif 

          enddo 
        enddo 
        close(fid)
      endif
      end subroutine ! internal subroutine readmwm
