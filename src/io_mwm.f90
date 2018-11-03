      subroutine io_mwm(io,mwm,mst,mend,iomst,iomend,isp,iq,irk,ierr)
!
!     I/O of M*W*M 
!
      use bands,      only: nspin,nbandsgw,nbands_c
      use constants,  only: size_cmplx 
      use task,       only: save_prefix
      implicit none
      character,intent(in) :: io   ! '[rR]/[wW]'
      integer,intent(in):: mst,mend,isp,iq,irk,iomst,iomend
      complex(8),intent(out):: mwm(mst:mend,nbandsgw,iomst:iomend)     
      integer,intent(out):: ierr

      character(len=2):: spflag
      integer :: ie1,ie2,iom
      integer :: irec
      integer :: fid=999
      integer(8):: recl_mwm 
      complex(8) :: xtmp(nbandsgw) 

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

      recl_mwm = size_cmplx*nbandsgw
      fn_mwm=trim(save_prefix)//".mwm"//trim(spflag)//"-q" &
     &       //trim(int2str(iq))//"-k"//trim(int2str(irk))

      if(io.eq.'w'.or.io.eq.'W') then  !! write 
        open(unit=fid,file=fn_mwm,action='write',form='unformatted',&
     &       recl=recl_mwm,access='direct',iostat=ierr)

        if(ierr.ne.0) return 
           
        do iom=iomst,iomend
          do ie2=mst,mend 
            irec = (iom-1)*nbands_c + ie2 

            do ie1=1,nbandsgw 
              xtmp(ie1) = mwm(ie2,ie1,iom) 
            enddo 
            write(fid,rec=irec,iostat=ierr) xtmp 
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

        do iom=iomst,iomend 
          do ie2=mst,mend
            irec = (iom-1)*nbands_c + ie2
             
            read(fid,rec=irec,iostat=ierr) xtmp

            if(ierr.ne.0) then 
              close(fid)
              return 
            endif 

            do ie1=1,nbandsgw
              mwm(ie2,ie1,iom)= xtmp(ie1) 
            enddo 

          enddo 
        enddo 
        close(fid)
      endif
      end subroutine ! internal subroutine readmwm
