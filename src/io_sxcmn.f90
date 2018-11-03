!BOP
!
! !ROUTINE: io_sxcmn
!
! !INTERFACE: 
      subroutine io_sxcmn(io,diag,iq,isxc,isym,ierr)

! !DESCRIPTION:
!
!  This subroutine reads (io=='r' )or writes (io == 'w' ) full selfenergy matrix 
!  represented by KS wavefunctions 
!
! !USES:

      use bands,      only: ibgw,nbgw,nspin,nbandsgw
      use constants,  only: hev
      use freq,       only: nomeg
      use kpoints,    only: nkp,nirkp,get_kvec
      use selfenergy, only: sigc,sigc_q,sigx,sigx_q,sigm,info_sxc
      use task,       only: casename,savdir,fid_outgw
      
      implicit none
     
      character,intent(in) :: io   ! r/R for reading and w/W for writing    
      character,intent(in) :: diag ! 'd'/'f' diagonal/full sxc matrix elements  
      integer,  intent(in) :: iq
      integer,intent(in)   :: isxc    ! the type of selfenergy to be read
                                  !   0 -- Read both exchange and correlation selfenergy 
                                  !   1 -- exchange-only
                                  !   2 -- static Coulomb hole screened exchange
                                  !   3 -- SEX-only 
      integer,intent(in)   :: isym ! 0/1 -- full/irreducible BZ 
      integer,intent(out)  :: ierr
      
      
! !LOCAL VARIABLES:
      integer :: fid 
      integer :: ib, jb
      integer :: ie, ie1, ie2
      integer :: iq_in    !
      integer :: isp    ! index for spin 
      integer :: iik,ik,irk,nktot
      integer :: iom, nom
      integer :: nsp_in,nkp_in,ibgw_in,nbgw_in,nom_in 
      logical :: ldiag

      complex(8),pointer:: sx_p(:,:,:),sc_p(:,:,:,:)

      real(8) :: kvec(3)
      real(8) :: rtmp
      logical :: lerr
      character(80)::fname,prefix
      character(3) :: diag_tag
      character(10) :: sxc_tag 

      character(10), external:: int2str

      character(20):: sname="io_sxcmn"

!EOP
!BOC
!

!
!     Set file name 
!
      if(diag.eq.'d') then 
        if(iq.gt.0) then 
          sx_p => sigx_q
          sc_p => sigc_q
        else
          sx_p => sigx
          sc_p => sigc 
        endif 
      endif 

      fid=999

      if(isxc.eq.1) then 
        sxc_tag = 'sx'
        nom=0
      elseif(isxc.eq.0) then 
        sxc_tag = "sxc"
        nom = nomeg 
      elseif(isxc.eq.2) then 
        sxc_tag = "sxc_chsx"
        nom = 1
      elseif(isxc.eq.3) then 
        sxc_tag = "sxc_sx"
        nom = 1
      endif 

      if(isym.eq.0) then
        nktot=nkp
      else
        nktot=nirkp
      endif 

      if(diag.eq.'d'.or.diag.eq.'D') then 
        ldiag=.true.
        diag_tag='_nn'
      else 
        ldiag=.false.
        diag_tag='_mn'
      endif 
      fname=trim(casename)//"."//trim(sxc_tag)//trim(diag_tag)
      if(iq.ne.0) then
        fname=trim(savdir)//trim(fname)//"-q"//trim(int2str(iq))
      endif 


      if(io.eq.'r'.or.io.eq.'R') then 

        open(unit=fid,file=fname,action='read',iostat=ierr)
        if(ierr.ne.0) return 

        read(fid,100,iostat=ierr) iq_in 
        if(ierr.ne.0 .or. iq_in.ne.iq ) then 
          !write(6,*) "WARNING: inconsistent iq!" 
          write(fid_outgw,*) "WARNING: inconsistent iq!" 
          close(fid) 
          return 
        endif 

        !! consistency check  
        read(fid,101) nsp_in,nkp_in,ibgw_in,nbgw_in,nom_in 
        lerr=nsp_in.ne.nspin.or.nkp_in.ne.nktot.or.ibgw_in.gt.ibgw  &
     &     .or.nbgw_in.lt.nbgw.or.nom_in.ne.nom

        if(lerr) then 
          ierr = 1 
          !write(6,*) "WARNING: inconsistent data in "//trim(fname)
          write(fid_outgw,*) "WARNING: inconsistent data in "//trim(fname)
          close(fid) 
          return 
        endif 

        do isp=1,nspin
          do iik=1,nktot
            read(fid,*)    !! the line with k-vector information 
            if(ldiag) then !! read diagonal elements
              do ie=ibgw_in,nbgw_in
                if( ie.lt.ibgw.or.ie.gt.nbgw ) then
                  read(fid,*)
                  cycle
                endif
                select case(isxc)
                case(0)   !! full GW
                  read(fid,302) sx_p(ie,iik,isp),sc_p(1:nom,ie,iik,isp)
                case(1)
                  read(fid,302) sx_p(ie,iik,isp)
                case(2)
                  read(fid,302) sx_p(ie,iik,isp),sc_p(1:3,ie,iik,isp)
                case(3)
                  read(fid,302) sx_p(ie,iik,isp),sc_p(1,ie,iik,isp)
                endselect
              enddo !! ie

            else  !! read whole matrix
              do ie2=ibgw_in,nbgw_in
                do ie1=ibgw_in,nbgw_in 
                  if(ie1.lt.ibgw.or.ie1.gt.nbgw.or. &
     &               ie2.lt.ibgw.or.ie2.gt.nbgw) then 
                    read(fid,*) 
                  else 
                    read(fid,310,advance='no') ib,jb
                    do iom=0,nom 
                      read(fid,320,advance='no')sigm(ie1,ie2,iom,iik,isp)
                    enddo 
                    read(fid,*) 
                  endif 
                enddo !! ie1
              enddo !! ie2
            endif 
          enddo !! irk
        enddo !! isp
        close(fid)
      endif !! io.eq.'r'

      if(io.eq.'w'.or.io.eq.'W') then
        open(unit=fid,file=fname,action='write',iostat=ierr)
        if(ierr.ne.0) then 
          !write(6,*) "WARNING: fail to open "//trim(fname)
          write(fid_outgw,*) "WARNING: fail to open "//trim(fname)
          return 
        endif  

        write(fid,200) iq,trim(info_sxc(isxc))
        write(fid,201) nspin,nktot,ibgw,nbgw,nom

        do isp=1,nspin
          do iik=1,nktot
            call get_kvec(isym,iik,ik,irk,kvec) 
            write(fid,301) ik,irk,kvec  
            if(ldiag) then !! write diagonal elements
              do ie=ibgw,nbgw
                selectcase(isxc)
                case(0)
                  write(fid,302) sx_p(ie,iik,isp),sc_p(:,ie,iik,isp)
                case(1)
                  write(fid,302) sx_p(ie,iik,isp)
                case(2)
                  write(fid,302) sx_p(ie,iik,isp),sc_p(1:3,ie,iik,isp)
                case(3)
                  write(fid,302) sx_p(ie,iik,isp),sc_p(1,ie,iik,isp)
                endselect
              enddo !! ie
            else   !! write whole matrix elements
              do ie2=ibgw,nbgw
                do ie1=ibgw,nbgw
                  write(fid,310,advance='no') ie1, ie2
                  do iom=0,nom
                    write(fid,320,advance='no')sigm(ie1,ie2,iom,iik,isp)
                  enddo
                  write(fid,330) ' '
                enddo !! ie1
              enddo !! ie2
            endif 
          enddo !! irk
        enddo !! isp
        close(fid)
      endif !! io.eq.'w'
      return

 100  format(5x,i5)
 101  format(5i5)
 200  format("# iq=",i5,"  Matrix elements of ",a)
 201  format(5i5)


 301  format(2i5,3F12.6,4X,"# ik,irk,kvec(1:3)")  
 302  format(6g18.10) 
 310  format(2I6)
 320  format(2g24.16)
 330  format(A)
      end subroutine io_sxcmn
!EOC          
            
