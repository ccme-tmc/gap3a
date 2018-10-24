!BOP
!
! !ROUTINE: kip_readeqp
!
! !INTERFACE:
      subroutine kip_readeqp

! !DESCRIPTION:
! 
! This subroutine reads the qp energies from file
!
! !USES:
!
      use kmeshintp, only: nvbm,ncbm,eks1,eqp1,nkp1,kvecs1,idvk1,klist1,&
     &                     ib0_kip,ib1_kip,nsp_kip,eqptag_kip,  &
     &                     init_kmeshintp,eferqp1
      use task,      only: casename
      implicit none  

! !LOCAL VARIABLES:
      
      integer(4) :: is,ib,ie,io,iu   !(Counter) Runs over bands
      integer(4) :: ik          !(Counter) Runs over k-points
      integer(4) :: itmp,idv
      integer(4) :: ierr
      integer(4):: fid
      integer(4) :: ikvec(1:3)
      
      real(8) :: ehf,elda,rtmp
      character(20):: sname="kip_readeqp"
      character(80):: fname
      logical :: lprt=.false.

!EOP
!BOC
      call linmsg(6,'-',sname) 

      if(trim(eqptag_kip).eq.'') then 
        fname=trim(casename)//".eqpH"
      else 
        fname=trim(casename)//".eqpH_"//trim(eqptag_kip)
      endif

      fid=999
      open(unit=fid,file=fname,action='read',iostat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to open case.eqpH") 

      read(fid,*) ib0_kip,ib1_kip,nkp1,nsp_kip,eferqp1
      if(lprt) write(6,'(a,4i5)') 'nkp1,ib0_kip,ib1_kip',nkp1,ib0_kip, & 
     &                   ib1_kip

      call init_kmeshintp(1)
 
      do is=1,nsp_kip
        if(nsp_kip.eq.2) write(6,*) "  Read QPE for spin ",is 
        if(lprt) write(6,'(2a5,2a12)') 'ik','ie','Enk(KS)','Enk(GW)'

        nvbm(is)=0
        ncbm(is)=ib1_kip

        do ik=1, nkp1
          read(fid,*) itmp,itmp,ikvec(1:3),idv
          klist1(1:3,ik)=ikvec
          kvecs1(1:3,ik)=dble(ikvec)/idv

          io=ib0_kip
          iu=ib1_kip
          do ib=ib0_kip,ib1_kip
            read(fid,2) itmp,eks1(ib,ik,is),eqp1(ib,ik,is)

            if(lprt) write(6,'(2i5,2f12.6)') ik,ib,eks1(ib,ik,is),eqp1(ib,ik,is)
            if(eks1(ib,ik,is).lt.0.0d0)then
              if(ib.gt.io) io=ib
            else
              if(ib.lt.iu) iu=ib
            endif
          enddo

          read(fid,*)
          if(lprt) write(6,*) 
          if(io.gt.nvbm(is)) nvbm(is) = io
          if(iu.lt.ncbm(is)) ncbm(is) = iu
        enddo ! ik
        write(6,'(a,2i5)') "nvbm,ncbm=",nvbm(is),ncbm(is) 

        if(nvbm(is).ge.ncbm(is) ) then 
          write(6,'(a,3i5)') " !! nvbm >= ncbm",nvbm(is),ncbm(is)
          write(6,*) "  nvbm = ncbm -1 is forced!"
          nvbm(is) = ncbm(is)-1
        endif   
      enddo  ! is
      close(fid)
          
    2 format(i4,3f20.15,60x)
    3 format(i4,2f20.15)
      return 
      end subroutine kip_readeqp
!EOC      
