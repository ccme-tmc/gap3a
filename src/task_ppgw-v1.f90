!BOP
!
! !ROUTINE: task_ppgw
!
! !INTERFACE:
       subroutine task_ppgw
       use bands,      only: bande,nspin,eferqp,ibgw,nbgw,efermi,eqp,   &
     &                       nbandsgw,nbmax,fspin 
       use bzinteg,    only: kwt_ibz
       use bzint,      only: bzint_smear,esmear 
       
       use constants,  only: hev
       use freq,       only: omega,nomeg,freq_eta 
  
       use kpoints,    only: kirvecs,nirkp
       use selfenergy, only: init_selfenergy,end_selfenergy,selfx,sigc
       use xcpot,      only: init_xcpot,vxcnn
       use task,       only: casename
       use liboct_parser
       implicit none
       
! !DESCRIPTION:
!
! This subroutine perform the analytic continuation of the selfenergy to
! real frequencies and calculates the quasi-particle energies. The
! selfenergy and exchange correlation potential are read from file, thus, a
! previous run of the GW cycle is needed.        
!

      integer :: ierr,fid
      integer :: isp,ie,ip,irk
      integer :: n_min,n_max
      integer :: iom
       
      real(8) :: enk,delta,omg,rtmp 
      character(20)::blk_ppgw="ppgw",sname="task_gw"

      integer :: iop_pp,irk_in  
      integer :: nom
      integer :: fout 
      real(8) :: om_min,om_max,exnk,om_step,efer
      real(8) :: sumdos_qp,sumdos_ks,sumspf_ef,sumspf
      real(8),allocatable :: oms(:),spf(:), spf_n(:),dos_qp(:),dos_ks(:)

! !EXTERNAL ROUTINES: 


      external calceqp
      character(10),external::int2str

      fout = 999

!

!      
! !REVISION HISTORY:
!
! Created 16.09.2005 by RGA
!
!EOP
!BOC

!
! read some input that control the output that tests the acont 
!
! iop_pp -- indicate which kind of post-processing to do 
!   0 -- redo the analytic continuation with the new options  
!   1 -- calculate the total spectral function for a given range of bands
!   2 -- calculate the spectral function for a partical k-point 
! 

      ierr = loct_parse_isdef(blk_ppgw)
      if(ierr.eq.1) then
        call loct_parse_block_int(blk_ppgw,0,0,iop_pp)

        if(iop_pp.eq.1) then 
          call loct_parse_block_int  (blk_ppgw,0,1,n_min)
          call loct_parse_block_int  (blk_ppgw,0,2,n_max)
          call loct_parse_block_int  (blk_ppgw,0,3,nom)
          call loct_parse_block_float(blk_ppgw,0,4,om_min)
          call loct_parse_block_float(blk_ppgw,0,5,om_max)
          call loct_parse_block_float(blk_ppgw,0,6,freq_eta)
          esmear = freq_eta 
        elseif(iop_pp.eq.2) then 
          call loct_parse_block_int  (blk_ppgw,0,1,irk_in)
          call loct_parse_block_int  (blk_ppgw,0,2,n_min)
          call loct_parse_block_int  (blk_ppgw,0,3,n_max)
          call loct_parse_block_int  (blk_ppgw,0,4,nom)
          call loct_parse_block_float(blk_ppgw,0,5,om_min)
          call loct_parse_block_float(blk_ppgw,0,6,om_max)
        endif 
      else
        iop_pp=0
      endif

      write(6,*) "task_ppgw: init_selfenergy"
      call init_selfenergy(0)

!     Read the selfenergz matrix elements for imaginary frequencies
!
      write(6,*) "task_ppgw: io_sxcmn"
      call io_sxcmn('r','d',0,0,1,ierr)

      write(6,*) "task_ppgw: io_vxcmn"
      call init_xcpot
      call io_vxcmn('r','d',1,ierr)


!
!     Calculate the quasiparticle energies
!
      write(6,*) "task_ppgw: calceqp"
      call calceqp(0,0,-1)
      write(6,*) ":E_Fermi(QP)=",eferqp 

      call io_eqp('w',0,0,"GW")
      call io_eqp('w',0,1,"GW")

      call bandanaly(ibgw,nbgw,nirkp,kirvecs,bande(ibgw:nbgw,:,:),&
     &               efermi,nspin,"KS")

      call bandanaly(ibgw,nbgw,nirkp,kirvecs,eqp,eferqp,nspin,"GW")

!
!     Calculate the spectural function
!
      if(iop_pp.eq.1) then 

        !! set the real frequency grid 
        allocate(oms(nom),spf(nom),spf_n(nom),dos_qp(nom),dos_ks(nom),&
     &   stat=ierr) 
        call errmsg(ierr.ne.0,sname,"fail to allocate oms,...") 
          
           
        om_step = (om_max - om_min)/(nom-1) 
        do iom=1,nom
          oms(iom) = om_min + (iom-1.0)*om_step 
        enddo 

        open(fout,file=trim(casename)//"-Ank.dat",action='write')
        spf = 0.d0 
        do isp=1,nspin 
          do irk=1,nirkp 
            do ie = n_min, n_max 
              exnk = bande(ie,irk,isp) - vxcnn(ie,irk,isp) + selfx(ie,irk,isp)

              write(fout,100) ie,irk, kwt_ibz(irk),bande(ie,irk,isp)

              call calc_spectfun(spf_n,oms+eferqp,nom,exnk,sigc(:,ie,irk,isp),&
     &         omega,nomeg,fout) 
              spf = spf + spf_n*kwt_ibz(irk)*fspin
            enddo
          enddo
        enddo 
        close(fout) 

        ! calculate DOS directly from QP and KS energies 
        dos_qp = 0.d0 
        dos_ks = 0.d0 
        do isp=1,nspin
          do irk=1,nirkp
            do ie = n_min, n_max
              do iom=1,nom 
                dos_qp(iom) = dos_qp(iom)+fspin*kwt_ibz(irk) &
     &            *bzint_smear(2,oms(iom)-eqp(ie,irk,isp)+eferqp)
                dos_ks(iom) = dos_ks(iom)+fspin*kwt_ibz(irk) &
     &            *bzint_smear(2,oms(iom)-bande(ie,irk,isp))
              enddo 
            enddo 
          enddo
        enddo 

!       Determine the Fermi energy in terms of the spectral function 

        open(fout,file=trim(casename)//"-specfun.dat",action='write') 

        sumspf_ef = 0.d0
        sumspf = 0.d0
        sumdos_qp = 0.d0 
        sumdos_ks = 0.d0 
        do iom=1,nom 
          sumspf = sumspf + om_step * spf(iom)
          if(oms(iom).le.0.d0) then 
            sumspf_ef = sumspf_ef + om_step * spf(iom)
            sumdos_ks = sumdos_ks + om_step * dos_ks(iom) 
            sumdos_qp = sumdos_qp + om_step * dos_qp(iom)
          endif

          write(fout,120) oms(iom)*HeV,spf(iom)/HeV,dos_qp(iom)/HeV, &
     &     dos_ks(iom)/HeV 
        enddo 
        write(6,110) "Integrating DOS_QP up to E_Fermi",sumdos_qp 
        write(6,110) "Integrating DOS_KS up to E_Fermi",sumdos_ks
        write(6,110) "Integrating spectr. func. up to E_Fermi",sumspf_ef
        write(6,110) "Integrating spectr. func. to infinity ",sumspf
        close(fout)  
        deallocate(oms,spf,spf_n) 
      endif 

      return
 100  format(/,"# Spectral function for n=",i4," irk=",i4," k-weight=",&
     & f10.6," Enk=",f12.6) 
 110  format(a40,F10.4) 
 120  format(4F12.6) 
      
      end subroutine task_ppgw
!EOC      

