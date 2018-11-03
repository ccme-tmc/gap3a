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
       use xcpot,      only: init_xcpot
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

      integer :: ierr
      integer :: n_min,n_max
       
      real(8) :: enk,delta,omg,rtmp 
      character(20)::blk_ppgw="ppgw",sname="task_gw"

      integer :: iop_pp,irk_in  
      integer :: nom
      real(8) :: om_min,om_max,efer

! !EXTERNAL ROUTINES: 

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
        call calc_spectfun(0,n_min,n_max,nom,om_min,om_max)
      endif 
      
      end subroutine task_ppgw
!EOC      

