!BOP
!
! !ROUTINE: task_acont
!
! !INTERFACE:
       subroutine task_acont
       use bands,      only: bande,nspin,eferqp,ibgw,nbgw,efermi,eqp,   &
     &                       nbandsgw,nbmax
       use constants,  only: hev
       use freq,       only: omega,nomeg
       use kpoints,    only: kirlist,idvkir,nirkp,nvelgw,nvel,nkp
       use selfenergy, only: sigx,sigc,npar_ac,sacpar,iop_ac,       &
     &                       init_selfenergy,end_selfenergy
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

      integer :: ierr,fidi,fidr,fid
      integer :: isp,irk,ie,ip
      integer :: itmp,ib0,nb0      
      integer :: ib,nbout,nwout,isout,ikout,iwzero
      integer :: iom
       
      real(8) :: enk,delta,omg,rtmp
      real(8) :: lwout,wzero           ! lwout is the length of real frequency for output 
      complex(8) :: sig,dsig,ein,ceqp
      character(20):: blk_acont
      integer,allocatable::ibout(:) 
      complex(8),allocatable :: a(:),poles(:)
      real(8),allocatable :: kvecs(:,:)
      character(3)::tag,gwtag

      integer ::npdos
      real(8) ::emin,emax

! !EXTERNAL ROUTINES: 


      external calceqp
      character(10),external::int2str

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
      blk_acont="acont"
      ierr = loct_parse_isdef(blk_acont)
      if(ierr.eq.1) then
         call loct_parse_block_int(blk_acont,0,0,nbout)
         if(nbout.gt.0) then 
           call loct_parse_block_int  (blk_acont,0,1,isout)
           call loct_parse_block_int  (blk_acont,0,2,ikout)
           call loct_parse_block_float(blk_acont,0,3,lwout)
           call loct_parse_block_int  (blk_acont,0,4,nwout)
           call loct_parse_block_int  (blk_acont,0,5,iwzero) 
           allocate(ibout(nbout))
           do ib=1,nbout
             call loct_parse_block_int(blk_acont,1,ib-1,ibout(ib))
           enddo 
         endif 
      else
         nbout= 0
      endif

      write(6,*) "task_acont: init_selfenergy"
      call init_selfenergy(0)


!     Read the selfenergz matrix elements for imaginary frequencies
!
      write(6,*) "task_acont: io_sxcmn"
      call io_sxcmn('r','d',0,0,1,ierr)

      write(6,*) "task_acont: io_vxcmn"
      call io_vxcmn('r','d',1,ierr)
!
!     Calculate the quasiparticle energies
!
      write(6,*) "task_acont: calceqp"
      call calceqp(0,0,-1)

      call io_eqp('w',0,0,"GW")
      call io_eqp('w',0,1,"GW")

      allocate(kvecs(3,nirkp))
      kvecs=dble(kirlist)/idvkir

      call bandanaly(ibgw,nbgw,nirkp,kvecs,bande(ibgw:nbgw,:,:),&
     &               efermi,nspin,"KS")

      call bandanaly(ibgw,nbgw,nirkp,kvecs,eqp,eferqp,nspin,"GW")

      if(nbout.lt.0) return 

      call calcscgw0
      call io_eqp('w',0,0,"GW0")
      call io_eqp('w',0,1,"GW0")
      call bandanaly(ibgw,nbgw,nirkp,kvecs,eqp,eferqp,nspin,"GW0")
      deallocate(kvecs)

      if(nbout.le.0) goto 99

      write(6,*) "--- output some acont information ----"
      isp=isout
      irk=ikout
      write(6,*) " AC parameters for selected bands"
      write(6,'(a,2i5)') "# isp,irk=",isp,irk 
      do ib=1,nbout
        ie=ibout(ib)
        enk=bande(ie,irk,isp)
        write(6,'(a,i5,f16.6)') "# ie,Enk(eV)=",ie,enk*heV
        do ip = 1, npar_ac
          write(6,'(i5,2f16.6)') ip,sacpar(ip,ie,irk,isp)
        enddo 
      enddo 

      allocate(a(npar_ac),poles(npar_ac/2))
      fidi=999
      fidr=1000

      do ib=1,nbout
        ie=ibout(ib)
        enk=bande(ie,irk,isp)
        a(1:npar_ac)=sacpar(1:npar_ac,ie,irk,isp)

        tag=trim(int2str(ie))
        open(fidi,file=trim(casename)//".outaci-n"//trim(tag))
        open(fidr,file=trim(casename)//".outacr-n"//trim(tag))

! write some information to the beginning of the files

        write(fidi,*) "# check Analytic Continuation on imag axis"
        write(fidr,*) "# check Analytic Continuation on real axis"
        write(fidi,'(a,3i5,f16.4)') "# isp,irk,ie,Enk(eV)=",isp,irk,ie, &
     &             enk*heV
        write(fidi,'(a1,a15,4a16)') '#','Omega','Re\Sigma','Im\Sigma',  &
     &  'Re\Sigma(fitted)', 'Im\Sigma(fitted)'
        write(fidr,'(a1,a15,6a16)') '#','Omega',' Re\Sigma','Im\Sigma', &
     &  "Re\Sigma'","Im\Sigma'"

        do iom=1,nomeg
          ein=cmplx(0.d0,sign(omega(iom),enk))
          call getsac(iop_ac,nomeg,npar_ac,enk,ein,omega,a,sig,dsig)
          if(enk.le.0.d0) sig=conjg(sig)
          write(fidi,'(5f16.6)') omega(iom),sigc(ie,irk,iom,isp),sig
        enddo

        if(iwzero.eq.1) then 
          wzero=enk
        else 
          wzero=0.d0
        endif 

        do iom=0,nwout
          omg=wzero+(iom-nwout/2)*(lwout/nwout)
          ein=cmplx(omg,0.d0) 
          call getsac(iop_ac,nomeg,npar_ac,enk,ein,omega,a,sig,dsig)
!          sig=sig +sigx(ie,irk,isp)   !-vxc(ie,irk,isp)

          write(fidr,'(5f16.4)') (omg-wzero)*hev,sig*hev,dsig
        enddo
        close(fidi)
        close(fidr)
      enddo

      deallocate(a,poles)

 99   return
      
      end subroutine task_acont
!EOC      

