!BOP
!
! !ROUTINE: task_emac
!
! !INTERFACE:
      subroutine task_emac
      
! !DESCRIPTION:
!
! This task calculates macroscopic diectric function 
!
! !USES:

      use constants, only: hev
      use freq,      only: nomeg,omega,womeg,nomeg_blk
      use bands,     only: metallic,nomax,numin,nbmaxpol,nspin
      use core,      only: iop_core 
      use kpoints,   only: nkp,nirkp
      use mixbasis,  only: matsiz,init_mixbasis,end_mixbasis
      use barcoul,   only: init_barcoul,end_barcoul,barcevtol,end_barcev,&
                           iop_coul_c
      use dielmat,   only: bandtype,c0_head,emac,init_dielmat,end_dielmat
      use freq,      only: nomeg
      use mommat,    only: iop_mommat
      use task,      only: taskname,casename
      use modmpi 
      use liboct_parser
      
! !LOCAL VARIABLES:

      implicit none
 
      integer :: iop_sym,iop_emac
      integer :: iq, ierr,ip,comm,iom,iom_blk,ist_blk,iend_blk,   &
     &           nom_blk,iom0,iom1,nom_p
      real(8) :: wpl
 
      character(20):: sname="task_emac"
      character(len=10):: blk_emac='emac'
      character(120) :: fname 

      complex(8),allocatable :: emac_a(:,:),emac_nlf(:)
      integer::rstype
      integer :: iomcnt(0:nproc_max),iomdsp(0:nproc_max)
!
!EOP
!BOC            

      call linmsg(6,"-","TASK:"//trim(taskname))

!
! Read task-specific input parameters from *.ingw 
!
!   iop_emac 
!     0 -- RPA without local field effect
!     1 -- RPA with local field effect
!     2 -- calculate plasmon frequency only 
!   iop_sym: whether to use the symmetry when summing over k
!     0 -- not use symmetry 
!     1 -- use the symmetry  
!   iop_mommat -- control how to treat moment matrix 
!     0 => calculate within the code 
!     1 => read from the external file 
!
!   bandtype 
!    'KS' -- use KS orbital energies 
!    'QP' -- use QP orbital energies 
        
      iq=1
      bandtype='KS'
      iop_emac = 1
      iop_mommat = 0
      iop_sym = 0
      ierr = loct_parse_isdef(blk_emac)
      if(ierr.eq.1) then 
        call loct_parse_block_int(blk_emac,0,0,iop_emac)
        call loct_parse_block_int(blk_emac,0,1,iop_sym)
      endif 

      if(iop_mommat.eq.1) then
        iop_core = 2
      endif 

      allocate(emac_a(2,nomeg),emac_nlf(nomeg),stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate emac_a")

      if(bandtype .eq. 'QP') then 
        write(6,*) "task_emac: GW-RPA using interpolated QP energies "
        write(6,*) "   !!! This feature not fully tested "
        return 
        call kip_readeqp
        call kip_readenk          !* Read the eigenenergies from the Wien2k "case.energy" file
        call kip_qpeintp
      endif 

      if(iop_emac.eq.2.and.metallic) then 
        call calcplasmon(wpl,iop_sym) 
        write(6,'(a,f12.4)')"Plasmon frequency (eV):",wpl*hev
        return 
      endif 

      !! Bare Coulumb matrix 
      if(iop_emac.eq.1) then 
        call init_mixbasis(iq) 
        call init_barcoul(iq)
      
        write(6,*) "task_emac: coul_barc"
        call coul_barc(iq)
        call coul_setev(iq,barcevtol,iop_coul_c)
      endif 
!
! to make calculations of large number of freq points possible, all freq. points are divided into blocks
!   nom_blk  - the real number of freq points within the current block 
!
      do iom_blk=1, nomeg,nomeg_blk 
        nom_blk=min(iom_blk+nomeg_blk-1,nomeg)-iom_blk+1

        if(iop_emac.eq.0) then  !! without local field effect
          call calcemac_nlf(emac_nlf(iom_blk:iom_blk+nom_blk-1),&
                iom_blk,nom_blk)

        else  !! with the local field effect 

          !! set the parallelization over frequency
          call mpi_set_range(nproc_col,myrank_col,nom_blk,iom_blk, &
     &        iom0,iom1,iomcnt,iomdsp)
          write(6,*) "iom0,iom1=",iom0,iom1 
          nom_p = iom1 - iom0 + 1
          call init_dielmat(iq,iom0,iom1)

          call calceps(iq,iom0,iom1,iop_sym,-1,1,.false.)

          if(nproc_col.gt.1) then
            write(6,*) sname//":collect emac data for different freq"
#ifdef MPI
            call MPI_Type_contiguous(2,MPI_DOUBLE_COMPLEX,rstype,ierr)
            call MPI_Type_commit(rstype,ierr)
            call MPI_GatherV(emac,nom_p,rstype,emac_a(:,iom0:iom1),&
     &                iomcnt,iomdsp,rstype,0,mycomm_row,ierr)
            call MPI_Type_free(rstype,ierr)
#endif
          else
            emac_a(:,iom0:iom1) = emac(:,iom0:iom1)
          endif
          call end_dielmat(iq)
        endif 

      enddo !! iom_blk

      if(iop_emac.eq.1) then 
        call end_barcev 
        call end_barcoul(iq)
        call end_mixbasis
      endif

      if(myrank_row.eq.0.and.myrank_col.eq.0) then 
        if(metallic) write(6,'(a,f12.4)')"Plasmon frequency (eV):",sqrt(c0_head)*hev
        write(6,*) sname//" write out emac data "

        if(iop_emac.eq.1) then 
          fname=trim(casename)//".emac_"//bandtype
          open(unit=999,file=fname,action='write')
          write(999,10)
          write(999,11)
          write(6,10)
          write(6,11)

          do iom=1,nomeg
            write(999,12) omega(iom)*hev,emac_a(1:2,iom),dimag(-1.d0/emac_a(2,iom))
            write(6  ,12) omega(iom)*hev,emac_a(1:2,iom),dimag(-1.d0/emac_a(2,iom))
          enddo ! iom

        else
          fname=trim(casename)//".emac_nlf_"//bandtype
          open(unit=999,file=fname,action='write')
          write(999,20)
          write(999,21)
          write(6,20)
          write(6,21)

          do iom=1,nomeg
            write(999,22) omega(iom)*hev,emac_nlf(iom),dimag(-1.d0/emac_nlf(iom))
            write(6  ,22) omega(iom)*hev,emac_nlf(iom),dimag(-1.d0/emac_nlf(iom))
          enddo ! iom
        endif
        close(999)
      endif 
! 
!     deallocate the momentum matrices         
!
      deallocate(emac_a,emac_nlf)
      return

 10   format("# emac with and without local field effects")
 11   format("# \omega(eV)",14x," \eps_M     ",10x,          &
     &                      10x," \eps_M(NLF)",10x,          &
     &                      "    Im(-1/eps_M)")
 12   format(f16.6,5e16.6)

 20   format("# emac without local field effects")
 21   format("# \omega(eV)",14x,"   \eps_M   ",10x,          &
     &                      "    Im(-1/eps_M)")
 22   format(f16.6,3e16.6)

      end subroutine task_emac
!EOC
