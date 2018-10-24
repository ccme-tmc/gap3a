!BOP
!
! !ROUTINE: readingw
!
! !INTERFACE:
      subroutine readingw()
      
!
! !DESCRIPTION:
!
! This subroutine reads the main input file for the gw calculation
! \verb"case.ingw"

! !USES:
 
      use bands,     only: nspin,fspin,emingw,emaxgw,ibgw,nbgw, &
     &                      iop_metallic,spinmom,band_scissor,       &
     &                      emaxpol,eminpol,emaxsc,eminsc,efermi
      use barcoul,   only: barcevtol,rcut_coul,stctol,&
     &                      iop_coul_x,iop_coul_c,iop_coulvm

      use bzint,     only: iop_bzint,iop_bzintq,n_gauq,eta_freq,     &
     &                      ztol_sorteq,tol_taylor,tol_unphys_weight,   &
     &                      esmear,nomg_ht,omgmax_ht,ztol_vol,iop_bcor 

      use core,      only: core_ortho,iop_core
      use constants, only: hev,pi
      use crpa,      only: iop_pln_phase,wf_centers_new,nlorb_wf,  &
     &                     tol_shift_centers 
      use dielmat,   only: q0_eps,wt_excl,noc_excl,nun_excl,ioc_excl,  &
     &                      iun_excl,iop_mask_eps,occ_win,unocc_win,&
     &                      iop_drude,omega_plasma,iop_epsw,eta_head 
      use eigenvec,  only: lsymvector,lcmplx
      use fouri,     only: rmax
      use freq,      only: omegmax,nomeg,nomeg_blk,iop_fgrid,omegmin,   &
     &                     iop_freq
      use kmeshintp, only: iop_kip,eqptag_kip 
      use kpoints,   only: nkdivs,k0shift,nvel
      use lapwlo,    only: lomax
      use minmmat,   only: mblksiz
      use mixbasis,  only: kmr,pwm,lmbmax,lblmax,wftol,nspin_mb,nlo_mb, &
     &                     mb_ludot,mb_emin, mb_emax 
      use modmpi,    only: myrank, nproc_row, nproc_col 
      use selfenergy,only: npol_ac,npar_ac,iop_es,iop_ac,beta_sc,&
     &                     iop_esgw0,iop_gw0
      use struk,     only: nat
      use xcpot,     only: lvorb,natorb,iatorb,nlorb,lorb,vorb,&
     &                      dmorb,iop_vxc,lhybrid
      use task,      only: taskname,casename,gwinp,gwout,&
     &                lrestart,spflag,iop_dftpkg,iop_scratch,savdir, &
     &                progname,l_usesavedminm,l_save_dielmat
      use liboct_parser
! !LOCAL VARIABLES:      
      implicit none
      
!
! !LOCAL VARIABLES:
!
      integer:: i
      integer:: ierr 
      integer:: il    ! index for l
      integer:: ia    ! index for inequivalent atoms
      integer:: iorb
      integer:: jatom
      integer:: l     ! the azimutal quatum number
      integer:: lind  ! l value for which mixopt is read.
      integer:: irow

      integer iarg,iargc

      logical :: opl,opr
     
      character(10):: sname='readingw' 
      character(4) :: chr,readop,csymvector
      character(5) :: bzcon, interp
      character(6) :: fdep
      character(3) :: corflag

! !name of the input blocks, used for easier modifications 
      character(20):: blk_bcoul = "BareCoul"
      character(20):: blk_bzcon = "BZConv"
      character(20):: blk_bzint = "bzint"
      character(20):: blk_epsmask = "eps_mask"
      character(20):: blk_fgrid = "FreqGrid"
      character(20):: blk_fouri =  "FourIntp"
      character(20):: blk_kip   = "kMeshIntp"
      character(20):: blk_kmesh = "kmesh"
      character(20):: blk_mixbas= "MixBasis"
      character(20):: blk_mpi   = "mpi"
      character(20):: blk_para  = "ParaOpt"
      character(20):: blk_q0eps = "q0eps"
      character(20):: blk_selfe = "SelfEnergy"
      character(20):: blk_vorb  = "Vorb" 
      character(20):: blk_wann  = "Wannier" 
      character(80):: default_casename
      character(80):: stdout="gw.log"
      character(len=10),external::Int2Str

!EOP
!BOC
      open(6,file=stdout,action='write')  
      open(999,file=gwinp,iostat=ierr,action='read') 
      call errmsg(ierr.ne.0,sname,"Fail to open main input file")
      close(999)

      write(6,'(a)') "Main input from "//trim(gwinp)   

!
! initialize the liboct_parser and open the input file 
!
      ierr = loct_parse_init(stdout)
      if(ierr .ne. 0) then
         write(6,'(a)') '*** Fatal Error (description follows)'
         write(6,'(2x,a)') 'Error initializing liboct'
         write(6,'(2x,a)') 'write permissions in this directory?'
         call outerr(sname,"")
      end if

      ierr = loct_parse_input(gwinp)
      if(ierr .ne. 0) then
        write(6,'(a)') '*** Fatal Error ***'
        write(6,'(2x,a)') 'Error initializing liboct'
        write(6,'(2x,a)') 'Can not open input file!'
        write(6,'(2x,a)') 'Please provide an input file with name inp'
        call outerr(sname,"")
      endif
      call boxmsg(6,'*',"GW Program : "//trim(progname))

      call linmsg(6,'-',sname) 
      call loct_parse_string('CaseName','gw',casename)
      write(6,100) "CaseName: "//casename

      call loct_parse_int("UseScratch",2,iop_scratch)
      if(iop_scratch.eq.0)then
        write(6,*) " -- Do not use scratch space"
      else if (iop_scratch.eq.1) then
        write(6,*) " -- Use different scratch files for different proc."
      else
        write(6,*) " -- Use same scratch files for different proc."
      endif

      !! check whether the case name matches the files in the current directory
      call check_file(trim(casename)//".struct",ierr) 
      call errmsg(ierr.ne.0,sname,trim(casename)//".struct not found!")

      call loct_parse_string('SavDir','tmp',savdir)
      write(6,100) "SavDir: "//trim(savdir)

      call io_initfiles

      gwout=trim(casename)//".outgw"
      if(myrank.eq.0) then
        close(6) 
        open(6,file=gwout,position='append')
      else 
        close(6) 
        stdout = trim(savdir)//'/'//trim(stdout)//'-p'//trim(int2str(myrank))
        open(6,file= stdout,action='write') 
      endif  
!
!     set the DFT package used to generate Kohn-Sham inputs
!
      call loct_parse_int("dftpkg",0,iop_dftpkg)

      !! Read the task to calculate
      call loct_parse_string('Task','gw',taskname)
      write(6,100) "Task: "//taskname

        
      call loct_parse_logical('Restart',.false.,lrestart) 
      if(lrestart) then 
        write(6,100) "Restart= .TRUE. (Restart from previous run"
      else  
        write(6,100) "Restart= .FALSE.(Start a new calculation)"
      endif 

      call loct_parse_logical("UseSavedMinm",.true.,l_UseSavedMinm)
 
      call loct_parse_logical('SaveDielmat',.true.,l_save_dielmat)

!
! Set spin polarization option 
!
      call loct_parse_int("nspin",1,nspin)
      if(nspin.ne.1 .and. nspin.ne.2) then 
        write(6,100) "Wrong options for nspin"
        write(6,100) "  nspin=1 --- spin-unpolarized caculations"
        write(6,100) "       =2 --- spin-polarized caculations"
        write(6,100) "  -- set to the default value (nspin=1)"
        nspin=1
      endif 
      fspin=2/nspin 
      write(6,101) "nspin:",nspin
      if(nspin.eq.2) then
        spflag(1)='up'
        spflag(2)='dn'
      else
        spflag=''
      endif


      call loct_parse_int("iop_core",0,iop_core)

      call loct_parse_int("LOmax",3,lomax) 


!
!     whether using vector files taking symmetry into account
! 
      call loct_parse_logical("SymVector",.false.,lsymvector)
      if(lsymvector) then 
        write(6,100) " Use symmetrized eigenvector file"
      else 
        write(6,100) " Use non-symmetrized eigenvector file"
      endif 

!
! Set the number of unoccupied bands for which GW band correction is to be calculated
! 
      call loct_parse_float("emaxgw", 1.0d4,emaxgw) 
      call loct_parse_float("emingw",-1.0d4,emingw) 

      !! Convert to Hartree unit
      emaxgw = emaxgw*0.5d0
      emingw = emingw*0.5d0    

      !! these two I/O are used to be compatable with old gw calculations
      call loct_parse_int("ibgw",-1,ibgw)
      call loct_parse_int("nbgw",-1,nbgw)
      call loct_parse_float("beta_sc",1.0,beta_sc)
!
! Block size used for Minm matrix operations
!
      write(6,100) "Options related to Minm:"
      call loct_parse_int("Minm_mblksiz",48,mblksiz)
      write(6,101) "block size for m-index(mblksiz):",mblksiz
!      write(6,*) "  whether scratch minm (iop_minm):",iop_minm

!
! Read the energy cut-off for polarization matrix and correlation selfenergies 
!  when emaxpol/emaxsc is negative, all unoccupied bands are used  
!
      call loct_parse_float("eminpol",-1.0D10,eminpol) 
      call loct_parse_float("emaxpol",-1.0D10,emaxpol) 
      call loct_parse_float("emaxsc", -1.0D10,emaxsc) 
      call loct_parse_float("eminsc", -1.0D10,eminsc) 
      eminpol = eminpol * 0.5D0
      emaxpol = emaxpol * 0.5d0  !! convert to unit of Ha. 
      emaxsc  = emaxsc  * 0.5d0  !! convert to unit of Ha. 
      eminsc  = eminsc  * 0.5d0  !! convert to unit of Ha. 

      call loct_parse_float("nvel",-1.d0,nvel) 

      call loct_parse_logical("Core_ortho",.false.,core_ortho)

!
! Is the input WIEN2K vector complex or not 
!
      call loct_parse_logical("ComplexVector",.false.,lcmplx) 
      write(6,100) "Complex or real KS vectors?"
      if(lcmplx) then 
        write(6,100) " -- Complex Vector"
      else 
        write(6,100) " -- Real Vector"
      endif 

      
!    iop_vxc : control how to treat vxc 
!          0 -- calculate from vxc data read from the wien2k case.r2v file 
!          1 -- directly read from an external file
      call loct_parse_int("iop_vxc",0,iop_vxc)
      call loct_parse_logical("lhybrid",.false.,lhybrid) 


      call loct_parse_float("nvel",0.d0,nvel)

      !! read the option regarding the additional phase factor for pln
      call loct_parse_int('iop_pln_phase',0,iop_pln_phase) 
      call loct_parse_float('tol_shift_centers',0.01,tol_shift_centers)

      ierr = loct_parse_isdef(blk_wann)
      if(ierr.eq.1) then 
        call loct_parse_block_int(blk_wann,0,0,nlorb_wf)
        allocate(wf_centers_new(3,nlorb_wf)) 
        do iorb=1,nlorb_wf 
          do i=1,3
            call loct_parse_block_float(blk_wann,iorb,i-1,&
     &             wf_centers_new(i,iorb))
          enddo
        enddo 
      else
        nlorb_wf = 0 
      endif 
!
! Set options for BZ integration  
!   iop_bzint --- control standard BZ integration 
!             = 0 -- the old way to calculate weights 
!             = 1 -- use the smearing
!   iop_bzintq -- control how to obtain q-dependent BZ integration weights
!       0 -- analytic (old)
!       1 -- q-BZint by simple summation with smearing 
!      -1 --  numerical integration for weights
!      -2 -- q-BZint by Hilbert transform (HT) approach 
!
      ierr = loct_parse_isdef(blk_bzint)
      if(ierr.eq.1) then
        call loct_parse_block_int(blk_bzint,0,0,iop_bzint)
        call loct_parse_block_int(blk_bzint,0,1,iop_bzintq)
        if(iop_bzintq.eq.0) then 
          write(6,*) "  q-BZint by analytic weights"

        elseif(iop_bzintq.eq.-1) then  !! numerical integration
          write(6,*) "  q-BZint by numerical weights"

        elseif(iop_bzintq.eq.-2) then 

          write(6,*) "  q-BZint by Hilbert transform (HT) approach"
          call loct_parse_block_int(blk_bzint,1,0,nomg_ht)
          call loct_parse_block_float(blk_bzint,1,1,omgmax_ht)

        elseif(iop_bzintq.eq.1) then 

          write(6,*) "  q-BZint by simple summation with smearing"

        else
          write(6,*) "ERROR: unsupported option for iop_bzintq!"
          stop
        endif 
      else
        iop_bzint  = 0 
        iop_bzintq = 0
      endif

      call loct_parse_float("eta_head",    0.01d0, eta_head    )
      call loct_parse_float("ztol_sorteq", 0.01d0, ztol_sorteq )
      call loct_parse_float("tol_taylor",  10.0d0, tol_taylor  )
      call loct_parse_float("esmear",      0.005,  esmear      )
      call loct_parse_float("eta_freq",    0.01d0, eta_freq    )
      call loct_parse_int  ("n_gauq",      8,      n_gauq      )
      call loct_parse_float("ztol_vol",    1.0d-10, ztol_vol   ) 
      call loct_parse_int(  "iop_bcor",    0,      iop_bcor    )

      call loct_parse_float("EFermi",1.0d4,efermi) 
      if(efermi.lt.1.e2) then 
        EFermi = EFermi*0.50
        write(6,200) "Read Fermi energy from the gw input:",efermi 
      endif   
!
!     Parametars to do "insulating" calculations from a metallic
!     starting point 
! 
      call loct_parse_int("iop_metallic",0,iop_metallic)
      call loct_parse_float("SpinMoment",0.d0,spinmom)
      call loct_parse_float("BandScissor",0.d0,band_scissor)

!  BZ convolution parameters are defined in the block 
!    %BZConv
!      bzcon |	fdep 
!    % 
      ierr = loct_parse_isdef(blk_bzcon)
      if(ierr.eq.1) then 
         call loct_parse_block_string(blk_bzcon,0,0,bzcon)
         call loct_parse_block_string(blk_bzcon,0,1,fdep) 
      else 
         bzcon='tetra'
         fdep='imfreq'
      endif

      select case (fdep)
        case ('nofreq')
          iop_freq = 0 
          if((taskname.eq.'gwsc'))then
            write(6,100)"WARNING: fdep=nofreq inconsistent with 'gwsc'"
            write(6,100)'   -- change to default value:imfreq'
            fdep='imfreq'
            iop_freq = 3
          endif  
        case ('refreq')
          iop_freq = 2
        case ('imfreq','IMFREQ')
          iop_freq = 3 
        case default  
          write(6,100) "WARNING: unsupported option for fdep!"
          write(6,100) "--Taking default value: imfreq"
          fdep='imfreq'
          iop_freq = 3 
      end select  
      write(6,100) "fdep="//trim(fdep) 

!
! Parameters for Fourier interpolation
!
     
      ierr=loct_parse_isdef(blk_fouri)
      if(ierr.eq.1) then 
         call loct_parse_block_float(blk_fouri,0,0,rmax)
      else 
         rmax=40.0
      endif

      ierr = loct_parse_isdef(blk_kip)
      if(ierr.ne.1) then
        iop_kip=0
        eqptag_kip=''
      else
        call loct_parse_block_int(blk_kip,0,0,iop_kip)
        call loct_parse_block_string(blk_kip,0,1,eqptag_kip)
      endif

      call loct_parse_int  ("iop_drude",1,iop_drude) 
      call loct_parse_float("Omega_Plasma",-1.d0,omega_plasma)
      omega_plasma = omega_plasma/hev

      call loct_parse_int  ("iop_epsw",0,iop_epsw) 

!
! this block reads parameters that define "eps_mask" 
! when calculating the dielect matrix, the transition from bands defined ioc_excl(:) 
! to those given in iun_excl(:) will be weighted by 
! wt_excl, if wt_excl=0.0, then they are completely exclued. 
!  the first parameter, iop_mask_eps, determines how the mask is defined 
!   iop_mask_eps == 1: define the mask in terms of the band index
!   iop_mask_eps == 2: define the mask in terms of the energy range 
!
      ierr= loct_parse_isdef(blk_epsmask)
      iop_mask_eps = 0
      if(ierr.eq.1) then 
        call loct_parse_block_int  (blk_epsmask,0,0,iop_mask_eps)
        call loct_parse_block_float(blk_epsmask,0,1,wt_excl) 
        if (iop_mask_eps.eq.1) then  ! define the mask in terms of the band index 
          call loct_parse_block_int(blk_epsmask,1,0,noc_excl) 
          if(noc_excl.gt.0) then 
            allocate(ioc_excl(noc_excl))
            call loct_parse_block_int(blk_epsmask,1,1,ioc_excl(1)) 
            do i=2,noc_excl
              ioc_excl(i)=ioc_excl(i-1)+1
            enddo 
          endif 

          call loct_parse_block_int(blk_epsmask,2,0,nun_excl)
          if(nun_excl.gt.0) then
            allocate(iun_excl(nun_excl))
            call loct_parse_block_int(blk_epsmask,2,1,iun_excl(1))
            do i=2,nun_excl
              iun_excl(i)=iun_excl(i-1)+1
            enddo 
          endif

        elseif(iop_mask_eps.eq.-1) then ! define the mask in terms of energy (in units of eV))
          call loct_parse_block_float(blk_epsmask,1,0,occ_win(1))
          call loct_parse_block_float(blk_epsmask,1,1,occ_win(2))
          call loct_parse_block_float(blk_epsmask,2,0,unocc_win(1))
          call loct_parse_block_float(blk_epsmask,2,1,unocc_win(2))
          occ_win   = occ_win/HeV             ! converted to Ha unit 
          unocc_win = unocc_win/HeV
          write(6,*) "The energy window to mask the eps matrix:"
          write(6,'(A,2F12.6)') "Occ. States:  ",occ_win(:)
          write(6,'(A,2F12.6)') "Unocc. States:",unocc_win(:)
        endif 
      endif 
!
! block::q0eps -- the direction q --> 0  
!
      ierr=loct_parse_isdef(blk_q0eps)
      if(ierr.eq.1) then 
        call loct_parse_block_float(blk_q0eps,0,0,q0_eps(1))
        call loct_parse_block_float(blk_q0eps,0,1,q0_eps(2))
        call loct_parse_block_float(blk_q0eps,0,2,q0_eps(3))
      else
        q0_eps=1.d0/sqrt(3.d0)
      endif  

     

!
!     Read the data for the frequecy grid
!
! Frequency grid data is given in the block 
!  %FreqGrid
!    iop_fgrid nomeg omegmax omegmin nomeg_blk
!  % 
!   iop_fgrid =  1 eqdist
!                2 Gauss-Laguerre 
!                3 double Gauss-Legende
!                4 fermion Matsubara frequency

      iop_fgrid=3
      nomeg=16
      omegmax=0.42
      omegmin=0.0    ! minimal Matsubara frequencies at 290.1 K 
      nomeg_blk=0
      ierr = loct_parse_isdef(blk_fgrid)
      if(ierr.eq.1) then 
         call loct_parse_block_int   (blk_fgrid,0,0, iop_fgrid  )
         call loct_parse_block_int   (blk_fgrid,0,1, nomeg  )
         call loct_parse_block_float (blk_fgrid,0,2, omegmax)
         call loct_parse_block_float (blk_fgrid,0,3, omegmin)
         call loct_parse_block_int   (blk_fgrid,0,4, nomeg_blk)
      endif
!
      if(nomeg_blk.gt.nomeg .or.nomeg_blk.le.0 ) then 
        write(6,100) "WARNING: nomeg_blk is out of range!"
        write(6,100) " -- use the default(==nomeg) instead"
        nomeg_blk=nomeg
      endif 

      if(trim(taskname).eq.'chsx'.or.trim(taskname).eq.'ex') then 
        nomeg=1
        omegmax=0.d0
      endif 
      write(6,103) iop_fgrid,nomeg,omegmax

!
!     MPI 
!
#ifdef MPI 
      ierr = loct_parse_isdef(blk_mpi) 
      if (ierr.eq.1) then 
        call loct_parse_block_int ( blk_mpi,0,0,nproc_col) 
        call loct_parse_block_int ( blk_mpi,0,1,nproc_row) 
      else
        nproc_col = 0
        nproc_row = 1
      endif
#endif  

!
!     k-mesh
!
      ierr = loct_parse_isdef(blk_kmesh)
      if(ierr.eq.1) then
        call loct_parse_block_int   (blk_kmesh,0,0,nkdivs(1))
        call loct_parse_block_int   (blk_kmesh,0,1,nkdivs(2))
        call loct_parse_block_int   (blk_kmesh,0,2,nkdivs(3))
        call loct_parse_block_int   (blk_kmesh,1,0,k0shift(1))
        call loct_parse_block_int   (blk_kmesh,1,1,k0shift(2))
        call loct_parse_block_int   (blk_kmesh,1,2,k0shift(3))
        call loct_parse_block_int   (blk_kmesh,1,3,k0shift(4))
      else
        nkdivs=1
        k0shift(1:3)=0
        k0shift(4)=1
      endif
      write(6,'(a,7i5)') "k-mesh: (nkx,nky,nkz,deltak)=",&
     &  nkdivs,k0shift

!
! Read options for self-energy 
!  %SelfEnergy
!    <maxexp>  |  <iop_es=0/1 > 
! %

      npol_ac=2
      iop_es=0
      iop_ac=1
      ierr=loct_parse_isdef(blk_selfe)
      if(ierr.eq.1) then 
        call loct_parse_block_int(blk_selfe,0,0,npol_ac,2)
        call loct_parse_block_int(blk_selfe,0,1,iop_es,2) 
        call loct_parse_block_int(blk_selfe,0,2,iop_ac,1) 
      endif 

      if(npol_ac.eq.0) then
        write(6,100) "The input npol_ac == 0: use default npol_ac setup"
        if(iop_ac.eq.1) then
          npol_ac=2
        else
          npol_ac=nomeg/2
        endif
      endif
      npar_ac=2*npol_ac
      if(nomeg.lt.npar_ac) then
        write(6,100)'WARNING: not enough freq for analytic continuation'
        write(6,*)'  - npar_ac .gt. nomeg'
        write(6,*)'  - npar_ac,nomeg=',npar_ac,nomeg
        write(6,*)'  - reset npar_ac =nomeg'
        npar_ac=nomeg
        npol_ac=npar_ac/2
      endif 

      write(6,*) '- Nr. of poles used in analytic continuation:',npol_ac
      write(6,*) '- Option for calculating selfenergy(iop_es): ',iop_es
      if(iop_es.eq.0) then 
        write(6,*) "  -- perturbative calculation"
      else 
        write(6,*) "  -- iterative calculation"
      endif 

      write(6,*) '- Option of analytic continuation (iop_ac):',iop_ac
      if(iop_ac.eq.1) then 
        write(6,*) "  -- RGN method(Rojas, Godby and Needs)"
      else 
        write(6,*) "  -- Pade's approximation "
      endif

      ! whether shift the Fermi energy during self-consistent GW0 
      call loct_parse_int("iop_esgw0",1,iop_esgw0) 

      ! how the do GW0 self-consistent iteration
      call loct_parse_int("iop_gw0",1,iop_gw0) 

!
!     Read the options for the mixed basis functions
!
      call loct_parse_int("MB_nlo",1,nlo_mb)    
      call loct_parse_logical("MB_ludot",.false.,mb_ludot) 
      call loct_parse_float("MB_emin",-1.E10, mb_emin) 
      mb_emin = mb_emin / 2.0    ! change the unit from Ry to Ha
      call loct_parse_float("MB_emax",20.0, mb_emax) 

      ierr=loct_parse_isdef(blk_mixbas)
      if(ierr.eq.1) then 
        call loct_parse_block_float(blk_mixbas,0,0,kmr  ) 
        call loct_parse_block_int  (blk_mixbas,1,0,lmbmax  ) 
        call loct_parse_block_float(blk_mixbas,1,1,wftol) 
        call loct_parse_block_int  (blk_mixbas,1,2,lblmax  ) 
      else 
        kmr = 1.0
        lmbmax = 3
        lblmax = lmbmax*2
        wftol = 1.D-4
      endif

      write(6,*)'Mixed basis parameters'
      write(6,*)'Interstitial: Maximum |G| of ipw in rkmax units'
      write(6,104) kmr
      write(6,*)'MT-Spheres: Maximum l, Linear dependence tolerance'
      write(6,205) lmbmax, wftol

      call loct_parse_int("nspin_mb",1,nspin_mb) 
      if(nspin.eq.2 ) then
        write(6,*) " MixBasis option for the spin-polarized case:" 
        if(nspin_mb.eq.1) then 
          write(6,*) " use only spin-up radial func."
        else 
          write(6,*) " use only both up and dn radial func."
        endif 
      endif 

      call loct_parse_float("barcevtol" ,-1.0d-10,barcevtol)

! Read the option 

!
!  Input for LDA+U based GW calculations 
!
!  %<blk_vorb>
!     <norbu>
!     <_i>  |  <lorb_i>
!     ---
!  %<blk_vorb> 
      ierr=loct_parse_isdef(blk_vorb) 
      if(ierr.eq.1) then 
        write(6,*) " Option for Vorb"
        call loct_parse_block_int(blk_vorb,0,0,natorb)
        if(natorb.gt.0) then 
          write(6,*) "  GW Calculations based on LDA+U"
          lvorb = .true. 
          allocate(iatorb(natorb),nlorb(natorb),lorb(2,natorb), &
     &             vorb(-3:3,-3:3,2,natorb,nspin),              &
     &             dmorb(-3:3,-3:3,2,natorb))
          do ia=1,natorb
            call loct_parse_block_int(blk_vorb,ia,0, iatorb(ia) )
            call loct_parse_block_int(blk_vorb,ia,1,  nlorb(ia) )
            do il=1,nlorb(ia) 
              call loct_parse_block_int(blk_vorb,ia,il+1,lorb(il,ia))
            enddo 
          enddo 
          
        endif 
      else 
        lvorb = .false. 
      endif 

! 
!     Read the parameters for the Bare coulomb potential
!

      ierr=loct_parse_isdef(blk_bcoul)
      if(ierr.eq.1) then 
        call loct_parse_block_float(blk_bcoul,0,0,pwm)
        call loct_parse_block_float(blk_bcoul,0,1,stctol)
      else 
        pwm=2.0
        stctol = 1.E-15
      endif 

      !! set the trancation radius for the bare Coulomb interaction, needed for finite systems
      call loct_parse_int  ("iop_coulvm",  0,iop_coulvm )
      call loct_parse_int  ("iop_coul_x",  0,iop_coul_x )
      call loct_parse_int  ("iop_coul_c",  0,iop_coul_c )

!      if(iop_coul.ne.0.and.iop_coulvm.eq.0) then
!        write(6,*) "WARNING: For truncated/screened Coulumb interaction,&
!     & Coulomb matrix must be calculated by the plane-wave expansion&
!     & scheme (i.e. iop_coulvm=1)!"
!        iop_coulvm = 1
!      endif    

      call loct_parse_float("rcut_coul",-1.0d0,rcut_coul)

      write(6,100) 'Parameters for Coulomb matrix:'
      write(6,200) "  Maximum |G| in kmr units = ",pwm
      write(6,200) "  Error tolerance for struc. const = ",stctol
      write(6,101) "  Coulomb interaction for exchange",   iop_coul_x 
      write(6,101) "  Coulomb interaction for correlation",iop_coul_c 

      write(6,*)'------------------------------------------------------'
      
   10 format(a4,' or ',a4,' ---> ',a)
   11 format(a5,' or ',a5,' ---> ',a)
   12 format(a6,' or ',a6,' ---> ',a)
  100 format(2x,a)
  101 format(2x,a,i5)
  102 format(2x,a,2i5)
  103 format(2x,i4,1x,i4,1x,2f10.5)
  104 format(f10.5)
  106 format(i4)
  108 format(a4,1x,a4)
  109 format(a5)
  110 format(a3)
  111 format(a3,1x,'electron calculation')
  200 format(4x,a,f12.6)
  205 format(i4,2x,1p,g12.5)
      end subroutine readingw
      
!EOC      
