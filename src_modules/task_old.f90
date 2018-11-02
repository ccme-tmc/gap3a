!--------------------------------------------------------------
!BOP
! 
! !MODULE: task
      module task
!
! !DESCRIPTION:
!  this module defines general parameters related to the computational task 
!
      
! !PUBLIC VARIABLES: 
        integer:: iop_dftpkg=0           !! control which DFT package is interfaced 
                                       !!  0 -- wien2k 
                                       !!  1 -- EXCITING  
        integer :: nmax_sc = 10               !! Maximum number of self-consistent iterations 
        real(8) :: eps_sc = 1.0d-3            !! convergence for self-consistent 

        character(len=10)::  progname = 'GAP 2.0'
        character(len=20) :: gwinp="gw.inp"      !! name of the main input file  
        character(len=20) :: gwout="gw.out"      !! name of the main input file  
        character(len=10) :: taskname       !! name of the running task 
        character(len=20)::  casename        !! name of the case to be calculated, 
                                            !! all input and output file will start with this name 
        character(len=120):: scrdir= './'   !! scratch directory
        character(len=120):: scrfn          !! the prefix for scratch files 
        character(len=120):: savdir="."     !! the name for the directory used to save restart data 
        character(len=120):: save_prefix    !! the prefix for files to be saved 

        character(2):: spflag(2)            !! spin flag ['up'/'dn',''] used for I/O files
        integer :: iop_scratch = 2     ! control how to use the scratch space
                                       !   0 -- no scratch
                                       !   1 -- different vectord file for different processes
                                       !   2 - a single vectord for all processes 


        logical :: lrestart=.false.     !! if true, output from previous calculations, if available, will be used  
        logical :: l_debug=.true.        !! 
        logical :: l_UseSavedMinm =.true.
        logical :: l_save_dielmat=.false.  !! indicate whether save the dielmatrix 
        logical :: l_newlo = .false. 

        logical :: ldbg_bzint = .true. 
        logical :: ldbg_qpe   = .false.

!
! record the cpu time used for the most time-consuming subroutines 
!
        real(8) :: time_wall   = 0.d0     ! total wall-time
        real(8) :: time_minm   = 0.d0     ! minm matrix
        real(8) :: time_mwm    = 0.d0     ! M*W*M matrix
        real(8) :: time_eps    = 0.d0     ! dielectrix matrix         (including the calculation of minm) 
        real(8) :: time_selfx  = 0.d0     ! calculation of selfenergy (from given dielectric matrix) 
        real(8) :: time_selfc  = 0.d0     ! calculation of selfenergy (from given dielectric matrix) 
        real(8) :: time_lapack = 0.d0     ! cpu time for calling LAPACK+BLAS subroutines 
        real(8) :: time_evec   = 0.d0     ! cputime used for expanding eigenvectors
        real(8) :: time_coul   = 0.d0     ! cputime used for bare Coulomn matrix 
        real(8) :: time_vxc    = 0.d0     ! cputime used for LDA/GGA Vxc matrix  

!
! memory usage information for largest matrices
!
        real(8) :: mem_minm=0.0
        real(8) :: mem_eps=0.0 
        real(8) :: mem_zzk=0.0

!
! set a set of file IDs to be used for global I/O
        integer,parameter ::fid_outgw   = 6,  & ! general output 
     &                      fid_outmb   = 7,  & ! mixed basis    
     &                      fid_outmom  = 8,  & ! momentum matrix
     &                      fid_outdbg  = 9,  & ! general debug info
     &                      fid_outkpt  =10,  & ! k-point and bzint
     &                      fid_outqp   =11     ! quasi-particle output

!EOP
! 

      end module task 
