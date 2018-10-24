!BOP
!
! !ROUTINE: set_struct
!
! !INTERFACE: 
      subroutine set_struct
!
! !DESCRIPTION:
!
! {\bf WIEN2k interface:}
!
! Reads the structure from file case.struct.
!
! 
! !USES:
      use constants, only: pi 
      use radwf,     only: rel
      use struk,     only: nat,vi,rmt,rmtmin,vmt
      use mixbasis,  only: lmbmax
      use task,      only: casename,iop_dftpkg

! !LOCAL VARIABLES:

      implicit none

      integer(4) :: i
      integer(4) :: iat      ! (Counter) runs over inequivalent atoms
      integer(4) :: idf      ! (Counter) runs over atoms (including  equivalent ones)
      integer(4) :: ieq      ! (Counter), run over equivalent atoms.
      integer(4) :: isym     !(Counter) Runs over symmetry operations
      integer(4) :: j
      integer(4) :: isplit
      integer(4) :: fid

      character(len=4)  :: irel   ! 'rela'/'nrel' for scalar-relativistic/non-relativistic  calculations,
      character(len=4)  :: luni   ! Units of length
      character(len=80) :: title  ! name of the job
      character(len=67) :: errmsg ! text of the error message.
!
! !DEFINED PARAMETERS:

      character(len=12), parameter :: sname = 'set_struct' ! name of the  subroutine, to include it in error messages.
!EOP
!BOC
  
      call linmsg(6,'-',"Set Struct")

!
!  Generate Bravais vectors 
!
      if(iop_dftpkg.eq.0) then 
        call w2k_readstruct
        call latgen
!      else 
!        call exc_readstruct
      endif   
!
!     determine minimum rmt
!
      rmtmin = minval(rmt(1:nat))
      do iat = 1, nat
        vmt(iat) = vi*4.0d+0*pi*rmt(iat)**3/3.0d+0
      enddo
      return

      end subroutine set_struct
!EOC
