!-----------------------------------------------------------------------
!BOP
!
! !MODULE: acfd
      module acfd
      implicit none 
! !PUBLIC VARIABLES:

      integer(4) :: ifxc     ! control xc kernel used in ACFD calculations  
                             ! ifxc = 0 RPA, integrate over $\lambda$ analytically
                             !      = 1 RPA, integrate over $\lambda$ numerically
                             !      = 2 ALDA kernel 
                             !      = 3 PGG  kernel 
      integer(4) :: ibzint   ! control the BZ integration scheme used in ACFD correlation energy 
                             ! ibzint = 0 

      real(8) :: ex_hf = 0.0D0     ! Hartree-Fock (exact) exchange energy
      real(8) :: ec_acfd = 0.0D0   ! ACFD correlation energy 
      real(8) :: etot_lda = 0.0D0  ! Total energy in LDA/GGA calculation
      real(8) :: exc_lda = 0.0D0   ! LDA/GGA exchange-correlation energy 

      real(8),allocatable :: exq(:)  !! q-dependent HF exchange energy
      real(8),allocatable :: ecwq(:,:)  !! w- and q-dependent ACFD correlation energy 

      integer,private::ierr

      
      contains
      
      subroutine init_acfd(iq0,iq1,iom0,iom1) 
      integer:: iq0,iq1,iom0,iom1  !! ranges for q and \omega index 

      ec_acfd = 0.d0
      ex_hf=0.d0
      allocate(ecwq(iom0:iom1,iq0:iq1),&
     &          exq(iq0:iq1),                 &
     &         stat=ierr)
      if(ierr.ne.0) then
        write(6,*) "init_acfd: Fail to allocate memory"
        stop
      endif

      end subroutine 

      subroutine end_acfd
        deallocate(ecwq,exq)
      end subroutine 

      end module 
!EOP
