!BOP
!
! !ROUTINE: task_exhf
!
! !INTERFACE:
      subroutine task_exhf
      
! !DESCRIPTION:
!
! This subroutine calculate the Fock exchange energy 
!
! !USES:

      use hfexch
      use kpoints,     only: nirkp,idikp,wkir
      use bzinteg,     only: kiw,singc2 
      use mixbasis,    only: init_mixbasis,end_mixbasis
      use barcoul,     only: init_barcoul, end_barcoul
      use minmmat,     only: init_minmmat,end_minmmat
      
      
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: iq,iqir        ! (Counter) Runs over q-points
      integer(4) :: iq0, iq2
      integer(4) :: ierr 
      
      real(8) :: tini         ! Initial CPU-time of the subroutine
      real(8) :: tend         ! Final CPU-time of the subroutine
      real(8) :: tsub(2)      ! Initial and final time of each called
!                               subroutine
      

! !EXTERNAL ROUTINES: 

      external addsingint
      external addtoselfe
      external addtoselfx
      external calcpolmat
      external calcvsq
      external calcvxcnn
      external epsq0
      external expand_prods
      external selfesing
      external selfxsing
      

! !INTRINSIC ROUTINES: 


      intrinsic cpu_time      

! !REVISION HISTORY:
!
! Created 16.09.2005 by RGA
!      
!EOP
!BOC     


      call linmsg(6,'-','task:exhf')
      call cpu_time(tini)            
      call init_hfexch(nirkp) 

!
!     Loop over q-points.
!
      do iqir=1,nirkp
        iq = idikp(iqir) 
!
!       Set the size of the basis for the corresponding q-point
!       Allocate the matrices needed for calculating the selfenergy
!       (matsiz dependent)
!
        write(6,*) "iq=",iq
        call init_mixbasis(iq) 
        call init_barcoul(iq)
        call init_minmmat(5)

        call coul_barc(iq, -1) ! no cut-off

!
!        Calculate the Minm matrix elements
!
        call expand_prods(iq)

        call calcexhf(iq) 

!
!       Deallocate arrays with q-dependent sizes
!
        call end_barcoul(iq)
        call end_minmmat
        call end_mixbasis
      enddo ! iq
!
!    Calculate acfd correlation energy by integrating over q-points  
!
      call linmsg(6,'-',"Summary of EXHF ")

      write(6,*) "Contribution from singular integral:",                &
     &           singc2*real(exq0s)

      write(6,'(2A6,2A12)') 'iqir',"wkir",'kiw','exq'
      do iqir=1,nirkp 
        write(6,'(2I6,F10.4,2g16.6)') iqir,wkir(irk),kiw(1,iqir,1),exq(iqir)
      enddo 
      write(6,*) " Fock exchange energy   =", exhf 

      call cpu_time(tend)
      write(6,*)'task_exhf',tend-tini
  101 format(10x,'Data for q-point nr.:',i4,//,10x,'Mixed basis:',/,10x, &
     &           'Number of atomic basis functions:       ',i4,/,10x,   &
     &           'Number of interstitial basis functions: ',i4,/,10x,   &
     &           'Total number of basis functions:        ',i4,/)
      
      return
 1000 format('CPUTIME for ',A20,F16.2,' seconds')
      
      end subroutine task_exhf
!EOC      
