!BOP
!
! !ROUTINE: w2k_readin1
!
! !INTERFACE:
      subroutine w2k_readin1

! !USES:

      use bands,   only: nspin
      use lapwlo
      use struk,  only: zz,mult,nat,atomname,rmtmin
      use task,   only: casename
      use modmpi, only: myrank
!
! !INPUT PARAMETERS:

      implicit none

      
! !DESCRIPTION:
!      
! {\bf WIEN2k interface:}
!
! This subroutine reads the file case.in1 and initializes the data
! needed for the radial wave functions
! 
!
! !LOCAL VARIABLES:

      integer(4) :: i,j,k    ! just a counter
      integer(4) :: isp  ! index for spin 
      integer(4) :: iapw ! = 1 for APW's, = 0 for LAPW's. Read from   file case.in1
      integer(4) :: iat  ! Counter: runs over inequivalent atoms
      integer(4) :: l    ! Counter: runs over angular momentum l
      integer(4) :: lalt !  Stores previous value of l (to recognize  repetitions for LO's)
      integer(4) :: nlr  ! number of orbitals for which the default (global) linearization energy (el, e(l,iat)) 
                         ! and wavefunction character (APW+lo or LAPW) are not used. Read from file case.in1
      integer(4) :: nt1  ! maximum l for radial functions      
      integer(4) :: fin ! file id for in1 file 
      integer(4) :: ierr
      
      real(8) :: de      ! Energy increment for searching the resonance  linearization energy (if de = 0 the value is fixed)
      real(8) :: ei      ! expansion energy E_l (Read from file case.in1gw and later stored in e(l,iat)a

      real(8) :: els(nspin) ! store the linearization energy obtained by scan  
      real(8) :: znuc 

      character(len=20), parameter :: sname = 'w2k_readin1'

      character(len=4) :: scanflag ! Used only if de .ne. 0,
                                != CONT: calculation continues, even if either E_{top} or  E_{botton} of the band are not found
                                != STOP: calculation stops if not both E_{top} and E_{botton} of the band are found 

!
! !REVISION HISTORY:
!  Created Feb. 5th. 2004
!
!EOP
!
!BOC      
!rga:
!     It is not necesary to read modus, since if the case.vector file does
!     exist, it will give an error when trying to open it.
!rga:

      call linmsg(6,'-',"w2k_readin1")

      fin=999
      open(unit=fin,file=trim(casename)//".in1",action='read',iostat=ierr) 
      call errmsg(ierr.ne.0,sname,"Fail to open case.in1")
     
      read(fin,*)
      read(fin,*) rkmax, lmax, lnsmax   !* r-mt times k-max,nt and lnsmax
      nt = lmax + 1
      kmax = rkmax/rmtmin

      if(myrank.eq.0) write(6,6020) rkmax, lmax, lnsmax, nat
!
!     Allocate loor,lapw,nlo
!
      call init_lapwlo(nat,nspin)

      do iat = 1, nat

        read(fin,*,err=901) ei,nlr,iapw  !* default parameters ei,iapw and number of exceptions nlr
        if ( iapw .eq. 1 ) lapw(:,iat) = .false.  !* Set default value for lapw(l,iat):

        if(myrank.eq.0) then 
          if(iapw.eq.1)then
            write(6,6000) atomname(iat), ei,'  apw', nlr
          else
            write(6,6000) atomname(iat), ei,' lapw', nlr
          endif
        endif 

!
! Set default value for the linearization energy and local orbital parameters
!
        do l = 0, lmax
          umt(1,l,iat,:) = ei
          if (l.le.lomax) then
            loor(l,iat) = .false.
            nlo(l,iat) = 0
          endif
        enddo 
!
! Read  exceptional cases. 
!   1) If the same l appears twice, the second one is a LO
!         set loor to true, nlo must be incremented by 1
!         and g_tot by (2l+1)*mult(iat)
!   2) If it is an APW+lo (iapw == 1): lapw(l,iat) is set to false, and
!         g_tot has to be incremented, nlo =1 , elo = ei...
!
        l = -1
        do j = 1, nlr
          lalt = l                      !* Set lalt to the previous value of  l
          read(fin,5000,err=901) l, ei, de, scanflag, iapw  

          if(l.eq.lalt)then    ! same l appear for the second time  
            loor(l,iat) = .true.
            nlo(l,iat) = nlo(l,iat) + 1
            nlo_tot = nlo_tot + (2*l+1)*mult(iat)
            if(myrank.eq.0) write(6,6010) l,'Local Orbital'
          else
            if(iapw.eq.1)then
              lapw(l,iat) = .false.
              nlo(l,iat) =  1
              nlo_tot = nlo_tot + (2*l+1)*mult(iat)
              if(myrank.eq.0) write(6,6010)l,'APW+lo'
            else 
              if(myrank.eq.0) write(6,6010) l,'LAPW'
            endif
          endif

        enddo  ! j

      enddo ! iat

      

      nLO_at = 0
      do iat=1, nat
        do l=0,lomax
          if(lapw(l,iat)) then  !! LAPW+LO
            nLO_at(l,iat) = nlo(l,iat)  
          else                  !! APW+lo+LO
            nLO_at(l,iat) = nlo(l,iat) - 1
          endif 
        enddo
      enddo
      nlomax = maxval(nLO_at) 

      !! set lvmax_at
      do iat = 1, nat
        znuc= zz(iat) 
        if (znuc <= 2.0) then 
          lvmax_at(iat) = 0 
        else if (znuc > 2.0 .and. znuc <= 18.0) then 
          lvmax_at(iat) = 1
        else if (znuc > 18.0 .and. znuc <= 54.0) then 
          lvmax_at(iat) = 2
        else 
          lvmax_at(iat) = 3 
        endif 
      enddo 

      if(nlomax.gt.1) then 
        l_newlo = .true. 
        nloat = maxval(nlo) 
        write(6,*) "More than one LO's are detected: set l_newlo"
      endif      

      write(6,*) "Max. number of LO's per l per atom  ( nlomax):",nlomax
      write(6,*) "Max. number of LO/lo's per l per atom (nloat):",nloat

      if(myrank.eq.0) then 
        write(6,*) "Information on Local Orbitals (LO or lo)  "
        write(6,*) "  lomax=",lomax
        write(6,100) (("l=",l,"LO?","#lo","#LO"),l=0,lomax)
        do iat=1,nat
          write(6,102) ((loor(l,iat),nlo(l,iat),nlo_at(l,iat)),l=0,lomax)
        enddo
        write(6,*) ' Total number of local orbitals (nlo_tot)= ',nlo_tot 
      endif 
 100  format(4(a2,i2,a4,a4,a4,8x))
 102  format(4(4x,l4,i4,i4,8x)) 

      close(fin)
      return
!
!     error handling
!      
!     error while readin case.in1
!      
  901 call outerr(sname," error reading file case.in1(c)")
!  
! formats for WIEN2k.03  
!
 5000 format (1x,i1,f10.5,f9.3,1x,a4,i2)
 6000 format (/,10x,'Atomic sphere dependent parameters for atom  ',a10, &
     &        /,10x,'Overall energy parameter is',f10.4,                 &
     &        /,10x,'Overall basis set on atom is',a5,                   &
     &        /,10x,'Number of nonglobal choices',i4)
 6010 format (10x,'l=',i2,4x,a15,2f10.4)
 6020 format(10x,'r-mt times k-max is',f5.2,/10x,'max l is',i3,5x, &
      'max l in nonspherical matrixelements:',i3,/,10x, &
             'number of inequivalent atoms is',i5)

      end subroutine w2k_readin1
!EOC      
