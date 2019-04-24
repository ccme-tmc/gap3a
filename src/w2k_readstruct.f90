!BOP
!
! !ROUTINE: w2k_readstruct
!
! !INTERFACE: 
      subroutine w2k_readstruct
!
! !DESCRIPTION:
!
! {\bf WIEN2k interface:}
!
! Reads the structure from file case.struct.
!
! 
! !USES:
      use constants, only: pi,cfein1,cfein2
      use radwf,     only: rel 
      use struk,     only: nat,ndf,alat,alpha,iatnr,mult,pos,vi,dh,nrpt,&
     &                     ro,rmt,rmtmin,vmt,zz,lattic,atomname,nsym,   &
     &                     imat,inddf,tau,izmat,rotloc,init_struk
      use task,      only: casename

! !LOCAL VARIABLES:

      implicit none

      integer(4) :: i
      integer(4) :: iat      ! (Counter) runs over inequivalent atoms
      integer(4) :: idf      ! (Counter) runs over atoms (including  equivalent ones)
      integer:: ierr
      integer(4) :: ieq      ! (Counter), run over equivalent atoms.
      integer(4) :: isym     !(Counter) Runs over symmetry operations
      integer(4) :: j
      integer(4) :: isplit
      integer(4) :: fid

      character(len=4)  :: irel   ! 'rela'/'nrel' for scalar-relativistic/non-relativistic  calculations,
      character(len=80) :: title  ! name of the job
      character(len=80) :: msg
!
! !DEFINED PARAMETERS:

      character(len=12), parameter :: sname = 'w2k_readstruct' ! name of the  subroutine, to include it in error messages.
!EOP
!BOC
  
      call linmsg(6,'-',sname)
      fid=999
      open(unit=fid,file=trim(casename)//".struct",action='read',iostat=ierr)
      call errmsg(ierr.ne.0,sname,"Failt to open struct file")

!
!     read the title
!
      read(fid,5000) title

      write(6,*) title 
!
!     read the type of lattic, number of equivalent atosm and
!     relativistic option
!
      read(fid,5010) lattic, nat, irel
      write(6,'(A,A4)') "Lattice type    ", lattic
      write(6,'(A,I4)') "Number of atoms ", nat
      write(6,'(A,A4)') "Rel. or not     ", irel

!
!     read  lattice constants and angles in degrees
!     transform angles to radians
!
      read(fid,5020) alat(1:3),alpha(1:3) 

      write(6,'(A,6F12.6)') " Lattice constants: ",alat(1:3),alpha(1:3) 
     
      if (alpha(1) .eq. 0) alpha(1) = 90.0d+0
      if (alpha(2) .eq. 0) alpha(2) = 90.0d+0
      if (alpha(3) .eq. 0) alpha(3) = 90.0d+0
      if(trim(lattic).eq.'H'.and.alpha(3).ne.120.0) then 
        alpha(3)=120.0
      endif 

      alpha(1) = alpha(1)*pi/180.0d0
      alpha(2) = alpha(2)*pi/180.0d0
      alpha(3) = alpha(3)*pi/180.0d0
!
!     set the flag for relativistic calculations  
!     and set the corresponding fine structure constants
!
      if (irel .eq. 'RELA') rel = .true.
      if (irel .eq. 'NREL') rel = .false.
      cfein1 = 1.
      if (rel) then
        cfein2 = 137.0359895**(-2)
      else
        cfein2= 1.e-22
      endif

!
!     Allocate necesary arrays
!
      call init_struk(0)

!
!     read crystal-structure (atompositions, symmetry-parameters,
!       muffin-tin radius, ...)
!
      idf = 0
      do iat = 1,nat    !! loop over inequivalent atoms
        idf = idf + 1
        inddf(idf) = iat
        read(fid,5030) iatnr(iat),(pos(i,idf),i=1,3),mult(iat),isplit  !! index and position of noneq. atoms
        write(6,'(I4,3F8.4,I2)') iatnr(iat), pos(1:3,idf),mult(iat)

        do ieq = 1, mult(iat) - 1
          idf = idf + 1
          inddf(idf) = iat
          read(fid,5040) iatnr(iat), (pos(i,idf),i=1,3)  !! index and position of each equivalent atom
        enddo
!
!       read name, number of radial mesh points (nrpt), minimum mesh 
!       point (ro), muffin tin radius (rmt) and number of electrons (yy)
!       of each inequivalent atom
!
        read(fid,5050) atomname(iat),nrpt(iat),ro(iat),rmt(iat),zz(iat)

        dh(iat)  = log(rmt(iat)/ro(iat))/(nrpt(iat) - 1)      !* set the logarithmic step for the radial mesh
        rmt(iat) = ro(iat)*exp(dh(iat)*(nrpt(iat)-1))         !* calc. the muffin tin radius

        read(fid,5060) ((rotloc(i,j,iat),i=1,3),j=1,3)        !* local rotation matrices
      enddo
      ndf = idf 
      write(6, *) "Number of all atoms: ", ndf
!
!     reallocate the arrays to their actual sizes
!
      call init_struk(1)
!
!     read symmetry operations and nonprimitive translations
!
      read(fid,5100) nsym
      allocate(imat(3,3,nsym),izmat(3,3,nsym),tau(3,nsym))

      do isym = 1, nsym
        read(fid,5110) ((imat(i,j,isym),i=1,3),tau(j,isym),j=1,3)
      enddo

      close(fid)

      return
!
!     error handling
!
  930 write(msg,9050) iat, mult(iat), idf
      call outerr(sname,msg)
!
! formats for Wien2k.03
!
 5000 format(a80)
 5010 format(a4,23x,i3,/,13x,a4)
 5020 format(6f10.6)
 5030 format(4x,i4,4x,f10.7,3x,f10.7,3x,f10.7,/,15x,i2,17x,i2)
 5040 format(4x,i4,4x,f10.7,3x,f10.7,3x,f10.7)
 5050 format(a10,5x,i5,5x,f10.9,5x,f10.5,5x,f10.5)
 5060 format(20x,3f10.8)
 5100 format (i4)
 5110 format (3 (3i2,f10.7,/))
 6000 format(///,3x,'error in ebnd  : mult(iat)=0 ...', &
             /, 20x,'iat =',i3,3x,'idf=',i3,3x,'mult=',i3)
 6010 format(a4,'lattice, nonequiv atoms',i3,/,'MODE OF CALC=',a4,      &
     & ' unit=',a4)
 6030 format('atom',i4,': x=',f10.7,' y=',f10.7,' z=',f10.7,/,10x,      &
     &'mult=',i2,10x,'isplit=',i2)
 6040 format('atom',i4,': x=',f10.7,' y=',f10.7,' z=',f10.7)
 6050 format(a10,'npt= ',i5,'  r0=',f10.9,' rmt=',f10.5,3x,'z:',f10.5)
 6060 format('Local rot matrix:',3x,3f10.8,2(/,20x,3f10.8))
 6100 format(i4,'      number of symmetry operations')
 6110 format (2(3i2,f10.7,/),3i2,f10.7)
 6111 format(4x,i4)
 9050 format('mult(',i3,'=',i3,', idf =',i3)

      end subroutine w2k_readstruct
!EOC
