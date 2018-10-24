!BOP
!
! !ROUTINE: w2k_readcore
!
! !INTERFACE:
      subroutine w2k_readcore
!
! !DESCRIPTION:
!
! {\bf WIEN2k interface:}
!
! Reads the core states energy and wave functions from the file case.core and sets up the
! corresponding quantum numbers. The number of valence electrons is also
! set here 
! It is assumed both core states energy and wave functions are stored in 
!    case.core 
! which is concatation of case.scfc[up/dn] and case.corewf[up/dn] as
! obtained from wien2k output 
!
! !USES:

      use bands,  only: nspin
      use core,   only: nclmmax,ncore,npqn,kappa,ocmax,eigcore,symbl,      &
     &                  init_core,lcore,corind,clmind,ncg,nclm,      &
     &                  ncoremax,lcoremax,ucore,uscore,nelcor,nelfc
      use kpoints,only: nvel
      use struk,  only: iatnr,nat,ndf,mult,zz,nrpt
      use task,   only: casename
!
! !OUTPUT PARAMETERS:
!

!EOP
!BOC
      implicit none

! 
!
! !LOCAL VARIABLES:
!
      integer:: pqn   ! principal quantum number
      integer:: iat   ! index unequivalent atoms
      integer:: ic      ! (Counter) runs over core states of an atom
      integer:: icg     ! (Counter) runs over all core states
      integer:: ieq     ! (Counter) runs over equivalent atoms.
      integer:: iclm    ! index for core states at one atom
      integer:: idf
      integer:: iorb  ! index unequivalent atoms
      integer:: isp    ! index for spin
      integer:: itmp   ! integer temporal 
      integer:: qn_n,qn_kappa  ! quantum number n and kappa
      integer:: mc      
      integer:: norb   ! number of core states 
      integer:: lc     ! maximum azimutal quantum number of the core states of this atom
      integer:: nclm1,lcoremax1
      integer:: fin, fout
      integer:: ierr
      real(8):: ec,occ,rtmp
      character(3) :: nm
      character(2) :: lsymbl

      character(len=10), allocatable :: atomname(:) ! name of the core state
      real(8),allocatable:: occ_inc(:,:)
      logical::lprt=.false.,l_oldcore
      character(len=20):: sname="w2k_readcore"

#ifdef DEBUG
      lprt = .true.
#endif 

      fin=999
      fout=6 
      call linmsg(fout,'-',trim(sname))
      allocate(atomname(nat),ncore(nat),occ_inc(20,nat))
      lcoremax=0
      open(unit=fin,file=trim(casename)//".core",action='read',iostat=ierr) 
      call errmsg(ierr.ne.0,"w2k_readcore","Failt to open case.core")
!
! Get ncoremax and init core states 
!

      !! check whether the core file is written in new or old style 
      !! the new core file contains a copy of case.inc
      l_oldcore=.false.
      read(fin,fmt=*,iostat=ierr) itmp,rtmp
      if(ierr.ne.0) then 
        write(6,*) "WARNING: old core file is detected!"
        l_oldcore=.true.
      endif 
      rewind(fin) 

      if(l_oldcore) then  !! to be compatible with gwxv1
        occ_inc=1.d0 
        do iat=1,nat
          read(fin,100) norb
          do iorb=1,norb
            read(fin,102) ec 
            if(norb.eq.1.and.abs(ec)<1.e-5) norb = 0  
          end do
          ncore(iat)=norb
          write(6,101) iat,ncore(iat) 
        end do
        rewind(fin) 

      else 
        if(lprt) write(6,*) "read core states occupuation information!"
        do iat=1,nat 
          read(fin,*) norb
          ncore(iat) = norb
          do iorb=1,norb
            read(fin,*) qn_n,qn_kappa,occ
            occ_inc(iorb,iat) = occ
            if(lprt) write(6,103) qn_n,qn_kappa,occ_inc(iorb,iat) 
            if(occ.lt.1.e-2) ncore(iat) = ncore(iat) - 1
          end do 
          if(lprt) write(6,101) iat,ncore(iat)
        enddo
        read(fin,*) itmp
      endif 
 100  format(/,41x,i2,12x)
 102  format(20x,f20.6)
 101  format(' Number of core states at ',i4,'-th atom:',i4)
 103  format(' n=',i4,' kappa=',i4,' occ=',f6.2)

      ncoremax=maxval(ncore)
      if(ncoremax.eq.0) return   ! no core states needs to be read 
      call init_core(nat,nspin)
!
!     Read core states energy 
!

      isp=1   !! spin up
      do iat=1,nat 
        read(fin,200) atomname(iat),norb
        ic = 0 
        do iorb=1,norb
          read(fin,202) pqn,lsymbl,nm,ec

          if(norb.eq.1.and.abs(ec)<1.e-5) cycle 

          if(lprt) write(6,*) pqn,lsymbl,nm,ec,occ_inc(iorb,iat)
          if (occ_inc(iorb,iat) < 1.e-2) then 
            write(6,*) " the core state ",nm," is frozen"
            cycle  
          endif
          ic = ic+1 
          eigcore(ic,iat,isp)=ec*0.5d0
          symbl(ic,iat)=nm 
          npqn(ic,iat)=pqn
          kappa(ic,iat)=l2kappa(lsymbl)
          ocmax(ic,iat)=iabs(kappa(ic,iat))*2
        enddo 
      end do

      if(nspin.eq.2) then !! spin down 
        isp=2
        do iat=1,nat
          read(fin,200) atomname(iat),norb
          ic = 0 
          do iorb=1,norb
            read(fin,204) nm,ec
            if(occ_inc(iorb,iat) .lt. 1.e-2 )   cycle  
            ic = ic + 1 
            eigcore(ic,iat,isp)=ec*0.5d0
          enddo
        enddo
      endif 
 200  format(/20x,a10,11x,i2,12x)              
 202  format(1x,i1,a2,5x,a3,8x,f20.6) 
 204  format(9x,a3,8x,f20.6) 

!
! setup nclm and lcoremax
!
      nclmmax=0
      lcoremax=0
      do iat=1,nat
        nclm1=0
        lcoremax1=0
        do ic=1,ncore(iat) 
          lc=iabs(kappa(ic,iat))-(1-isign(1,kappa(ic,iat)))/2
          nclm1=nclm1+2*lc+1
          lcoremax1=max(lcoremax1,lc)
          lcore(ic,iat) = lc
        enddo 
        nclmmax  = max(nclmmax,nclm1)
        lcoremax = max(lcoremax,lcoremax1)
      enddo

!
! setup corind
!
      allocate(corind(5,nclmmax*ndf),clmind(nclmmax,ndf),nclm(ndf))
      idf=0
      icg=0
      do iat=1,nat
        do ieq=1,mult(iat)
          idf=idf+1
          iclm = 0
          do ic = 1, ncore(iat)
            lc = lcore(ic,iat)
            do mc=-lc,lc
              iclm=iclm+1
              icg=icg+1
              clmind(iclm,idf) = icg
              corind(1,icg)=iat
              corind(2,icg)=idf
              corind(3,icg)=ic
              corind(4,icg)=lc
              corind(5,icg)=mc
            enddo
          enddo
          nclm(idf) = iclm
        enddo
      enddo
      ncg=icg
!
!     Read core states wave functions 
!
      do isp=1,nspin
        do iat = 1, nat
          call sub_readcorewf(iat,isp)
        enddo
      enddo
      close(fin)

!
!    Set the number of valence electrons 
!
      call sub_calcnvel 

!
!     writing out some information 
!
      write(fout,*) " ncoremax= ", ncoremax
      write(fout,*) " Total Number of Core states(ncg): ",ncg 
      write(fout,*) ""
      do iat=1,nat 
        write(fout,721)iat,atomname(iat),ncore(iat)
        do iorb=1,ncore(iat)
          write(fout,731) symbl(iorb,iat),kappa(iorb,iat),            &
     &              lcore(iorb,iat),ocmax(iorb,iat),                  &
     &              eigcore(iorb,iat,1:nspin)
        enddo 
      end do

      deallocate(atomname)
      
 721  format(7x,i2,'.atom ',5x,a10,11x,i2,' core states')
 731  format(a3,' kappa =',i4,'l=',i4,': ocmax =',i2,4x,2f20.6) 

      return

      contains 
      integer(4) function l2kappa(ls)
      character(2)::ls 
      integer(4) ::kp

      select case(trim(ls)) 
      case ("S")
        kp=-1
      case("P")
        kp=-2
      case("PP") 
        kp=1
      case("D")
        kp=-3
      case("DD")
        kp=2
      case("F")
        kp=-4
      case("FF")
        kp=3
      endselect 
      l2kappa=kp
      end function 

      subroutine sub_calcnvel
      implicit none
      integer(4) :: iat
      integer(4) :: ic
      real(8) :: occoreat
      real(8) :: nvelat
      real(8) :: ntot
      ntot=0.0d0
      nelcor=0.d0
      do iat=1, nat
        ntot=ntot+zz(iat)*mult(iat)
        occoreat=0.0d0
        do ic=1,ncore(iat)
          occoreat=occoreat+ocmax(ic,iat)
        enddo
        nelcor= nelcor + mult(iat)*occoreat
      enddo

      if(nvel.lt.0.1) then 
        nvel = ntot - nelcor
      else 
        nelfc = ntot - nvel - nelcor
      endif 

      write(6,*)
      write(6,1) ntot,nvel,nelcor,nelfc
      write(6,*)

    1 format(10x,'Total number of electrons:      ',f8.4,/,&
     &       10x,'Number of valence electrons:    ',f8.4,/,&
     &       10x,'Number of core electrons:       ',i5,/,  &
     &       10x,'Number of frozen core electrons:',i5)
      return

      end subroutine 


      subroutine sub_readcorewf(iat,isp)
!
!     Read core states wave functions
!
      implicit none
      integer,intent(in):: iat,isp

      integer:: lc
      integer:: m     ! indexes the spherical grid point
      integer:: iorb,norb
      integer:: npt
      real(8):: norm

      real(8), allocatable :: uc(:),usc(:)
      character(len=67) :: errmsg ! error message texta

      if(lprt) write(6,*) "read core states on atom ",iat

      npt=nrpt(iat)
      allocate(uc(1:npt),usc(1:npt))

      read(fin,*,err=910) norb
      if(lprt)  write(6,*) " total number of core states ",norb
      ic=0
      do iorb=1,norb
        read(fin,2020,err=910)  nm
        if(lprt) write(6,*) " read orbital ",nm,"occ=",occ_inc(iorb,iat) 
        read(fin,2021,err=910) ( uc(m),m=1,npt)
        read(fin,2021,err=910) (usc(m),m=1,npt)

        if(occ_inc(iorb,iat).lt.1.e-2) then 
          write(6,*) "  exclude core state ", nm 
        else 
          ic = ic + 1 
          lc= lcore(ic,iat)
          norm=sqrt(0.5d0*ocmax(ic,iat)/dble(2*lc+1))
          ucore(1:npt,ic,iat,isp)=uc(:)*norm
          uscore(1:npt,ic,iat,isp)=usc(:)*norm
        endif
      enddo
      deallocate(uc,usc)
      return 

  910 call outerr(sname,'error reading core states wave functions')
 2020 format(17x,a3)
 2021 format(3x,4e19.12)
      end subroutine 
  
      end subroutine w2k_readcore
!EOC
