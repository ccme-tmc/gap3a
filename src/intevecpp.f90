!BOP
!
! !ROUTINE: intevecpp
!
! !INTERFACE:
      subroutine intevecpp(ikp,jkp,ib1,ib2)

! !DESCRIPTION:
!
! This subroutine calculates the integral of the square
! modulus of the product of  two LAPW eigenvectors
! !USES:

      use angintegrals, only: n_angrid, w_angrid, sphar
      use bands,        only: nv
      use constants,    only: czero,pi,cfein1,cfein2
      use lapwlo,       only: nt, lomax, nlomax, lmax
      use eigenvec,     only: zzk, zzq
      use mixbasis,     only: ipwint
      use radwf,        only: init_radwf, end_radwf
      use recipvec,     only: gindex, indgk, ig0
      use rspevec,      only: evecmt,evecmts,rg,rd
      use struk,        only: nat, nrpt, mult, vi, ndf
      use task,         only: errfn

! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: ikp ! The k-point for which the (L)APW+lo
!                                     function is ploted
      
      integer(4), intent(in) :: jkp ! The k-point for which the (L)APW+lo
!                                     function is ploted
      
      integer(4), intent(in) :: ib1 ! THe band index of the function

      integer(4), intent(in) :: ib2 ! THe band index of the function

      
! !LOCAL VARIABLES:

      integer(4) :: iat,idf,ngr,ileb
      
      integer(4) :: i,i1
      integer(4) :: ieq,ik1,ik2,ik3,ik4
      
      integer(4) :: ig(3)
      
      real(8) :: epsq,epsqs
      
      real(8) :: epii,intep
      real(8) :: epint(ndf),ig1(3),ig2(3)
      real(8), allocatable :: rs(:),epai(:),epsai(:)
      
      complex(8) :: intmix
      complex(8), allocatable :: yl(:)
      complex(8), allocatable :: evecprod(:)
      complex(8), allocatable :: evecprods(:)
      complex(8) :: fac1,fac2,cepii
      complex(8) :: evecprod1,evecprods1
      
      character(len=67) :: errmsg

! !DEFINED PARAMETERS:
 
        character(len=9) , parameter :: sname = 'intevecpp'        
!
! !EXTERNAL ROUTINES: 
!
      external calcrwf
      external expand_evec
      external evec2rspmt
      external outerr
      external prep_ang_int
      external readvector
      external radmesh
      external rint13
      external rotate
! 
! !INTRINSIC ROUTINES: 
!
      intrinsic conjg
! 
! !REVISION HISTORY:
!
! Created: 9th. July 2004 by RGA
! Last modified: 15th. July 2004 by RGA
!
!EOP
!BOC
!
!     Set the indexes iat and ieq of the two atoms
!
      ngr=16*nt*nt/3
!
!       Allocate the arrays for the Alm, Blm and Clm coefficients
!
      allocate(yl((nt+1)*(nt+1)))
!
!       Read the reciprocal vectors corresponding to the given k-point
!
      call readvector(ikp,1,1)
      call expand_evec(ikp,1,.true.,1)
      call readvector(jkp,2,1)
      call expand_evec(jkp,2,.true.,1)
      idf=0
      intep=0.0d0
      do iat = 1,nat
        allocate(evecprod(nrpt(iat)))
        allocate(evecprods(nrpt(iat)))
!        
!       allocate the necesary arrays for atom 1
!
        call init_radwf(iat)
        allocate(evecmt(nrpt(iat)))
        allocate(evecmts(nrpt(iat)))
        allocate(rg(nrpt(iat)))
        allocate(rs(nrpt(iat)))
        allocate(epai(nrpt(iat)))
        allocate(epsai(nrpt(iat)))
!
!       calculate the radial mesh for atom iat
!       
        call radmesh(iat,rg)        
!
!       calculate the radial wavefunctions of atom 1
!            
        call calcrwf(iat,1)

        do ieq=1,mult(iat)
          idf=idf+1
          epai(:)=0.0d0
          epsai(:)=0.0d0
          do ileb=1,n_angrid
            yl(:)=sphar(ileb,:)
!        
!           calculate the wave function inside the Muffin tin sphere of atom 1
! 
            call evec2rspmt(iat,idf,ib1,yl,1)
            do i=1,nrpt(iat)
              evecprod(i)=evecmt(i)
              evecprods(i)=evecmts(i)
            enddo
!        
!           calculate the wave function inside the Muffin tin sphere of atom 1
! 
            call evec2rspmt(iat,idf,ib2,yl,2)
            do i=1,nrpt(iat)
              evecprod1=evecprod(i)*evecmt(i)
              evecprods1=evecprods(i)*evecmts(i)
              epsq=real(evecprod1*conjg(evecprod1))
              epsqs=real(evecprods1*conjg(evecprods1))
              epai(i)= epai(i)+4.0d0*pi*w_angrid(ileb)*epsq
              epsai(i)= epsai(i)+4.0d0*pi*w_angrid(ileb)*epsqs
            enddo
          enddo
          do i=1,nrpt(iat)
            rs(i)=rg(i)*rg(i)
          enddo
          call rint13(cfein1,cfein2,epai,epsai,rs,rs,epint(idf),iat)
          intep= intep+epint(idf)
        enddo 
        deallocate(evecprod)
        deallocate(evecprods)
        call end_radwf
        deallocate(evecmt)
        deallocate(evecmts)
        deallocate(rg)
        deallocate(rs)
        deallocate(epai)
        deallocate(epsai)
      enddo   ! iat
      idf=0
      i1=0
      intmix=czero
      do iat=1,nat
        do ieq=1,mult(iat)
          idf=idf+1
          write(71,11)ib1,ib2,idf,epint(idf)
        enddo
      enddo     
      cepii=czero
      do ik1=1,nv(ikp)
        do ik2=1,nv(jkp)
          ig1(1:3)=gindex(1:3,indgk(ik1,ikp))-gindex(1:3,indgk(ik2,jkp))
          fac1=vi*zzk(ik1,ib1)*conjg(zzq(ik2,ib2))
          do ik3=1,nv(ikp)
            ig2(1:3)=ig1(1:3)-gindex(1:3,indgk(ik3,ikp))
            fac2=fac1*conjg(zzk(ik3,ib1))
            do ik4=1,nv(jkp)
              ig(1:3)=ig2(1:3)+gindex(1:3,indgk(ik4,jkp))
              cepii=cepii+fac2*zzq(ik4,ib2)*ipwint(ig0(ig(1),ig(2),ig(3)))
            enddo
          enddo        
        enddo
      enddo
      epii=real(cepii)
      intep=intep+epii
      write(72,10)ib1,ib2,epii
      write(73,10)ib1,ib2,intep
      deallocate(yl)

   10 format(2i4,1d15.7)
   11 format(3i4,1d15.7)
      return
      
  900 call outerr(sname,'The kpoint nr. given is not in the list')
  
      end subroutine intevecpp  
!EOC



