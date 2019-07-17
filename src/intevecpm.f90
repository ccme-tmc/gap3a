!BOP
!
! !ROUTINE: intevecpm
!
! !INTERFACE:
      subroutine intevecpm(jkp,ib1,ib2,maxmt,releri)

! !DESCRIPTION:
!
! This subroutine calculates the integral of the product of two LAPW
! eigenvectors by integration in the MT sphere and by expandding it in the
!mixed basis
!
! !USES:

      use constants,  only: czero,pi,cfein1,cfein2
      use minmmat,    only: minm
      use mixbasis,   only: maxbigl,nmix,bigl,umix,mbsiz,locmbsiz
      use rspevec,    only: evecmt,evecmts,rg,rd
      use struk,       only: nat, mult, ndf

! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: ib1 ! THe band index of the function

      integer(4), intent(in) :: ib2 ! THe band index of the function

      integer(4), intent(in) :: jkp ! The k-point for which the (L)APW+lo
!                                     function is ploted
      
      real(8), intent(out) :: releri, maxmt 
! !LOCAL VARIABLES:

      integer(4) :: iat,idf
      
      integer(4) :: i1,irm,l1,m1, ii1, ii2, jatom
      integer(4) :: ieq
      
      real(8) :: err, reler, relermt
      
      real(8) :: epii,intep
      real(8) :: epint(ndf)
      
      complex(8) :: intmix,emii
      complex(8) :: emint(ndf)
      
!      character(len=67) :: errmsg

! !DEFINED PARAMETERS:
 
        character(len=9) , parameter :: sname = 'intevecpm'        
                                     
 
!
! !EXTERNAL ROUTINES: 
!
      external calcrwf
      external expand_evec
      external outerr
      external prep_ang_int
      external readvector
      external radmesh
      external rint13
      external rotate

 
! 
! !INTRINSIC ROUTINES: 
!
! 
! !REVISION HISTORY:
!
! Created: 9th. July 2004 by RGA
! Last modified: 15th. July 2004 by RGA
!
!EOP
!BOC
      idf=0
      i1=0
      intmix=czero
      intep=0.0d0
      maxmt=0.0d0
      do iat=1,nat
        do ieq=1,mult(iat)
          idf=idf+1
          emint(idf)=czero
          do irm=1,nmix(iat)
            l1=bigl(irm,iat)  
            do m1=-l1,l1
              i1=i1+1
              emint(idf)=emint(idf)+minm(i1,ib1,ib2)*          &
     &                     conjg(minm(i1,ib1,ib2))
            enddo
          enddo 
          read(71,13)ii1,ii2,jatom,epint(jatom)
          intep=intep+epint(idf)
          if(ii1.ne.ib1)then
            write(*,*)'ii1 = ',ii1,'ib1 = ',ib1
            write(*,*)'ii2 = ',ii2,'ib2 = ',ib2
            write(*,*)'jatom = ',jatom,'idf = ',idf
            call outerr(sname,"ii1.ne.ib1")
          endif  
          if(ii2.ne.ib2)then
            write(*,*)'ii1 = ',ii1,'ib1 = ',ib1
            write(*,*)'ii2 = ',ii2,'ib2 = ',ib2
            write(*,*)'jatom = ',jatom,'idf = ',idf
            call outerr(sname," ii2.ne.ib2")
          endif  

          if(jatom.ne.idf)then
            write(*,*)'ii1 = ',ii1,'ib1 = ',ib1
            write(*,*)'ii2 = ',ii2,'ib2 = ',ib2
            write(*,*)'jatom = ',jatom,'idf = ',idf
            call outerr(sname,"jatom.ne.idf")
          endif  

          err=abs(epint(idf)-real(emint(idf)))
          relermt=err/epint(idf)
          if(abs(relermt).gt.maxmt)maxmt=abs(relermt)
          write(74,11)ib1,ib2,idf,epint(idf),real(emint(idf)),          &
     &                err,relermt   
          intmix=intmix+emint(idf)
        enddo
      enddo     
      read(72,14)ii1,ii2,epii
          if(ii1.ne.ib1)then
            write(*,*)'ii1 = ',ii1,'ib1 = ',ib1
            write(*,*)'ii2 = ',ii2,'ib2 = ',ib2
            call outerr(sname,'ii1.ne.ib1')
          endif  
          if(ii2.ne.ib2)then
            write(*,*)'ii1 = ',ii1,'ib1 = ',ib1
            write(*,*)'ii2 = ',ii2,'ib2 = ',ib2
            call outerr(sname,'ii2.ne.ib2')
          endif  
      intep=intep+epii
      emii=czero
      do i1=locmbsiz+1,mbsiz
        emii=emii+minm(i1,ib1,ib2)*conjg(minm(i1,ib1,ib2))
      enddo  
      err=abs(epii-real(emii))
      releri=err/epii
      write(75,10)ib1,ib2,epii,real(emii),err,releri
      intmix=intmix+emii
      err=abs(intep-real(intmix))
      reler=err/intep
      write(76,10)ib1,ib2,intep,real(intmix),err,reler

   10 format(2i4,4d15.7)
   11 format(3i4,4d15.7)
   13 format(3i4,d15.7)
   14 format(2i4,d15.7)
      return
  
      end subroutine intevecpm  
!EOC



