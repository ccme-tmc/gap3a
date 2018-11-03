!BOP
!
! !ROUTINE: task_evec
!
! !INTERFACE:
      subroutine task_evec

! !DESCRIPTION:
!
! This subroutine calculates the real space representation of the
! (L)APW+lo basis functions in the line joining the two given atoms 
! for ploting
!
! !USES:

      use bands
      use constants,  only: czero,pi
      use eigenvec
      use kpoints
      use lapwlo,     only: nt,lomax,nlomax
      use mixbasis,   only: maxbigl,nmix,bigl,umix
      use param
      use radwf
      use recipvec,   only: indgk, gindex
      use rspevec,     only: evecmt,evecmts,rg,rd
      use selfenergy
      use struk
      use liboct_parser
      use task


! !LOCAL VARIABLES:

      implicit none
      integer(4) :: at1,at2,iik,ibmin,ibmax,ierr
      
      integer(4) :: iat,j,latom

      integer(4) :: ib ! THe band index of the function
      
      integer(4) :: i,igp,ib1,ib2     
      integer(4) :: iat1,iat2,idf,idf1,idf2
      
      integer(4) :: ikvec(3),igvec(3)
      
!      integer(4) :: ikvec(3)

!      real(8) :: kvec(3)    ! Coordinates of the k-point.
      
      real(8) :: rdlen,ri(3),gvec(3),phs,rr(99),rtem(3),rdrot(3)
      
      real(8) :: kvec(3),phsat
      real(8) :: rdi(3),irlen,rd2(3),kgveci(3),kgvec(3)
      
      complex(8) :: pw(99)
      complex(8), allocatable :: yl(:)
      
      character(len=1)  :: extik,extat1,extat2
      character(len=2)  :: extib
      character(len=17) :: filename
      character(len=67) :: errmsg

!     real(8), dimension(3) ::  qvec(3) ! q-vector for which the matrix 
!                                        has to be calculated.
! !DEFINED PARAMETERS:

 
      character(len=8) , parameter :: sname = 'task_evec'        
                                     
 
!
! !EXTERNAL ROUTINES: 
!



      external calcrwf
      external evec2rspmt
      external expand_evec
      external k2cart
      external outerr
      external radmesh
      external readvector
      external rotate
      external unrotate
      external ylm

 
! 
! !INTRINSIC ROUTINES: 
!


      intrinsic achar
      intrinsic dcos
      intrinsic dsin
      intrinsic mod

! 
! !REVISION HISTORY:
!
! Created: 9th. July 2004 by RGA
! Last modified: 15th. July 2004 by RGA
!
!EOP
!BOC

      call linmsg(6,'-',"task:evec")
      ierr=loct_parse_isdef(taskname)
      if(ierr.ne.0) then
        call loct_parse_block_int(taskname,0,0,iik)
        call loct_parse_block_int(taskname,0,1,ibmin)
        call loct_parse_block_int(taskname,0,2,ibmax)
        call loct_parse_block_int(taskname,0,3,at1)
        call loct_parse_block_int(taskname,0,4,at2)
      else
        iik=1
        ibmin=nomax
        ibmax=numin
        at1=1
        at2=2
      endif

!
!     Set the name of the output file
!
      extik=achar(iik+48)
      extat1=achar(at1+48)
      extat2=achar(at2+48)
!
!     Loop over eigenvectors
!
      do ib = ibmin, ibmax     
      
        ib1=mod(ib,10)
        ib2=(ib-ib1)/10
        extib=achar(ib2+48)//achar(ib1+48)
        filename='evec-'//extik//'-'//extib//'-'//extat1//'-'//         &
     &           extat2//'.out'

!
!       open output file
!       
        open(unit=71,file=filename,status='unknown')
        
!
!       Set the indexes iat and idf of the two atoms
!
        iat1=0
        iat2=0
        if(at1.eq.at2)then
          rdi(1)=1.0d0
          rdi(2:3)=0.0d0
          latom=0
          do iat = 1,nat
            do idf=1,mult(iat)
              latom=latom+1
              if(latom.eq.at1)then
                idf1=idf
                iat1=iat
              endif  
            enddo 
          enddo  
          iat2=iat1
          idf2=idf1
        else  
          rdi(1:3)=pos(1:3,at2)-pos(1:3,at1)
          latom=0
          do iat = 1,nat
            do idf=1,mult(iat)
              latom=latom+1
              if(latom.eq.at1)then
                idf1=idf
                iat1=iat
              endif  
              if(latom.eq.at2)then
                idf2=idf
                iat2=iat
              endif
            enddo
          enddo 
        endif  
        if((iat1.gt.0).and.(iat2.gt.0))then
!
!         Allocate the arrays for the Alm, Blm and Clm coefficients
!
          allocate(yl((nt+1)*(nt+1)))
!
!         Read the reciprocal vectors corresponding to the given k-point
!
          call readvector(iik,1)
!          
!          calculate the coordinates of the reciprocal vector
!         
          ikvec(1:3)=klist(1:3,iik)
          call k2cart(ikvec,idvk,kvec)
!
!         calculate the coefficients Alm, Blm, Clm
!
          call expand_evec(iik,1,.true.)
! 
!         Set the cartesian coordinates of the vector
!          
!          do i=1,3
!            if(abs(rdi(i)).gt.0.5)rdi(i)=rdi(i)-sign(1.0d0,rdi(i))
!          enddo  
          do i=1,3
            rd(i)=rbas(1,i)*rdi(1)+rbas(2,i)*rdi(2)+rbas(3,i)*rdi(3)
          enddo  
          rdlen=sqrt(rd(1)*rd(1)+rd(2)*rd(2)+rd(3)*rd(3))
          write(6,*)'rdi,rd,rlen'
          write(6,*)pos(:,at1)
          write(6,*)pos(:,at2)
          write(6,*)rdi
          write(6,*)rd
          write(6,*)rdlen
!          
!         allocate the necesary arrays for atom 1
!
          call init_radwf(nrpt(iat1),lomax,nt,nlomax,lmax)
          allocate(evecmt(nrpt(iat1)))
          allocate(evecmts(nrpt(iat1)))
          allocate(rg(nrpt(iat1)))
!
!         rotate rd into the internal coordinates of at1
!
          call unrotate(rd,rotloc(1:3,1:3,iat1),rtem)
          call unrotate(rtem,rotij(1:3,1:3,at1),rdrot)
          write(6,*)'rd/rdrot at1'
          write(6,*)rd
          write(6,*)rdrot
!          
!         calculate the values of the spherical harmonics for atom 1
!
          call ylm(rdrot,nt,yl)
!
!         calculate the radial wavefunctions of atom 1
!              
          call calcrwf(iat1)
!          
!         calculate the wave function inside the Muffin tin sphere of atom 1
! 
          call radmesh(iat1,rg)
          call evec2rspmt(iat1,at1,ib,yl,1)
!          
!          write to file
!
          do i=1,nrpt(iat1)
            write(71,*)rg(i),real(evecmt(i)),aimag(evecmt(i))
          enddo
          if(iat1.ne.iat2)then
            call end_radwf
            deallocate(evecmt)
            deallocate(evecmts)
            deallocate(rg)
          endif  
!
!         Calculate the phase of the plane waves due to the change of origin
!         
          irlen=rdlen-rmt(iat1)-rmt(iat2)
          do i=1,99
            rr(i)=rmt(iat1)+dble(i)*irlen/1.0d+2
          enddo  
          pw(1:99)=czero
          do igp = 1, nv(iik)
            igvec(1:3)=gindex(:,indgk(igp,iik)) 
            call k2cart(igvec,1,gvec)
            do i=1,3
              kgveci(i)=dble(igvec(i))+dble(ikvec(i))/dble(idvk)
              kgvec(i)=gvec(i)+kvec(i)
            enddo  
            do i=1,99
              phsat=2.0d0*pi*(kgveci(1)*pos(1,at1)+kgveci(2)*         &
     &                       pos(2,at1)+kgveci(3)*pos(3,at1))
              ri(1:3)=rd(1:3)*rr(i)/rdlen
              phs=phsat+(kgvec(1)*ri(1)+kgvec(2)*ri(2)+kgvec(3)*ri(3))
              pw(i)=pw(i)+zzk(igp,ib)*sqrt(vi)*cmplx(dcos(phs),dsin(phs),8)
            enddo  
          enddo  
          do i=1,99
            write(71,*)rr(i),real(pw(i)),aimag(pw(i))
          enddo  
!          
!         allocate the necesary arrays for atom 2
!
          rd2(1:3)=-1.0d0*rd(1:3)
!
!         rotate rd into the internal coordinates of at1
!
          call unrotate(rd2,rotloc(1:3,1:3,iat2),rtem)
          call unrotate(rtem,rotij(1:3,1:3,at2),rdrot)
          write(6,*)'rd/rdrot at2'
          write(6,*)rd2
          write(6,*)rdrot
          
!
!         calculate the radial wavefunctions of atom 1
!         
          if(iat1.ne.iat2)then
            call init_radwf(nrpt(iat2),lomax,nt,nlomax,lmax)
            allocate(evecmt(nrpt(iat2)))
            allocate(evecmts(nrpt(iat2)))
            allocate(rg(nrpt(iat2)))
            call calcrwf(iat2)
          endif  
!          
!         calculate the values of the spherical harmonics for atom 1
!
         
!          call ylm(rd,nt,yl)
!          
!          do i=1,nt*nt
!            write(72,*)i,real(yl(i)),aimag(yl(i))
!          enddo  
!          
          call ylm(rdrot,nt,yl)
          
!          do i=1,nt*nt
!            write(73,*)i,real(yl(i)),aimag(yl(i))
!          enddo  
!          
!         calculate the wave function inside the Muffin tin sphere of atom 1
! 
          call radmesh(iat2,rg)
          call evec2rspmt(iat2,at2,ib,yl,1)
!          
!          write to file
!          
          do i=1,nrpt(iat2)
            j=nrpt(iat2)-i+1
            write(71,*)rdlen-rg(j),real(evecmt(j)),aimag(evecmt(j))
          enddo
          write(71,*)
          call end_radwf
          deallocate(yl)
          deallocate(evecmt)
          deallocate(evecmts)
          deallocate(rg)
          deallocate(almr,almi)
          deallocate(blmr,blmi)
          deallocate(clmr,clmi)
        else   
          write(errmsg,*)'Atom not found'
          call outerr(sname,errmsg)
        endif
        
        close(71)
      
      enddo ! ib  
      
      return
  
      end subroutine task_evec  
!EOC      
