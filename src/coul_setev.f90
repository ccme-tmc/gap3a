!BOP
!
! !ROUTINE: coul_setev
!
! !INTERFACE:
      subroutine coul_setev(iq,evtol,icoul)

! !DESCRIPTION:
!
! This subroutine diagonalizes the coulomb matrix and reset the eigenvectors and eigenvalues of 
! in terms of evtol 
!
! !USES:

      use barcoul,    only: vmat,ev,barcvm,barcev,barcevsq,rcut_coul,im_g0
      use constants,  only: cone, czero, imag, pi,twopi,fourpi 
      use mixbasis,   only: mbsiz,matsiz,wi0
      use recipvec,   only: ngq
      use task,       only: fid_outdbg, fid_outgw
!
! !INPUT PARAMETERS: 

      implicit none
      integer, intent(in) :: iq    ! index of the q-point
      real(8), intent(in) :: evtol ! tollerance to screen eigenvalues
      integer, intent(in) :: icoul ! option to treat Coulomb interaction
!
! !LOCAL VARIABLES:

      integer :: im,jm,iipw,ipw      ! Indexes the mixed wave function for the columns of barc (barc(im,jm))
      integer :: immax
      real(8) :: test1,test2,qgeff,rc2inv,ev_new
      integer :: im_kept(mbsiz)      ! indicate which barc eigenvectors are kept as basis functions 
      complex(8) :: wi0new(mbsiz)    ! eigenvalues of sqrt(barc)

!! debug 
      integer :: ierr
      integer :: fdb = fid_outdbg
      logical :: ldbg = .false.
      character(len=20)::sname='coul_setev'
   
!
! !REVISION HISTORY:
! Created July 31,2009 by Hong Jiang
!
!EOP
!BOC

      if(ldbg) call linmsg(6,'-','coul_setev')

!
!    Reduce the basis size by choosing eigenvectors of barc with eigenvalues larger than 
!    evtol
!
      if(ldbg) write(fid_outgw,'(a,e12.4)') " Choose eigenvectors of bare Coulomb matrix &
    & with eigenvalues larger than evtol=",evtol
!
!     Construct new basis set by removing all eigenvectors 
!     with eigenvalues smaller than evtol  
!
      rc2inv = 1.d0/rcut_coul**2
      if(icoul.gt.-1) write(fid_outgw,*) "rc2inv=",rc2inv 

      im_kept = 1
      matsiz=mbsiz
      do im=1,mbsiz
        if( ev(im).lt.evtol ) then 
          im_kept(im)= 0
          matsiz = matsiz - 1
        endif 
      enddo 

      if(iq.eq.1) then
        call coul_wmix0
        call zgemv('c',mbsiz,mbsiz,cone,vmat,mbsiz,wi0,1,czero,wi0new,1)

        !! find the index of the diagonalized barc eigenvector that has maximal overlapw with
        !! G=0 (constant) plane wave
        test2=0.0d0
        do im=1,mbsiz
          test1= real(wi0new(im)*conjg(wi0new(im)))
          if(test1.gt.test2)then
            immax=im
            test2=test1
          endif
        enddo
        if(ldbg) then 
          write(fid_outgw,*)'- Maximum singular eigenvector ###'
          write(fid_outgw,100) immax,test2,ev(immax)
        endif 

        if(icoul.eq.-1) then         !! bare Coulomb interaction 
          if(im_kept(immax).eq.1) then 
            im_kept(immax) = 0
            matsiz = matsiz-1
          endif 
        else                           !! truncated/screened Coulomb interaction 
          if(im_kept(immax).eq.0) then 
            im_kept(immax) = 1
            matsiz = matsiz + 1
          endif 
        endif 
      endif

  100 format("immax,max(wi0new),ev(immax)=",i4,f8.3,e10.3)



      write(fid_outgw,*) "  - Old/New basis set size =",mbsiz,matsiz    
      allocate(barcvm(mbsiz,matsiz),                  &
      &        barcev(matsiz),                        &
      &        barcevsq(matsiz),                      &
      &        stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate barcvm etc.")
      barcvm = 0.d0 
      barcev = 0.d0 
      barcevsq = 0.d0 

      im=0
      do jm=1,mbsiz 
        if(im_kept(jm).eq.1) then 
          im=im+1
          ! - zmy start
          ! We cannot decompose qgeff to q_para and q_perp.
          ! We should not rely on this kind of update for 1D 2D
          ! - zmy end
          if(icoul.gt.-1) then 
            ! TODO how im_g0 is set for icoul = -1
            !! reset eigenvalues for icoul > -1 (truncated/screened Coulomb interaction)
            if(iq.eq.1.and.jm.eq.immax) then 
              if(icoul.eq.0) then     !! truncated Coulomb interaction for 0D
                ev_new = twopi*rcut_coul**2
              elseif(icoul.eq.3) then !! Thomas-Fermi screened Coulomb interaction
                ev_new = fourpi*rcut_coul**2
              elseif(icoul.eq.4) then !! erfc-screened Coulomb interaction
                ev_new = pi*rcut_coul**2
              endif
              im_g0 = im 
            else 
              qgeff = sqrt(fourpi/ev(jm))
              if(icoul.eq.0) then 
                ev_new = ev(jm)*(1.d0-cos(qgeff*rcut_coul))
              elseif(icoul.eq.3) then 
                ev_new = fourpi/(qgeff**2+rc2inv)
              elseif(icoul.eq.4) then 
                ev_new = ev(jm)*(1.d0-exp(-(qgeff*rcut_coul/2.0)**2))
              endif 
              if(ldbg) write(fid_outdbg,'(a,i5,4f12.6)') "i,ev,ev_new,qg",&
             &           im,ev(jm),ev_new,qgeff,cos(qgeff*rcut_coul)  
            endif 
          else 
            ev_new = ev(jm)
          endif ! icoul.gt.-1
          barcev(im)=ev_new 
          barcevsq(im)=sqrt(barcev(im))
          barcvm(:,im)=vmat(:,jm)*barcevsq(im) 
        endif ! im_kept(jm).eq.1
      enddo ! jm  

      end subroutine 
