!BOP
!
! !ROUTINE: set_mixbasis
!
! !INTERFACE:
      subroutine set_mixbasis

! !DESCRIPTION:
!
! This subroutine generate the mixed basis functions 
!
! !USES:

      use bands,      only: nspin 
      use core,       only: lcoremax,ncoremax
      use lapwlo,     only: nt,lmax,lomax,nlomax,lvmax_at
      use mixbasis,   only: lmbmax,lblmax,nmixmax,maxbigl,lmixmax,nmix,bigl, &
     &                      mbl,locmatsiz, init_mixfun,reinit_mixfun,   &
     &                      locmixind,nmixlm,init_mixmat,  &
     &                      init_s3r,nspin_mb,lmbmax_at,lblmax_at 
      use radwf,      only: nrad,init_radwf
      use struk,      only: nat,mult,ndf
      use task,       only: fid_outmb
      use liboct_parser
!
! !LOCAL VARIABLES:
!
      implicit none
      
      integer:: iat,ieq,idf  ! Index for inequivalent/equvalent/all atoms
      integer:: irm,imix,im  ! Index mixed basis functions.
      integer:: isp
      integer:: jlo      ! (Counter) runs over kind of local  orbitals (nlo)
      integer:: jri      ! Number of radial mesh points (for each atom  the corresponding nrpt(iat) is stored here).
      integer:: l,m        ! (Counter) runs over angular momentum l
      integer:: lms   
      integer:: maxnt    ! Maximum l for gaunt coefficients
      integer:: ierr
      logical :: lprt =.false.
      character(len=20):: sname="set_mixbasis"
      character(len=20):: blk_lmbmax="MB_lmbmax"


!
! !EXTERNAL ROUTINES: 
!

!
! !REVISION HISTORY:
! Created Nov. 2003, 
! Last Modified, 11. August 2005 by RGA 
!
!EOP
!BOC
!
!
!
      call linmsg(6,'-',sname)

!
!     set atom-specific lmbmax and lblmax (lmbmax_at and lblmax_at)
!
      allocate(lmbmax_at(nat),lblmax_at(nat),stat=ierr)
      call errmsg(ierr.ne.0,sname,"fail to allocate lmbmax_at...")

!      !! if the blk_lmbmax exists, then set the lmbmax and lblmax in
!      !  terms of the block input  
!      ierr=loct_parse_isdef(blk_lmbmax)
!      if(ierr.eq.1) then 
!        do iat=1,nat 
!          call loct_parse_block_int(blk_lmbmax,0,iat-1,lmbmax_at(iat))
!          call loct_parse_block_int(blk_lmbmax,1,iat-1,lblmax_at(iat))
!        enddo 
!      else a

!* set lmbmax_at
! 
      if(lmbmax > 10) then
        if( lmbmax < 99 ) then                 
          lmbmax_at(1) = lmbmax/10
          lmbmax_at(2) = mod(lmbmax,10)
        elseif (lmbmax < 999) then
          lmbmax_at(1) = lmbmax/100 
          lmbmax_at(2) = (lmbmax-lmbmax_at(1)*100)/10
          lmbmax_at(3) = mod(lmbmax-lmbmax_at(1)*100,10) 
        elseif (lmbmax < 9999) then
          lmbmax_at(1) = lmbmax/1000
          lmbmax_at(2) = (lmbmax-lmbmax_at(1)*1000)/100
          lmbmax_at(3) = (lmbmax-lmbmax_at(1)*1000 -lmbmax_at(2)*100)/10
          lmbmax_at(4) = mod(lmbmax-lmbmax_at(1)*1000-lmbmax_at(2)*100,10) 
        else
          write(6,*) "WARNING: unsupported option for lmbmax=",lmbmax
          write(6,*) " -- set it to the default lmbmax = -1)"
          lmbmax = -1 
        endif 

      endif 
         
      if(lmbmax < 10) then 
        if(lmbmax.le.0) then 
          lmbmax_at = lvmax_at + abs(lmbmax) 
        elseif(lmbmax.gt.0) then 
          lmbmax_at = lmbmax 
        endif  
      endif 

      !! set lblmax 
      if(lblmax > 10) then
        if( lblmax < 99 ) then
          lblmax_at(1) = lblmax/10
          lblmax_at(2) = mod(lblmax,10)
        elseif (lblmax < 999) then
          lblmax_at(1) = lblmax/100
          lblmax_at(2) = (lblmax-lblmax_at(1)*100)/10
          lblmax_at(3) = mod(lblmax-lblmax_at(1)*100,10)
        elseif (lblmax < 9999) then
          lblmax_at(1) = lblmax/1000
          lblmax_at(2) = (lblmax-lblmax_at(1)*1000)/100
          lblmax_at(3) = (lblmax-lblmax_at(1)*1000 -lblmax_at(2)*100)/10
          lblmax_at(4) = mod(lblmax-lblmax_at(1)*1000-lblmax_at(2)*100,10)
        else
          write(6,*) "WARNING: unsupported option for lblmax=",lblmax
          write(6,*) " -- set it to the default lblmax = 0)"
          lblmax = 0
        endif
      endif


      if(lblmax.le.0) then
        lblmax_at = lmbmax_at*2
      else
        lblmax_at = lblmax
      endif
      
      lmbmax = maxval(lmbmax_at)
      lblmax = maxval(lblmax_at) 
      write(6,*) "  Atom-specific lmbmax and lblmax:"
      write(6,'(4x,a4,a8,a8)') "iat","lmbmax","lblmax"
      do iat=1,nat
        write(6,'(4x,i4,i8,i8)') iat,lmbmax_at(iat),lblmax_at(iat)
      enddo 
      write(6,*) "lmbmax= ",lmbmax
      write(6,*) "lblmax= ",lblmax 

!
!     calculate Gaunt coefficients
!
      maxnt=max(lmax+1, 2*(lmbmax+1))
      call calcgcoef(maxnt)

!
!     Initial estimation of the maximum number of mixed basis functions       
!
      nmixmax=(lmax+1+(lomax+1)*nlomax)*(lmbmax+1)*(lmbmax+1)*nspin_mb
!
!     allocate the arrays to store the radial part of the mixed basis functions
!
      if(lprt) write(6,*) " init_mixfun"          
      call init_mixfun(nrad,nmixmax,nat)        

      ! calculate the radial part of the mixed basis functions
      do iat = 1, nat
        call mb_setumix(iat)
      enddo 
      maxbigl = maxval(mbl)      
      write(6,*) "maxbigl=",maxbigl

!
!     Calculate the total number of mixed wave functions (including M)    
!     = size of the local part of the matrices
!
      locmatsiz=0
      lmixmax=0
      do iat=1,nat
        lms=0
        do irm=1,nmix(iat)
          lms=lms+(2*bigl(irm,iat)+1)
        enddo  
        if(lms.gt.lmixmax) lmixmax=lms
        locmatsiz=locmatsiz+lms*mult(iat)
      enddo
      write(6,101) lmixmax,locmatsiz  

      nmixmax = maxval(nmix)

      !!reallocate the radial mixed functions
      call reinit_mixfun(nrad,nmixmax,nat)      

!
! set an array that stores the general index of the mixed
! function for a given mixed function of a given atom
!
      allocate(locmixind(lmixmax,ndf),nmixlm(ndf))
      locmixind(:,:)=0
      nmixlm=0
      idf=0
      imix=0
      do iat=1,nat
        do ieq=1,mult(iat)
          idf=idf+1
          im=0
          do irm = 1, nmix(iat)
            l = bigl(irm,iat)
            do m=-l,l
              im=im+1
              imix=imix+1
              locmixind(im,idf)=imix
            enddo
          enddo
          nmixlm(idf)=im
        enddo
      enddo

!----------------------------------------------------------------------!
!      Calculated all matrices related to radial mixbasis functions    !
!----------------------------------------------------------------------!

!
!  Some general matrices (independent of eigen-vectors ) related to
!  bare Coulomb interaction:
!    - rtl
!    - rrint
!    - tilg
!

      if(lprt) write(6,*) "set_mixbasis: calc mixmat"
      call init_mixmat(nat)
      do iat = 1, nat
        call mb_calcrlamint(iat)
      enddo

!
!     Calculate the matrix elements <NL,lambda|lambda'>:
!
      call init_s3r()
      do isp=1,nspin
        do iat = 1, nat
          call mb_calcs3r(iat,isp)
        enddo ! iat
      enddo ! isp

!
!     Calculate the rotation matrices for the spherical harmonics
!
      call gendjmm

  101 format(' Max. nr. of MT-sphere wavefunctions per atom ',i6,/,      &
     &       ' Total  nr. of MT-sphere wavefunctions        ',i6)

      return

      end subroutine set_mixbasis
!EOC
