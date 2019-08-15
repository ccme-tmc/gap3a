!BOP
!
! !ROUTINE: mb_setumix
!
! !INTERFACE:
      subroutine mb_setumix(iat)
      
! !DESCRIPTION:
!
! Set up the radial part of the mixed basis functions ($\upsilon_{aNL}$)used for the matrix
! expansion of the non local operators (Polarization, Bare and Screened Coulomb
! potential, and the Self energy) for atom \verb"iat". The procedure is as follows:
!
!\begin{itemize}
!\item For each $L$ we take the product of radial functions $u_{l}(r)u_{l'}(r)$ which fullfil the
!condition $|l-l'|\le L \le l+l'$.
!
!\item We alculate the overlap matrix of the products of radial functions:
!\begin{equation}
!\mathbb{O}_{(l,l');(l_1,l'_1)}=\int\limits_0^{R^a_{MT}} {u_{al}u_{al'}(\vec{r})u_{al_1}u_{al'_1}(\vec{r})
!r^2 dr}
!\end{equation}
!
!\item Diagonalize the matrix $\mathbb{O}_{(l,l');(l_1,l'_1)}$
!
!\item Discard the eigenvectors corresponding to eigenvalues with absolute value lower than a given
!tolerance (usually $10^{-5}$)
!
!\item The rest of the eigenvectors are normalized and stored for a grid of $\vec{r}$ that constitute the
!new basis $\left\{\upsilon_{NL}\right\}$
!
!\end{itemize}
!
!So defined the set of functions $\left\{\gamma
!_{aNLM}=\upsilon_{aNL}Y_{LM}\right\}$ constitute and orthonormal basis set,
!that is:
!
!\begin{equation}\label{mixborth}
!\int\limits_{V^a_{MT}}\gamma _{aNLM}(\vec{r})\gamma
!_{aN'L'M'}(\vec{r})d^3r=\delta_{N,N'}\delta_{L,L'}\delta_{M,M'}
!\end{equation}
!
!

! !USES:      

      use bands,     only: nspin
      use core,      only: ncore
      use constants, only: cfein1, cfein2
      use lapwlo,    only: lmax
      use mixbasis,  only: lblmax_at,bigl,mbl,&
     &                     nmix,umix,usmix,wftol
      use prodfun,   only: eles, nup, rp, umat, uprod, usprod,           &
     &                     init_prodfun, end_prodfun
      use radwf,     only: a,b
      use struk,     only: atomname,nrpt 
 
! !INPUT PARAMETERS:      
      implicit none
      
      integer(4), intent(in) :: iat

! !LOCAL VARIABLES:      

      character(10)::sname="mb_setumix"
      character(120)::msg
      integer:: fid
      integer(4) :: i,j,k,l,m,ir
      integer(4) :: iup           !! Index for u-product 
      integer(4) :: npt,nli,nor,nll
      integer(4) :: lammax, lammin, nwf
      integer(4) :: i1,i2,info,j2
      
      integer(4), allocatable :: ind(:)
      integer(4), dimension(0:2*lmax) :: nl      
      real(8) :: norm
      real(8), allocatable :: uol(:),work(:)
      real(8), allocatable :: uml(:,:)

      logical::lprt=.false.


!
! !EXTERNAL ROUTINES: 
!

 
      external dsyev
      external radmesh
      external rint13
      external ylm
!
! !REVISION HISTORY:
!
! Created May. 19th. 2004 by RGA
! Last Modified May 25th 2004 by RGA
!
!EOP
!BOC
!  ---------------------------------------------------------------------            
!
      fid=6
      write(fid,100) atomname(iat)
      npt = nrpt(iat)
      nmix(iat)=0

!
!     count the maximum number of product functions 
!
      call mb_setuprod(iat,0)

!
!     allocate the arrays for the product functions and the l,l' indexes      
!
      call init_prodfun(npt)
!
!     calculate the radial mesh points
!
      call radmesh(iat,rp)
!
!     Calculate the radial product functions and their overlap matrix
!
      call mb_setuprod(iat,npt)
!
!     Calculate the number of radial product functions for each L and the 
!     total number of orbitals
!
      nl(:)=0   
      do iup=1,nup
        lammax=eles(1,iup)+eles(2,iup)
        lammin=abs(eles(1,iup)-eles(2,iup))

        do l=0,lblmax_at(iat)       
          if((l.le.lammax).and.(l.ge.lammin))then
            nl(l) = nl(l) + 1
          endif 
        enddo
      enddo

      umix(:,:,iat)=0.0d0
      usmix(:,:,iat)=0.0d0
      bigl(:,iat)=0
!
!     diagonalize the overlap matrix of the product functions for each  L-block
!
      j2=0
      nor=0
      do l=0,lblmax_at(iat)
        nll=nl(l) 

        if(nll.gt.1)then
!
!         allocate the L-block overlap matrix (uml), the eigenvalues
!         vector (uol), the working space and the indexes
!
          allocate(uml(nll,nll))
          allocate(uol(nll))
          allocate(work(3*nll-1))
          allocate(ind(nll))
!
!         generate the L-block overlap matrix

          if(lprt) write(6,*) " - generate the L-block overlap matrix"
          j=0
          do i1=1,nup
            lammax = eles(1,i1)+eles(2,i1)
            lammin = abs(eles(1,i1)-eles(2,i1))
            if((l.le.lammax).and.(l.ge.lammin))then
              j=j+1
              ind(j) = i1
              uml(j,j) = umat(i1,i1)
              k=j
              do i2=i1+1,nup
                lammax = eles(1,i2)+eles(2,i2)
                lammin = abs(eles(1,i2)-eles(2,i2))
                if((l.le.lammax).and.(l.ge.lammin))then
                  k=k+1
                  uml(j,k)=umat(i1,i2)
                  uml(k,j)=uml(j,k) 
                endif
              enddo  !i2  
            endif
          enddo  ! i1

!
!         diagonalize the L-block overlap matrix
!
          call dsyev('v','u',nll,uml,nll,uol,work,3*nll-1,info)
          call errmsg0(info,sname,"calling dsyev")

          ! determine how many eigenvectors whose eigenvalues are larger than wftol 
          i2=0
          do i1=1,nll
            if(uol(i1).gt.wftol) i2 = i2+1
          enddo 
          nli=i2
          

          if(lprt) write(6,*) " - transform the umix"
          do i1=1,nll
            if(uol(i1).le.wftol) cycle 
            j2=j2+1
            do i2=1,nll
              umix( :,j2,iat) = umix(:,j2,iat)+uml(i2,i1)*uprod( :,ind(i2))
              usmix(:,j2,iat) =usmix(:,j2,iat)+uml(i2,i1)*usprod(:,ind(i2))
            enddo 
          enddo
!
!         normalize the radial mixed wave functions
!
          if(lprt) write(6,*) " - normalize the rad mixed func"
          do i1=nor+1,nor+nli
            bigl(i1,iat)=l
            a(1:npt)=umix( 1:npt,i1,iat)
            b(1:npt)=usmix(1:npt,i1,iat)
            call rint13(cfein1,cfein2,a,b,a,b,norm,iat)
            umix( 1:npt,i1,iat)=umix( 1:npt,i1,iat)/sqrt(norm)
            usmix(1:npt,i1,iat)=usmix(1:npt,i1,iat)/sqrt(norm)
          enddo

          write(6,101) l,nll,nli

          nor=nor+nli
          nll=nli
!
!         deallocate the temporary arrays
!
          deallocate(uml)
          deallocate(uol)
          deallocate(work)
          deallocate(ind)
!
!       in case the L-block is just one wavefunctions
!
        elseif(nll.eq.1)then
          if(lprt) write(6,*) " nll .eq. 1"
          do i1=1,nup
            lammax=eles(1,i1)+eles(2,i1)
            lammin=abs(eles(1,i1)-eles(2,i1))
            if((l.le.lammax).and.(l.ge.lammin))then
              j2=j2+1
              bigl(j2,iat)=l
              a(1:npt)=uprod(1:npt,i1)
              b(1:npt)=usprod(1:npt,i1)
              call rint13(cfein1,cfein2,a,b,a,b,norm,iat)
              umix(1:npt,j2,iat)=uprod(1:npt,i1)/sqrt(norm)
              usmix(1:npt,j2,iat)=usprod(1:npt,i1)/sqrt(norm)
            endif  
          enddo
          write(6,101) l,1,1
          nor=nor+1
        endif
      enddo  ! l 
      nmix(iat) = nor
      mbl(iat)  = maxval(bigl(:,iat))

      nwf=0
      do i=1,nor
        nwf=nwf+2*bigl(i,iat)+1
      enddo
      write(fid,102) nmix(iat),mbl(iat)
      write(fid,105) nwf
          
      call end_prodfun

  100 format(/,10x,'Mixed basis functions for atom',1x,a10,/)
  101 format(10x,'L =',i2,' Nr. of products:',i4, &
     &      ' Nr. of basis functions',i4)  
  102 format(59x,'----',/,10x,'Total number of radial functions',17x,   &
     &       i4,4x,'Maximum L',i4)
  103 format(26x,'N',3x,'L',2x,'deg.')
  104 format(23x,3i4)
  105 format(10x,'Total number of basis functions',i4)
      return
      
      end subroutine mb_setumix
!EOC
