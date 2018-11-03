!BOP
!
! !ROUTINE: evec2rspmt
!
! !INTERFACE:
      subroutine evec2rspmt(iat,idf,ie,shy,ij)
      
! !DESCRIPTION:      
!
! This subroutine calculates the real space representation of one
! eigenverctor inside the muffin tin sphere of one atom by using equation
!
!\begin{equation}\label{expevecat}
!\Psi^{a}_{n\vec{k}}(\vec{r})=\sum\limits_{l=0}^{lmax}{%
!\sum\limits_{m=-l}^{l}{\left[\mathcal{A}_{lm}^{na}(\vec{k})u_{al}(r)+%
!\mathcal{B}_{lm}^{na}(\vec{k})\dot{u}_{al}(r)+\mathcal{C}_{lm}^{na}(\vec{k})u_{al}(r,E_2)%
!\right]Y_lm(\hat{r})}}
!\end{equation}
!
! where the coefficients
!$\mathcal{A}_{lm}^{na}$, $\mathcal{B}_{lm}^{na}$ and
!$\mathcal{C}_{lm}^{na}$ are defined as in equation \ref{almevec}; and
!$r=|\vec{r}-\vec{r}_a|$ and  $ \widehat{r}=\widehat{\vec{r}-\vec{r}_a}$
!
! !USES:
  
      use constants,   only: czero
      use eigenvec,    only: alfa,beta,gama,alfap,betap,gamap
      use lapwlo,      only: lmax,loor,lomax,lmax
      use radwf,       only: u,udot,ulo,us,usdot,uslo
      use rspevec,     only: evecmt,evecmts,rg,rd
      use struk,       only: nrpt
      
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: iat    ! index of the inequivalent atom
      integer(4), intent(in) :: idf    ! index of the atom (including
!                                        equivalent ones).
      integer(4), intent(in) :: ie     ! band index of the eigenvector
      complex(8), intent(in) :: shy(*) ! Values of the spherical
!                                         harmonics for a given direction
      integer(4), intent(in) :: ij

      
! !LOCAL VARIABLES:

      integer(4) :: i,l,ml,llm       ! counters
      
      complex(8), dimension(0:lmax) :: aa, bb, cc 
!
 
! !INTRINSIC ROUTINES: 


      intrinsic conjg

!
!
!EOP
!BOC
!
!
!     calculate the eigenfunction
!      
      llm=0
      evecmt(1:nrpt(iat))=czero
      evecmts(1:nrpt(iat))=czero
      do l=0,lmax
        aa(l)=czero
        bb(l)=czero
        cc(l)=czero
        select case (ij)
        case (1)
          do ml=-l,l
            llm=llm+1
            aa(l)=aa(l)+alfa(ie,llm,idf)*shy(llm)
            bb(l)=bb(l)+beta(ie,llm,idf)*shy(llm)
            cc(l)=cc(l)+gama(ie,llm,idf)*shy(llm)
          enddo
        case (2)
          do ml=-l,l
            llm=llm+1
            aa(l)=aa(l)+alfap(ie,llm,idf)*conjg(shy(llm))
            bb(l)=bb(l)+betap(ie,llm,idf)*conjg(shy(llm))
            cc(l)=cc(l)+gamap(ie,llm,idf)*conjg(shy(llm))
          enddo
        end select
!
!       apw's and lapw's
!        
        do i=1,nrpt(iat)  
          evecmt(i)=evecmt(i)+aa(l)*u(i,l,iat)+bb(l)*udot(i,l,iat)
          evecmts(i)=evecmts(i)+aa(l)*us(i,l,iat)+bb(l)*usdot(i,l,iat)
        enddo
!        
!       local orbitals
!        
        if(l.le.lomax)then
          if(loor(l,iat))then
            do i=1,nrpt(iat)  
              evecmt(i)=evecmt(i)+cc(l)*ulo(i,l,iat)
              evecmts(i)=evecmts(i)+cc(l)*uslo(i,l,iat)
            enddo
          endif
        endif    
      enddo
!
!     divide by the distance to the nucleus
!      
      do i=1,nrpt(iat)
        evecmt(i)=evecmt(i)/rg(i)
        evecmts(i)=evecmts(i)/rg(i)
      enddo
      
      end subroutine evec2rspmt
!EOC      
