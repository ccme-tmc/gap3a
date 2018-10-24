!BOP
!
! !ROUTINE: calceffmass
!
! !INTERFACE:

      subroutine calceffmass(mass,n,irk,isp,lcore,iop_tran)
! !DESCRIPTION:
! This subroutine calculates the effective mass at given n and k-points 
!
! !USES:
      use kpoints,     only: kpirind,idikp
      use bands,       only: bande,nbmax
      use core,        only: ncg,eigcore,corind
      use mommat,      only: lhermit 
      implicit none 
      
! !INPUT PARAMETERS:

!EOP
!BOC
      integer, intent(in) :: n,irk,isp
      integer, intent(in) :: iop_tran  ! control whether to use transformed vectors 
      logical, intent(in) :: lcore     ! control whether include core states contributions  
      real(8), intent(out):: mass(3) 

! !LOCAL VARIABLES:
      logical :: ldbg =.true.
      character(20):: sname="calceffmass" 
      integer :: m
      integer :: i,j, ik,iat,ic 
      integer :: ierr
      real(8) :: en,em,minv(3,3) 
      complex(8):: pnm(3),pmn(3) 
      real(8) :: work(8) 
      
      complex(8), allocatable:: pmat(:,:) 
      

      ik=idikp(irk) 
      en=bande(n,irk,isp) 

      minv=0.0 
      do i=1,3
        minv(i,i) = 1.0 
      enddo 

      call readvector(ik,1,isp,iop_tran)
      call expand_evec(ik,1,.false.,isp)
      
      do m=1,nbmax
        if(m.eq.n) cycle 
        em = bande(m,irk,isp) 
        if(abs(em-en).lt.1.e-5) then 
          write(6,*) "WARNING: near degeneracy occurs" 
          write(6,*) " -- neglect the conbitutions from ",m
          cycle 
        endif 

        call calcmmatvv(n,m,ik,isp,pnm) 
        if(lhermit) then 
          call calcmmatvv(m,n,ik,isp,pmn)
          pnm =  0.5d0*(pnm+conjg(pmn))
        endif 

        do j=1,3
          do i=1,3
            minv(i,j) = minv(i,j) + 2.0*real(pnm(i)*conjg(pnm(j)))    &
     &                 /(en-em)
          enddo 
        enddo 
      enddo 
      if(ldbg) then 
        write(6,*) "1/m* tensor:"
        write(6,'(3f10.4)') ( (minv(i,j),i=1,3),j=1,3) 
      endif 

      if(lcore) then 
        allocate(pmat(3,ncg))
        call calcmmatcv(n,ik,isp,pmat)
        do m=1,ncg
          iat= corind(1,m)
          ic = corind(3,m)
          em = eigcore(ic,iat,isp)
          do j=1,3
            do i=1,3
              minv(i,j)=minv(i,j)+2.0*real(conjg(pmat(i,m))*pmat(j,m)) &
     &                  /(en-em)
            enddo 
          enddo
        enddo 
        if(ldbg) then 
          write(6,*) "1/m* tensor after adding core states:"
          write(6,'(3f10.4)') ( (minv(i,j),i=1,3),j=1,3) 
        endif 
        deallocate(pmat)
      endif
 
      !! diagonalize 1/m* tensor  
      call dsyev('v','l',3,minv,3,mass,work,8,ierr)
      call errmsg(ierr.ne.0,sname,"Fail to call dsyev") 

      do i=1,3
        mass(i)=1.0/mass(i) 
        write(6,100) i,mass(i),minv(1:3,i)
      enddo

 100  format("m(",i,")=",f8.3,"v=(",3f8.3,")") 
    
      end subroutine 
