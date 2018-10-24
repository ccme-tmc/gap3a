!BOP
!
! !ROUTINE: calcmommat
!
! !INTERFACE:
      subroutine calcmommat(iop_tran,nst,nend,mst,mend,iop_core)
      
! !DESCRIPTION:
!
! This subroutine calculates the momentum matrix elements 
!
! !USES:
      
      use bands,       only: bande,bande0,nspin
      use kpoints,     only: kpirind,idikp,nirkp
      use mommat,      only: mmatvv,mmatcv,lhermit,lrenorm
      use task,        only: fid_outmom
      
! !INPUT PARAMETERS:

!EOP  
!BOC  
      implicit none
      integer, intent(in) :: iop_tran            !! control whether (>0) or not (=0) to use transformed wave functions 
      integer, intent(in) :: nst,nend,mst,mend   !! the range of 1st band index 
      integer, intent(in) :: iop_core             !!  0/1 include/exclude moment matrix involving core states 
         
! !LOCAL VARIABLES:

      integer :: ik,irk  ! index of the k point.
      integer :: ie1,ie2 ! Counter: run over bands.
      integer :: isp     ! index for spin 
      integer :: nbk     ! k-dependent number of bands
      integer :: fout 

      complex(8) :: p12(3),p21(3)
      real(8) :: renorm
      logical :: ldbg=.true.
      character(20):: sname="calcmommat"

! !EXTERNAL ROUTINES: 
      external expand_evec
      external readvector
      external calcmmatvv
      external calcmmatcv
 
      if(ldbg) call linmsg(6,'-',trim(sname))  
      lhermit=.true.

      fout = fid_outmom

      do isp=1,nspin 
        do irk=1,nirkp
          !! Generate the expansion coefficient of the eigenvectors 
          ik=idikp(irk)
          call readvector(ik,1,isp,iop_tran)
          call expand_evec(ik,1,.false.,isp)
          do ie2=mst,mend
            do ie1=nst,nend

              if(lrenorm) then 
                renorm= ( bande(ie2,irk,isp) -  bande(ie1,irk,isp) )    &
     &              /( bande0(ie2,irk,isp) - bande0(ie1,irk,isp) )
              else 
                renorm = 1.d0
              endif 

              call calcmmatvv(ie1,ie2,ik,isp,p12)
              if(lhermit) then 
                call calcmmatvv(ie2,ie1,ik,isp,p21)    !! contribution from MT spheres
                p12=0.5d0*(p12+conjg(p21))
              endif 
              mmatvv(1:3,ie1,ie2,irk,isp) = p12(1:3)*renorm
            enddo ! ie1

            if(iop_core.eq.0) then 
              call calcmmatcv(ie2,ik,isp,mmatcv(:,:,ie2,irk,isp)) 
            endif 
          enddo ! ie2

      !! write out mommat 
          write(fout,*) "#momentum matrix for ik=",ik,"irk=",irk
          write(fout,'(a1,a5,a6,3a12)')'#',"ie1","ie2","px","py","pz"
          do ie1=nst,nend
            do ie2=mst,mend
              p12 = mmatvv(:,ie1,ie2,irk,isp)
              write(fout,'(2i6,4f12.6)')ie1,ie2,real(p12), &
     &          real(sum(p12*conjg(p12))/3.d0)
            enddo ! ie2
          enddo 
        enddo ! irk
      enddo ! isp
      end subroutine calcmommat
!EOC
            
