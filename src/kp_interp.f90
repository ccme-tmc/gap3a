!BOP
! 
! !ROUTINE: kp_interp

!
! !INTERFACE:
      subroutine kp_interp(eband_kp,klist_kp,nk_kp,iopout)  

! !DESCRIPTION:

! This subroutine calculate gw qp energies at a dense mesh using kp-interpolation 
 

! !USES:

      use bands,     only: bande,nspin,nbmax
      use kpoints,   only: kirvecs,nirkp
      use mommat,    only: mmatvv,init_mommat,end_mommat
      use struk,     only: br2
      use task,      only: casename 
!
! !LOCAL VARIABLES:

      implicit none
      integer, intent(in) :: nk_kp                    !! the number of k-points to interpolated 
      real(8), intent(in) :: klist_kp(3,nk_kp)        !! k-vectors to be interpolated
      real(8), intent(out):: eband_kp(nbmax,nk_kp,nspin)  !! band energies to be interpolated 
      integer, intent(in) :: iopout                       !! the option to control output 
                                                          !!  0 -- minimal output 
                                                          !!  1 -- write out coefficients  
                                                          !!  2 -- write out interpolated vectors  

      integer(4) :: ib,jb,nst,nend,nn
      integer(4) :: ik_kp,irk,ik0
      integer(4) :: isp
      integer(4) :: ierr 
      integer(4) :: ig(3),k0ind(4)
      integer :: lwork
      
      real(8) :: k0vec(3),kvec(3),gvec(3),k(3),k0(3),deltak(3) 
      
      complex(8), allocatable :: hmat(:,:),eig(:),mmat(:,:,:,:)
      complex(8), allocatable :: work(:),rwork(:)
      
      character(13) :: sname = 'kp_interp'
      character(80) :: fname 
      character(2)  :: spflag
   
! !REVISION HISTORY:
! Created on March 06, 2010 by H. Jiang

!EOP

!BOC

! The index for the highest band to be interpolated 
      nst=1 
      nend=nbmax 
      nn=nend-nst+1
      write(6,'(a,2i5)') ' nst,nend=',nst,nend 
      allocate(hmat(nst:nend,nst:nend),eig(nst:nend),    &
     &         work(2*nn),rwork(3*nn), stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate work arrays") 

      allocate(mmat(3,nst:nend,nst:nend,nirkp))
      call init_mommat(nst,nend,nst,nend,nirkp,nspin)
      call calcmommat(0,nst,nend,nst,nend,1)

      do isp=1,nspin 
        do ik_kp=1,nk_kp

          k=klist_kp(:,ik_kp)
          call kfrac2cart(k,kvec)
          call kp_k0index(k,kirvecs,nirkp,k0ind)
          ik0=k0ind(1)
          k0=kirvecs(:,ik0) 
          call kfrac2cart(k0,k0vec)

          ig(1:3)=k0ind(2:4)
          gvec(1:3)=dble(ig(1))*br2(1:3,1)+dble(ig(2))*br2(1:3,2)+      &
     &                    dble(ig(3))*br2(1:3,3)
          k0vec(1:3)=k0vec(1:3)+gvec(1:3)
          deltak=kvec-k0vec

          do jb=nst,nend
            hmat(jb,jb)=bande(jb,ik0,isp)+sum(deltak**2)/2.d0
            do ib=nst,jb-1
              hmat(ib,jb)=sum(mmatvv(1:3,ib,jb,ik0,isp)*deltak)
            enddo
          enddo    
          call zheev('v','l',nn,hmat,nn,eig,work,nn*2,rwork,ierr)
          call errmsg0(ierr,sname,"Fail to call zheev")
          eband_kp(nst:nend,ik_kp,isp)=eig
        enddo !! ik_kp

      enddo  !! isp
      deallocate(hmat,eig,work,rwork)

      return

      end subroutine kp_interp

!EOC
