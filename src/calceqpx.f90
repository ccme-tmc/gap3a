!
! Calculate exchange-only self-energies 
!
      subroutine calceqp_ex(eqpx,iop) 

      use bands,      only: bande0,nspin,ibgw,nbgw
      use kpoints,    only: nirkp 
      use xcpot,      only: vxcnn
      use selfenergy, only: sigx,selfx,sigm,vhmn,vh0mn,iop_dvh,vxcmn, &
     &                      znorm,beta_sc,sxcmn,   &
     &                      dvmn_old
      implicit none
      real(8),intent(out):: eqpx(ibgw:nbgw,nirkp,nspin) 
      integer,intent(in):: iop   !! 0/1 - diagonal-only/full Sx matrix  

      integer:: isp,irk,ie

      real(8):: delta
      real(8),allocatable:: rwork(:),eig(:),eig0(:)
      complex(8),allocatable:: work(:)
      complex(8),allocatable:: hmat(:,:),dvmn(:,:)

      character(len=10):: sname="calceqp_ex"


      call linmsg(6,'-','calceqp:exchange-only')
      if(iop.eq.0) then 
        do isp=1,nspin
          do irk=1,nirkp
            do ie=ibgw,nbgw
              delta=real(sigx(ie,irk,isp) - vxcnn(ie,irk,isp)
              eqpx(ie,irk,isp) = bande0(ie,irk,isp) + delta
            enddo ! ie
          enddo ! irk
        enddo ! isp

      else 

        n=nbgw-ibgw+1   !! the size of matrix is n*n
        lwork=2*n
        allocate(hmat(ibgw:nbgw,ibgw:nbgw), &
     &           eig(ibgw:nbgw),            &
     &           eig0(ibgw:nbgw),            &
     &           dvmn(ibgw:nbgw,ibgw:nbgw),  &
     &           work(2*n),                 &
     &           rwork(3*n),stat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to allocat hmat etc.")

        do isp=1,nspin
          do irk=1,nirkp
            
            !! set up the Hamiltonian matrix 
            do ie2=ibgw,nbgw
              do ie1=ibgw,ie2-1
                hmat(ie1,ie2)=0.5d0*( sigm(ie1,ie2,0,irk,isp) &
     &                          +conjg(sigm(ie2,ie1,0,irk,isp)) )
                hmat(ie2,ie1)=conjg(hmat(ie1,ie2))
              enddo
              hmat(ie2,ie2)=bande0(ie2,irk,isp)+real(sigm(ie2,ie2,0,irk,isp))
            enddo
            



      end subroutine

