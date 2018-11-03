!BOP
!
! !ROUTINE: task_ldau
!
! !INTERFACE:
      subroutine task_ldau 
      
! !DESCRIPTION:
!
! This subroutine performs perturbative LDA+U calculation based on LDA eigen-energies and vectors  
!
! !USES:
      use bands,       only: nspin,bande,ibgw,nbgw,eqp
      use constants,   only: cone, czero,hev
      use kpoints,     only: wkir,idvkir, kirlist, nirkp, idikp
      use xcpot,       only: vorbnn,dmorb,vorb,iatorb,natorb,nlorb,lorb
      use struk,       only: nat,mult
      use task,        only: casename
      
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: isp,iat,idf,ieq,ikp,irkp,ia,l,il,m1,m2,ie        
      integer(4) :: fid
      integer(4), dimension(3) :: ikvec
      
      complex(8) :: trvd 

! !INTRINSIC ROUTINES: 


      intrinsic cpu_time      

! !REVISION HISTORY:
!
! Created 16.09.2005 by RGA
!      
!EOP
!BOC     

      allocate(vorbnn(ibgw:nbgw,nirkp,nspin),eqp(ibgw:nbgw,nirkp,nspin))
      vorbnn=0.d0

      do isp=1,nspin
        dmorb=0.d0
        do irkp=1, nirkp                      !* Loop over the k-points:
          ikp=idikp(irkp)
          call readvector(ikp,1,isp,0)          !* Read the eigenvector corresponding to the k-point ikp
          call expand_evec(ikp,1,.true.,isp)  !* Calculate the expansion coeficients of the eigenvectors

          idf=0
          do iat = 1, nat                       !*  Loop over inequivalent atoms:
            do ieq = 1, mult(iat)               !* Loop over equivalent atoms:
              idf = idf + 1
              call calcvorbnn(iat,idf,irkp,isp)
            enddo 
          enddo 
        enddo 
        do ia=1,natorb
          iat=iatorb(ia)
          do il=1,nlorb(ia)
            l=lorb(il,ia)

            write(6,*)
            write(6,*) '--- dm for l=',l,'iat=',iat
            write(6,*)
            trvd=czero
            do m2=-l,l
              do m1=-l,l
                trvd = trvd+vorb(m2,m1,il,ia,isp)*dmorb(m1,m2,il,ia)
                write(6,'(2e16.8)') dmorb(m1,m2,il,ia)
              enddo
            enddo
            write(6,*)
            write(6,'(a,2e16.8)') "Tr(rho.V)=",trvd
          enddo
        enddo
        write(6,*)
      enddo 

      write(6,*)'-----------------------------------------------------'
      write(6,*)'             Perturbative LDA+U  Output'
      write(6,*)'------------------------------------------------------'

!
!     Loop over spin
!

      fid=999
      open(unit=fid,  file=trim(casename)//".eqpH_ldau")

      do isp=1,nspin
        do irkp=1,nirkp
          ikvec(1:3)=kirlist(1:3,irkp)
          write(fid,1) ikvec,idvkir,irkp,nirkp,nbgw-ibgw+1
          write(6,2)
          do ie=ibgw,nbgw
            eqp(ie,irkp,isp)=bande(ie,irkp,isp)+vorbnn(ie,irkp,isp)
            write(fid,3) ie,bande(ie,irkp,isp),eqp(ie,irkp,isp),0.d0,0.d0,0.d0,0.d0 
            write(6,4) irkp,ie,bande(ie,irkp,isp)*hev,eqp(ie,irkp,isp)*hev
          enddo ! ie
          write(fid,*) 
          write(6,*) 
        enddo ! irkp
      enddo ! isp
      close(fid)

      return
    1 format(7i6)
    2 format(2x,'ikp ie',7x,'E_lda',9x,'E_lda+u',11x)
    3 format(i4,6f20.15)
    4 format(2i4,2f14.6)
      return
      
      end subroutine task_ldau
!EOC      
