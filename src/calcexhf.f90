!BOP
!
! !ROUTINE: calcexhf
!
! !INTERFACE:
      subroutine calcexhf(iq)

!
! !DESCRIPTION:
!
!This subroutine calculates the q-point contribution to Hartree-Fock exchange energy 
!

! !USES:
      use acfd,        only: ex_hf
      use bands,       only: nomaxs,nspin,fspin
      use bzinteg,     only: kiw,singc2ex,kwt_ibz,kwt_bz
      use constants,   only: cone,czero,twopi
      use core,        only: ncg,nclm,corind,nclmmax
      use kpoints,     only: nkp,wkir,kpirind,nqp,idikp,kqid
      use mixbasis,    only: matsiz,lmixmax,locmixind,nmixlm,nmix,bigl, &
     &                       locmatsiz
      use struk,       only: ndf,inddf,neltot,vi
     
      
      implicit none 
! ! !INPUT PARAMETERS:
      integer(4),intent(in) :: iq



! !LOCAL VARIABLES:
      integer(4) :: ik,jk,irk,jrk  ! indices for k-vectors in BZ

      integer(4) :: imix,im,jm,lmix,irm ! counter for mixed basis function 
      integer(4) :: im1,im2,m1

      integer(4) :: ie1,ie2,ie12 ,ie12max! counter for bands 
      integer(4) :: icg,idf,iat,ic1,ic2,iclm,jclm 
      integer(4) :: ierr
      integer(4) :: isp  ! index for spin 
      integer(4) :: nmdim,ccdim,nomx 

      real(8)::exq=0.0D0,exq0s=0.0D0
      real(8) :: kwt
      real(8):: time1,time2
      complex(8):: mm 
 
      complex(8),allocatable::minm(:,:) 

! ! external functions
      complex(8), external:: zdotc 

!! Created by Hong Jiang on Sept. 20, 2009

!
! Setup the mm mattrix 
!
      
      nmdim = max(maxval(nomaxs),ncg)**2 

      allocate(minm(matsiz,nmdim))

      if(iq.eq.1) then
        exq0s= -singc2ex*twopi*neltot*vi  !! singular term of the HF exchange energy
        write(6,'(a,f10.4)') "  -- singular (q=0) term contr.=",exq0s
        ex_hf=ex_hf+exq0s
      endif

      exq=0.d0 
      do isp=1,nspin 
        nomx=nomaxs(isp)

        do ik=1,nkp
          jk=kqid(ik,iq)
          irk=kpirind(ik)
          jrk=kpirind(jk)

          call readvector(ik,1,isp,0)          !! Read the eigenvector corresponding to the k-point ik
          call expand_evec(ik,1,.true.,isp)  !! Calculate the expansion coeficients of the eigenvectors
          call readvector(jk,2,isp,0)          !! Read the eigenvector corresponding to the k'-point jk
          call expand_evec(jk,2,.true.,isp)  !! Calculate the expansion coeficients of the eigenvectors

          !! m,n both are band states 
          call calcminm(ik,iq,1,nomx,1,nomx,isp,minm)
          ie12=0
          do ie2=1,nomx
            do ie1=1,nomx
              ie12=ie12+1
!              write(6,*) "minm=",minm(1:matsiz:10,ie12)

              if(kiw(ie1,irk,isp)*kiw(ie2,jrk,isp).lt.1.d-20) cycle 
              mm= zdotc(matsiz,minm(:,ie12),1,minm(:,ie12),1)
!              write(6,101) ie1,ie2,kiw(ie1,ik,isp),kiw(ie2,jk,isp),mm
              exq=exq - 0.5d0*fspin*kiw(ie1,irk,isp)*real(mm)
            enddo 
          enddo 
 101  format(2i5,4f10.4)
       
          !! n => band m => core 
          if(ncg.eq.0) cycle 

          ccdim=sum(nclm(1:ndf)**2)*ndf 
          call calcminc(ik,iq,1,nomx,1,ncg,isp,minm) 
          ie12=0
          do icg=1,ncg
            do ie1=1,nomx
              ie12=ie12+1
              if(kiw(ie1,irk,isp).lt.1.d-10) cycle
              mm=zdotc(matsiz,minm(:,ie12),1,minm(:,ie12),1)
              exq=exq - 0.5d0*fspin*kiw(ie1,irk,isp)*real(mm)
            enddo
          enddo

          ! n => core, m => band
          call calcmicm(ik,iq,1,ncg,1,nomx,isp,minm)
          ie12=0
          do ie2=1,nomx
            kwt=kiw(ie2,jrk,isp)
            do icg=1,ncg
              ie12=ie12+1
              if(kiw(ie2,jrk,isp).lt.1.d-10) cycle
              mm=zdotc(matsiz,minm(:,ie12),1,minm(:,ie12),1)
              exq=exq - 0.5d0*fspin*kwt*real(mm)
            enddo
          enddo

          !! n => core, m => core
          call calcmicc(isp,ccdim,minm)
          kwt = kiw(1,1,1)
          ie12=0
          do idf=1,ndf
            do jclm=1,nclm(idf)
              do iclm=1,nclm(idf)
                ie12 = ie12+1
                mm = zdotc(matsiz,minm(:,ie12),1,minm(:,ie12),1)
                exq = exq - 0.5d0*fspin*kwt*real(mm)
              enddo 
            enddo
          enddo
        enddo !! ikp 
      enddo !! isp 

      deallocate(minm)
      ex_hf=ex_hf+kwt_bz(iq)*exq

      return 

      end subroutine 
