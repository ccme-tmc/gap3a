!BOP
!
! !ROUTINE: calcselfxmn
!
! !INTERFACE: 
      subroutine calcselfxmn(iq,isc,isym)
      
! !DESCRIPTION:
!
!This subroutine calculate the full exchange self-energy matrix elements
!

! !USES:

      use bands,      only: bande,nspin,nomaxs,numins,ibgw,nbgw,nbandsgw
      use barcoul,    only: barcev,barcevsq
      use bzinteg,    only: singc2,kiw
      use constants,  only: cone,czero,pi,fourpi,sqrt4pi
      use core,       only: corind, eigcore,ncg,iop_core 
      use kpoints,    only: nirkp,nkp,kqid,kpirind,get_kvec
      use mixbasis,   only: matsiz
      use selfenergy, only: sigm,sigm_f,qpwf_coef,qpwf_dmat
      use struk,      only: vi
      use task,       only: time_selfx,casename,savdir

      use modmpi

! !INPUT PARAMETERS:
      implicit none
      
      integer, intent(in) :: iq   ! index for q-point 
      integer, intent(in) :: isc  ! indicate the stage of self-consistency, which in the meanwhile have partial 
                                     ! effects on the treatment of Minm   
                                     !  < 0 -- non-self-consistency attempted, Minm are calculated and discarded  
                                     !    0 -- the first step during the self-consistency, Minm are calculated and saved
                                     !  > 0 -- continue self-consistency, Minm read from scratch space 
      integer, intent(in):: isym     !  0/1 -- calculate sigx in the full/irreducible BZ 
                              

! !REVISION HISTORY:
!
! Created Nov. 05, 2009 by H. Jiang 
!
!EOP

!BOC
! !LOCAL VARIABLES:            


      integer:: fid_minm
      integer:: imu,inu,ie2  ! index for bands.
      integer:: ik,irk         ! index for k-point related to n
      integer:: jk,jrk         ! index for k-point related to m  
      integer:: ierr
      integer:: isp            ! Index for spin 
      integer:: nktot,iik
      integer:: nomx,nbands
      integer:: mcounts
      integer:: iop_minm 

      real(8) ::    time1,time2,tstart,tend
      
      character(len=120) :: fn_minm
      character(len=20)::sname='calcselfxmn'

      complex(8), allocatable:: sx(:,:)     ! local exchange self-energy
      logical :: ldbg=.false.

! !EXTERNAL ROUTINES: 
      complex(8), external :: zdotc,zdotu
      external zhemm
      character(10),external:: int2str

      call cpu_time(tstart)

      allocate(sx(ibgw:nbgw,ibgw:nbgw),stat=ierr)
      call errmsg(ierr.ne.0,sname,"fail to allocate sx")

      if(isc.ge.0) then 
        fid_minm=300
        fn_minm=trim(savdir)//trim(casename)//".minm-sx-q"//trim(int2str(iq))
        if(isc.eq.0) then
          open(fid_minm,file=fn_minm,action='write',form='unformatted',&
     &             iostat=ierr)
          call errmsg(ierr.ne.0,sname,"Fail to open minm file to write")
        else
          open(fid_minm,file=fn_minm,action='read',form='unformatted',&
     &             iostat=ierr)
          call errmsg(ierr.ne.0,sname,"Fail to open minm file to read")
        endif
      endif 

      if(isym.eq.0) then
        nktot = nkp
      else
        nktot = nirkp
      endif
      
      do isp=1,nspin 


        !! set the number of bands to be summed 
        nomx=nomaxs(isp)        ! the index of the highest occupied band
        if(iop_core.le.1) then
          nbands=ncg+nomx
        else
          nbands=nomx
        endif

        do iik=1, nktot 

          ! set ik, irk, jk and jrk 
          call get_kvec(isym,iik,ik,irk)   
          jk=kqid(ik,iq)
          jrk=kpirind(jk)

          !! set up expansion coefficients of KS vectors 
          if(isc.le.0 ) then 
            call readvector(ik,1,isp,0)          
            call expand_evec(ik,1,.true.,isp)    
            call readvector(jk,2,isp,0)          
            call expand_evec(jk,2,.true.,isp)    
          endif
!
!         In the first iteration, first calculate the contributions from the frozen space.
!         During the self-consistency, only the contributions from the active space 
!         need to be calculated 
!
          sx=0.d0 
          if(isc.eq.-1) then  
            call sub_calcselfx_blk(1,nbands,-1,1)
          else 
            if(isc.eq.0) then 
              call sub_calcselfx_blk(1,ibgw-1,-1,0)       !! low lying semicore states 
              call sub_calcselfx_blk(nomx+1,nbands,-1,0)  !! core states
              sigm_f(:,:,0,iik,isp) = sigm_f(:,:,0,iik,isp)+sx(:,:)
              call sub_calcselfx_blk(ibgw,nomx,0,1)
            else  
              call sub_calcselfx_blk(ibgw,nomx,1,1)       !! contributions of updated active states 
              !! add the contributions of frozen  states 
              if(iq.eq.1) sx = sx + sigm_f(:,:,0,iik,isp)
            endif 
          endif 
          sigm(:,:,0,iik,isp) = sigm(:,:,0,iik,isp) + sx(:,:)

        enddo ! irk
      enddo ! isp

      if(isc.ge.0) close(fid_minm)

      call cpu_time(tend)
      time_selfx = time_selfx + tend - tstart 

      return
    
      contains 

        subroutine sub_calcselfx_blk(mst,mend,iop_minm,isxs) 
        integer, intent(in):: mst, mend     !! range of the index for the band to be summed over
        integer, intent(in):: iop_minm      !! control the treatment of Minm
                                            !!  -1  -- calculate Minm directly and discard 
                                            !!   0  -- calculate directly and then write to the file 
                                            !!   1  -- read from the disk file and then transform Minm by qpwf_coef
        integer, intent(in):: isxs          !! 0/1 - without/with the singular term 

        integer:: mm,nn,kk 
        integer:: imu,inu,ie2,icg,iat,ic
        complex(8):: mvm
        real(8):: wt,sxs2
        complex(8), allocatable:: cmat(:,:),mmat(:,:,:)
        complex(8), allocatable:: minm(:,:,:)

        if(mend.lt.mst) return 

        allocate(minm(matsiz,mst:mend,ibgw:nbgw),stat=ierr)
        call errmsg(ierr.ne.0,sname,"fail to allocate minm")

        if(iop_minm.eq.-1) then   ! calculate once and for all 
          call sub_calcminm(mst,mend,minm)

        elseif(iop_minm.eq.0) then  ! calculate and save 
          allocate(mmat(matsiz,ibgw:nbgw,ibgw:nbgw))
          call sub_calcminm(ibgw,nbgw,mmat)
          write(fid_minm) mmat

          do imu=ibgw,nbgw 
            do ie2=mst,mend
              minm(:,ie2,imu) = mmat(:,ie2,imu) 
            enddo 
          enddo 
          deallocate(mmat)

        elseif(iop_minm.eq.1) then  ! read from the file and transform 
          mm=matsiz
          nn=mend-mst+1
          kk=nbgw-ibgw+1
          allocate(cmat(kk,nn),mmat(matsiz,ibgw:nbgw,ibgw:nbgw))

          read(fid_minm) mmat

          cmat=conjg(qpwf_coef(:,mst:mend,jrk,isp))
          
          do imu=ibgw, nbgw 
            call zgemm('n','n',mm,nn,kk,cone,mmat(:,:,imu),mm,&
     &       cmat,kk,czero,minm(:,:,imu),mm)
          enddo 
          deallocate(mmat,cmat)
        endif 

        sxs2= - fourpi*vi*singc2
        do inu=ibgw,nbgw
          do imu=ibgw,nbgw

            do ie2=mst,mend
              if(ie2.le.nomx) then
                wt=kiw(ie2,jrk,isp)
              else
                wt=kiw(1,jrk,isp)
              endif

              mvm=zdotc(matsiz,minm(:,ie2,imu),1,minm(:,ie2,inu),1)
              sx(imu,inu)=sx(imu,inu)-mvm*wt
            enddo

            if(iq.eq.1.and.isxs.eq.1) then
              sx(imu,inu)=sx(imu,inu)+sxs2*qpwf_dmat(imu,inu,irk,isp) 
            endif
          enddo
        enddo
        deallocate(minm)
        end subroutine 

        subroutine sub_calcminm(mst,mend,minm) 
        !* calculate minm and/or minc
        integer,intent(in):: mst, mend
        complex(8),intent(out):: minm(matsiz,mst:mend,ibgw:nbgw)
        integer:: cst, cend, ie1

        if(mend.le.nomx) then
          call calcminm(ik,iq,ibgw,nbgw,mst,mend,isp,minm)
        elseif(mst.gt.nomx) then
          cst=mst-nomx
          cend=mend-nomx
          call calcminc(ik,iq,ibgw,nbgw,cst,cend,isp,minm)
        else
          cst=1
          cend=mend-nomx
          do ie1=ibgw,nbgw 
            call calcminm(ik,iq,ie1,ie1,mst,nomx,isp,&
     &                    minm(:,mst:nomx,ie1))
            call calcminc(ik,iq,ie1,ie1,cst,cend,isp, &
     &                    minm(:,nomx+1:mend,ie1))
          enddo 
        endif

        end subroutine 
              
      end subroutine calcselfxmn
!EOC        
