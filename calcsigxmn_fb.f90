!BOP
!
! !ROUTINE: calcsigxmn_fb
!
! !INTERFACE: 
      subroutine calcsigxmn_fb(iq,isc)
      
! !DESCRIPTION:
!
! This subroutine calculate the full exchange self-energy matrix elements
! using the fixed basis approach 
!

! !USES:

      use bands,      only: bande,nbmax,nspin,      &
     &                      nomaxs,numins,ibgw,nbgw,nbandsgw
      use barcoul,    only: barcev,barcevsq
      use bzinteg,    only: singc2,kiw,ciw
      use constants,  only: cone,czero,pi,nbcmplx,fourpi,sqrt4pi
      use core,       only: corind, eigcore,ncg,iop_core 
      use kpoints,    only: nirkp,nqp,kqid,idikp,     &
     &                      kpirind
      use minmmat,    only: mblksiz
      use mixbasis,   only: matsiz
      use selfenergy, only: sigmx,sigmx_frz,qpwf_coef
      use struk,      only: vi
      use task,       only: iop_scratch,time_selfx,casename,savdir

      use modmpi

! !INPUT PARAMETERS:
      implicit none
      
      integer(4), intent(in) :: iq   ! index for q-point 
      integer(4), intent(in) :: isc  ! indicate the stage of self-consistency, which in the meanwhile have partial 
                                     ! effects on the treatment of Minm   
                                     !  < 0 -- non-self-consistency attempted, Minm are calculated and discarded  
                                     !    0 -- the first step during the self-consistency, Minm are calculated and saved
                                     !  > 0 -- continue self-consistency, Minm read from scratch space 

! !REVISION HISTORY:
!
! Created Nov. 05, 2009 by H. Jiang 
!
!EOP

!BOC
! !LOCAL VARIABLES:            

      logical ::    ldbg=.false.

      integer:: fid_minm
      integer:: ie1m,ie1n,ie2  ! index for bands.
      integer:: ik,irk         ! index for k-point related to n
      integer:: jk,jrk         ! index for k-point related to m  
      integer:: ierr
      integer:: isp            ! Index for spin 
      integer:: nomx,numn,nbands
      integer :: mcounts

      real(8) ::    time1,time2,tstart,tend
      character(len=120) :: fn_minm
      character(len=10)::sname='calcsigxmn_fb'

      complex(8), allocatable:: sx(:,:)     ! local exchange self-energy
      complex(8),allocatable:: minm(:,:,:)

! !EXTERNAL ROUTINES: 
      complex(8), external :: zdotc,zdotu
      external zhemm
      character(10),external:: int2str

      if(ldbg) call linmsg(6,'-',sname)
      call cpu_time(tstart)

      allocate(sx(ibgw:nbgw,ibgw:nbgw),stat=ierr)
      call errmsg(ierr.ne.0,sname,"fail to allocate sx")

      if(iop_scratch.gt.0.and.isc.ge.0) then 
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
      
      do isp=1,nspin 
        nomx=nomaxs(isp) 
        numn=numins(isp)

        do irk=1, nirkp
          ik=idikp(irk)  
          jk=kqid(ik,iq)
          jrk=kpirind(jk)

          if(iop_core.le.1) then
            nbands=ncg+nomx
          else
            nbands=nomx
          endif

          if(isc.le.0.or.iop_scratch.eq.0) then !! first iteration in the self-consistency
            call readvector(ik,1,isp,0)          !! Read the eigenvector corresponding to the k-point ik
            call expand_evec(ik,1,.true.,isp)  !! Calculate the expansion coeficients of the eigenvectors
            call readvector(jk,2,isp,0)          !! Read the eigenvector corresponding to the k'-point jk
            call expand_evec(jk,2,.true.,isp)  !! Calculate the expansion coeficients of the eigenvectors
          endif
!
!           In the first iteration, first calculate the contributions from the frozen space.
!           During the self-consistency, only the contributions from the active space 
!           need to be calculated 
!
          sx=0.d0 
          if(isc.eq.-1) then 

            call sub_calcselfx_blk(1,nbands,-1)

          elseif(isc.eq.0) then 
            call sub_calcselfx_blk(1,ibgw-1,-1)       !! low lying semicore states 
            call sub_calcselfx_blk(nomx+1,nbands,-1)  !! core states
            sigmx_frz(:,:,irk,isp)=sigmx_frz(:,:,irk,isp)+sx(:,:)
            call sub_calcselfx_blk(ibgw,nomx,0)       !! contributions of active states 
          else  
            call sub_calcselfx_blk(ibgw,nomx,1)       !! contributions of updated active states 
          endif 
          sigmx(:,:,irk,isp)=sigmx(:,:,irk,isp)+sx(:,:)

        enddo ! irk
      enddo ! isp

      if(iop_scratch.gt.0.and.isc.ge.0) close(fid_minm)
      call cpu_time(tend)
      time_selfx = time_selfx + tend-tstart
      return
 101  format(2i5,3f12.4) 
    
      contains 

        subroutine sub_calcselfx_blk(mst,mend,iop_minm) 
        integer, intent(in):: mst, mend     !! range of the index for the band to be summed over
        integer, intent(in):: iop_minm      !! control the treatment of Minm
                                            !!  -1  -- calculate Minm directly and discard 
                                            !!  0  -- calculate directly and then write to the file 
                                            !!  1 -- read from the disk file and then transform Minm by qpwf_coef
        integer:: mm,nn,kk 
        complex(8),allocatable:: cmat(:,:),mmat(:,:,:)

        if(mend.lt.mst) return 

        allocate(minm(matsiz,ibgw:nbgw,mst:mend),stat=ierr)
        call errmsg(ierr.ne.0,sname,"fail to allocate minm")

        if(iop_minm.eq.-1) then 

          call sub_setminm(mst,mend,minm)

        elseif(iop_minm.eq.0) then 

          if(iop_scratch.gt.0) then 
            allocate(mmat(matsiz,ibgw:nbgw,ibgw:nbgw))
            call sub_setminm(ibgw,nbgw,mmat)
            write(fid_minm) mmat
            minm=mmat(:,:,mst:mend) 
            deallocate(mmat)
          else 
            call sub_setminm(mst,mend,minm)
          endif 

        elseif(iop_minm.eq.1) then 
          mm=matsiz*(nbgw-ibgw+1)
          nn=mend-mst+1
          kk=nbgw-ibgw+1
          allocate(cmat(kk,nn),mmat(matsiz,ibgw:nbgw,ibgw:nbgw))

          if(iop_scratch.gt.0) then 
            read(fid_minm) mmat
          else 
            call sub_setminm(ibgw,nbgw,mmat)
          endif 

          cmat=conjg(qpwf_coef(:,mst:mend,jrk,isp))
          call zgemm('n','n',mm,nn,kk,cone,mmat,mm,cmat,kk,czero,minm,mm)
          deallocate(mmat,cmat)

        endif 

        if(ldbg) write(6,*) "setselfx"
        call sub_setselfx(mst,mend)
        deallocate(minm)

        end subroutine 

        subroutine sub_setminm(mst,mend,minm) 
        !* calculate minm and/or minc
        integer,intent(in):: mst, mend
        complex(8),intent(out):: minm(matsiz,ibgw:nbgw,mst:mend)
        integer:: cst, cend 

        if(ldbg) write(6,*) "calcminm"
        if(mend.le.nomx) then
          call calcminm(ik,iq,ibgw,nbgw,mst,mend,isp,minm(:,:,mst:mend))
        elseif(mst.gt.nomx) then
          cst=mst-nomx
          cend=mend-nomx
          call calcminc(ik,iq,ibgw,nbgw,cst,cend,isp,minm(:,:,mst:mend))
        else
          call calcminm(ik,iq,ibgw,nbgw,mst,nomx,isp,minm(:,:,mst:nomx))
          cst=1
          cend=mend-nomx
          call calcminc(ik,iq,ibgw,nbgw,cst,cend,isp,minm(:,:,nomx+1:mend))
        endif

        end subroutine 

        subroutine sub_setselfx(mst,mend)
        integer(4),intent(in)::mst,mend
        integer::ie1m,ie1n,ie2,icg,iat,ic
        complex(8) :: mvm
        real(8) :: wt,sxs2

        sxs2= - fourpi*vi
        do ie2=mst,mend
          if(ie2.le.nomx) then 
            wt=kiw(ie2,jrk,isp)
          else 
            icg=ie2-nomx
            iat=corind(1,icg)
            ic=corind(3,icg)
            wt=ciw(iat,ic,isp)
          endif 
          do ie1n=ibgw,nbgw
            do ie1m=ibgw,nbgw
              mvm=zdotc(matsiz,minm(:,ie1m,ie2),1,minm(:,ie1n,ie2),1)
              sx(ie1m,ie1n)=sx(ie1m,ie1n)-mvm*wt
              if(iq.eq.1.and.ie1m.eq.ie1n.and.ie1m.eq.ie2.and.ie1m.le.nomx) then 
                sx(ie1m,ie1n)=sx(ie1m,ie1n)+sxs2*singc2
              endif 
            enddo 
          enddo
        enddo
        end subroutine ! internal subroutine setselfx
              
      end subroutine calcsigxmn_fb
!EOC        
