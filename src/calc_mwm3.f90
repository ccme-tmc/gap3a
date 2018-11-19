!BOP
!
! !ROUTINE: calc_mwm3
!
! !INTERFACE: 
      subroutine calc_mwm3(iq,iom_f,iom_l,isc,isym)
      
! !DESCRIPTION:
!
! This subroutine calculate the three-index matrix elements of 'M*W*M',
! i.e. X_{nn';m} in the GAPNOTES
!

! !USES:

      use bands,      only: nbmaxsc,nspin,ibgw,nbgw,nbandsgw,nbands_c
      use bzinteg,    only: singc1co,singc2co
      use constants,  only: cone,czero,pi,fourpi,sqrt4pi
      use core,       only: eigcore,ncg,iop_core 
      use dielmat,    only: eps,head,epsw1,epsw2
      use kpoints,    only: nirkp,nkp,nqp,kqid,kpirind,weightq,get_kvec
      use minmmat,    only: mblksiz
      use mixbasis,   only: matsiz
      use selfenergy, only: qpwf_coef 
      use struk,      only: vi 
      use task,       only: time_lapack,time_selfc,     &
     &                      casename,spflag,savdir,fid_outgw
      use modmpi

! !INPUT PARAMETERS:

      implicit none
      
      integer, intent(in) :: iq      ! index for q-point 
      integer, intent(in) :: iom_f, iom_l ! the first and last index for the frequency points considered in this work 
      integer, intent(in) :: isc     ! control which stage during the self-consistency  
                                     ! -1 -- one-shot calculation, without saving Minm
                                     !  0 -- the first run of selfconsistent GW calculations
                                     !       all states are divided into active and frozen space (AS and FS, respectively)
                                     !       FS contributions are calculated and stored in sigmc_frz, 
                                     !       minm for AS are written to the minm file  
                                     ! >0 -- only AS contributions are re-calculated, Minm are read from the file 
                                     !       and transformed by qpwf_coef 
      integer, intent(in):: isym     !  0/1 -- calculate sigx in the full/irreducible BZ
                                     
! !REVISION HISTORY:
!
! Created Oct. 18, 2017 by H.Jiang
!
!EOP

!BOC
! !LOCAL VARIABLES:            

      integer:: im      ! (Counter) Runs over bands.
      integer:: ik,irk  ! (Counter) Runs over k-points
      integer:: iom     ! (Counter) Runs over frequencies.
      integer:: jk,jrk  ! (Counter) Runs over k-points
      integer:: iik,nktot
      integer:: ierr
      integer:: isp     ! Index for spin 
      integer:: iop_minm 
      integer:: mst, mend 
      integer:: fout=999

      logical:: ldbg = .true.
      real(8):: time1,time2,tstart,tend
      real(8):: kvec(3) 
      
      character(len=20)  :: sname='calc_mwm3'
!
!     Auxiliary arrays 
!
      complex(8), allocatable :: minm(:,:,:)           ! 
      complex(8), allocatable :: mwm(:,:,:,:)

! !EXTERNAL ROUTINES: 
      complex(8), external :: zdotc,zdotu
      external zhemm

      call cpu_time(tstart)

      write(fid_outgw,*) "Calculate M*W*M by ", sname
      if(isc.eq.-1) then 
        iop_minm = -1 
      else if(isc.eq.0) then
        iop_minm = 0
      else
        iop_minm = 1
      endif

      if(isym.eq.0) then
        nktot = nkp
      else
        nktot = nirkp
      endif

      allocate(mwm(ibgw:nbgw,ibgw:nbgw,nbands_c,iom_f:iom_l),stat=ierr)
      call errmsg(ierr.ne.0,sname,"fail to allocate mwm")

      do isp=1,nspin 
        do iik=1, nktot
          ! set ik, irk
          call get_kvec(isym,iik,ik,irk)
          
          do im=1,nbands_c,mblksiz
             mst = im
             mend = min(nbands_c,mst+mblksiz) 
             call sub_setminm(mst,mend) 
     
             do iom = iom_f, iom_l 
               call sub_calcmwm(mwm(ibgw:nbgw,ibgw:nbgw,mst:mend,iom), &
     &                   mst,mend,iom)
             enddo ! iom
       
             deallocate(minm)
           enddo 
           call io_mwm3('w',mwm,1,nbands_c,iom_f,iom_l,isp,iq,irk,ierr)
        enddo ! irk
      enddo ! isp

      call cpu_time(tend)
      time_selfc = time_selfc + tend - tstart
      return

 100  format(2I6,3F10.5,4x,"# irk, ik, k(1:3)") 
 101  format(4g24.16) 
 102  format(' ') 
      contains 
!
!       Calculate the contributions to sigm for m falling the mst..mend block 
!
        subroutine sub_setminm(mst,mend) 
        implicit none 
        integer,intent(in):: mst, mend
        integer:: iom,ie1,imu,inu,cst,cend

        allocate(minm(matsiz,mst:mend,ibgw:nbgw),stat=ierr)
        call errmsg(ierr.ne.0,sname,"fail to allocate minm")

        if(mend.le.nbmaxsc) then
          call get_minm(iop_minm,'nm',minm,ibgw,nbgw,mst,mend,ik,iq,isp)
        elseif(mst.gt.nbmaxsc) then
          cst=mst-nbmaxsc; cend=mend-nbmaxsc
          call get_minm(iop_minm,'nc',minm,ibgw,nbgw,cst,cend,ik,iq,isp)
        else
          do ie1=ibgw,nbgw 
            call get_minm(iop_minm,'nm',minm(:,mst:nbmaxsc,ie1),ie1,ie1, &
     &            mst, nbmaxsc, ik,iq,isp)

            cst=1; cend = mend-nbmaxsc 
            call get_minm(iop_minm,'nc',minm(:,nbmaxsc+1:mend,ie1),ie1,ie1,&
     &            cst, cend,ik,iq,isp)
          enddo
        endif
        end subroutine 

!
!     This subroutine is used as a generic interface to calculate M*W*M
! 
      subroutine sub_calcmwm(xnm,m0,m1,iom)
      implicit none 
      integer,intent(in)::m0,m1,iom
      complex(8),intent(out)::xnm(ibgw:nbgw,ibgw:nbgw,m0:m1)
      integer::nmdim,imu,inu,ie2
      real(8):: coefs1,coefs2,t1,t2
      real(8):: wkq     ! Weight of the k-q combinatiion
      real(8):: vi4pi,sqvi4pi
      complex(8):: xs
      complex(8), allocatable:: wm(:,:,:) 

      vi4pi=fourpi*vi
      sqvi4pi=sqrt(vi4pi)
      wkq = dble(weightq(iq))/dble(nqp)
      coefs2=singc2co*vi4pi
      coefs1=singc1co*sqvi4pi

      nmdim=nbandsgw*(m1-m0+1)
      if(nmdim.le.0) then 
        write(fid_outgw,*) &
     &      trim(sname)//":WARNING -- nmdim <=0 in sub_calcmwm"
        return 
      endif 
      allocate(wm(matsiz,m0:m1,ibgw:nbgw))
      call cpu_time(t1)
      call zhemm('l','u',matsiz,nmdim,cone,eps(:,:,iom),matsiz,  &
     &           minm,matsiz,czero,wm,matsiz)

      do ie2=m0,m1
        do inu=ibgw,nbgw
          do imu=ibgw,nbgw 
            xnm(imu,inu,ie2)=wkq*zdotc(matsiz,minm(:,ie2,imu),1,    &
     &                         wm(:,ie2,inu),1)
            if(iq.eq.1.and.ie2.ge.ibgw.and.ie2.le.nbgw) then   !! add singular contributions 
              xs = coefs2*head(iom)*qpwf_coef(imu,ie2,irk,isp)      &
     &              *conjg(qpwf_coef(inu,ie2,irk,isp))              &
     &           + coefs1*qpwf_coef(imu,ie2,irk,isp)                &
     &              *zdotu(matsiz,minm(:,ie2,inu),1,epsw2,1)        &
     &           + coefs1*conjg(qpwf_coef(inu,ie2,irk,isp))         &
     &              *zdotc(matsiz,minm(:,ie2,imu),1,epsw1,1)
              xnm(imu,inu,ie2)=xnm(imu,inu,ie2)+xs
            endif  
          enddo
        enddo
      enddo
      deallocate(wm)
      call cpu_time(t2)
      time_lapack=time_lapack+t2-t1

      end subroutine 
              
      end subroutine calc_mwm3
!EOC        
