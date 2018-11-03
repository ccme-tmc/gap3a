!BOP
!
! !ROUTINE: calcselfcmn
!
! !INTERFACE: 
      subroutine calcselfcmn(iq,isxc,isc,isym)
      
! !DESCRIPTION:
!
! This subroutine calculate full correlation self-energy matrix using
! the initial Kohn-Sham vectors as the basis  
!

! !USES:

      use bands,      only: bande,nbmax,nspin,nomaxs,numins,ibgw,nbgw,  &
     &                      nbandsgw
      use bzinteg,    only: singc1,singc2
      use constants,  only: cone,czero,pi,fourpi,sqrt4pi
      use core,       only: corind, eigcore,ncg,iop_core 
      use dielmat,    only: eps,head,epsw1,epsw2
      use freq,       only: nomeg,omega,womeg
      use kpoints,    only: nirkp,nkp,nqp,kqid,kpirind,weightq,get_kvec
      use mixbasis,   only: matsiz
      use selfenergy, only: sigm,sigm_f,qpwf_coef 
      use struk,      only: vi 
      use task,       only: time_lapack,time_selfc,     &
     &                      casename,spflag,savdir
      use modmpi

! !INPUT PARAMETERS:

      implicit none
      
      integer, intent(in) :: iq      ! index for q-point 
      integer, intent(in) :: isxc    ! which approximation to the selfenergy
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
! Created Nov. 05, 2009 by H.Jiang
!
!EOP

!BOC
! !LOCAL VARIABLES:            

      integer:: ie1     ! (Counter) Runs over bands.
      integer:: ik,irk  ! (Counter) Runs over k-points
      integer:: iom     ! (Counter) Runs over frequencies.
      integer:: jk,jrk  ! (Counter) Runs over k-points
      integer:: iik,nktot
      integer:: nhomo,nbands
      integer:: ierr
      integer:: isp     ! Index for spin 
      integer:: iop_minm 
      integer:: fout=999

      logical:: ldbg = .true.
      real(8):: wkq     ! Weight of the k-q combinatiion
      real(8):: vi4pi,sqvi4pi
      real(8):: time1,time2,tstart,tend
      real(8):: kvec(3) 
      
      character(len=20)  :: sname='calcselfcmn'
      character(len=80)  :: fn_scq
!
!     Auxiliary arrays 
!
      complex(8), allocatable :: sc(:,:,:)             ! local scmn
      complex(8), allocatable :: minm(:,:,:)           ! 
      complex(8), allocatable :: mwm(:,:,:,:)

! !EXTERNAL ROUTINES: 
      character(len=10),external::int2str
      complex(8), external :: zdotc,zdotu
      external zhemm

      call cpu_time(tstart)

      if(ldbg) then
        fn_scq=trim(savdir)//trim(casename)//".sc-q"//trim(int2str(iq)) 
        open(fout,file=fn_scq,action='write') 
      endif 

      vi4pi=fourpi*vi
      sqvi4pi=sqrt(vi4pi)
      wkq = dble(weightq(iq))/dble(nqp)

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

      if(iop_core.eq.0) then
        nbands = ncg + nbmax
      else
        nbands = nbmax
      endif

      allocate(sc(ibgw:nbgw,ibgw:nbgw,nomeg),stat=ierr)
      call errmsg(ierr.ne.0,sname,"fail to allocate sc")

      do isp=1,nspin 
        nhomo = nomaxs(isp)

        do iik=1, nktot
          ! set ik, irk, jk and jrk 
          call get_kvec(isym,iik,ik,irk)
          jk=kqid(ik,iq)
          jrk=kpirind(jk)

          sc = 0.d0
          if(isc.eq.-1) then 
            call sub_calcselfc_blk(1,nbands)
          else 
            if(isc.eq.0) then 
              call sub_calcselfc_blk(1,ibgw-1)      ! low lying valence states 
              call sub_calcselfc_blk(nbgw+1,nbands) ! higher states and core states 
              if(myrank_row.eq.0) then 
                sigm_f(:,:,1:nomeg,iik,isp) = sigm_f(:,:,1:nomeg,iik,isp) + sc
              endif 
              call sub_calcselfc_blk(ibgw,nbgw)
            else
              call sub_calcselfc_blk(ibgw,nbgw)
              if(iq.eq.1) then 
                sc = sc + sigm_f(:,:,1:nomeg,iik,isp)
              endif 
            endif 
          endif 

          if(myrank_row.eq.0) then 
            sigm(:,:,1:nomeg,iik,isp) = sigm(:,:,1:nomeg,iik,isp)+sc
          endif 

          if(ldbg) then 
            write(fout,100) irk,ik,kvec(1:3) 
            write(fout,101) sc
            write(fout,102) 
          endif 

        enddo ! irk
      enddo ! isp
      deallocate(sc)

      if(ldbg) close(fout) 

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
        subroutine sub_calcselfc_blk(mst,mend) 
            implicit none 
            integer,intent(in):: mst, mend

            integer:: iom,ie1,imu,inu,cst,cend
            integer:: ncount,ierr,utype,usiz
            real(8):: time1,time2
            complex(8), allocatable :: mwm_p(:,:,:,:)
            complex(8), allocatable :: rbuf(:,:,:,:)

            if(mst.gt.mend) return 
            allocate(minm(matsiz,mst:mend,ibgw:nbgw),stat=ierr)
            call errmsg(ierr.ne.0,sname,"fail to allocate minm")

            if(mend.le.nbmax) then
              call get_minm(iop_minm,'nm',minm,ibgw,nbgw,mst,mend,ik,iq,isp)
            elseif(mst.gt.nbmax) then
              cst=mst-nbmax; cend=mend-nbmax
              call get_minm(iop_minm,'nc',minm,ibgw,nbgw,cst,cend,ik,iq,isp)
            else
              do ie1=ibgw,nbgw 
                call get_minm(iop_minm,'nm',minm(:,mst:nbmax,ie1),ie1,ie1, &
     &                mst, nbmax, ik,iq,isp)

                cst=1; cend = mend-nbmax 
                call get_minm(iop_minm,'nc',minm(:,nbmax+1:mend,ie1),ie1,ie1,&
     &                cst, cend,ik,iq,isp)
              enddo
            endif

            if(myrank_row.eq.0) then
              allocate(mwm(mst:mend,ibgw:nbgw,ibgw:nbgw,nomeg),stat=ierr)
              call errmsg(ierr.ne.0,sname,"fail to allocate mwm")
            endif
    
            if(nproc_row.eq.1.or.nomeg.eq.1) then 
              do iom=1,nomeg
                call sub_calcmwm(mwm(mst:mend,ibgw:nbgw,ibgw:nbgw,iom), &
     &                       mst,mend,iom)
              enddo ! iom
       
!! parallel in terms of frequency points has been used 
#ifdef MPI  
            else 
              allocate(rbuf(mst:mend,ibgw:nbgw,ibgw:nbgw,nomeg)) 
              allocate(mwm_p(mst:mend,ibgw:nbgw,ibgw:nbgw,iom_first:iom_last))

              do iom=iom_first,iom_last
                call sub_calcmwm(mwm_p(:,:,:,iom),mst,mend,iom)
              enddo ! iom

              !!* collect data from all processes in the same group 
              ncount=iom_last-iom_first+1
              usiz=(nbgw-ibgw+1)**2*(mend-mst+1)
              call MPI_Type_contiguous(usiz,MPI_DOUBLE_COMPLEX,utype,ierr)
              call MPI_Type_commit(utype,ierr)
              call MPI_GatherV(mwm_p,ncount,utype,rbuf,iom_cnts,iom_dspl,utype,&
     &                         0,mycomm_row,ierr)
              if(myrank_row.eq.0) then 
                do iom=1,nomeg
                  do inu=ibgw,nbgw 
                    do imu=ibgw,nbgw 
                      mwm(:,imu,inu,iom)=rbuf(:,imu,inu,iom) 
                    enddo 
                  enddo
                enddo 
              endif 
              call MPI_Type_free(utype,ierr)

              deallocate(rbuf,mwm_p)
#endif
            endif  
            deallocate(minm)

            if(myrank_row.eq.0) then
              if(isxc.eq.0) then
                call sub_calc_freqconvl(mst,mend)
              elseif(isxc.eq.2) then
                call sub_calc_cohsex(mst,mend)
              endif
              deallocate(mwm)
            endif

        end subroutine 


!======================================================================#
!     calculate correlation self-energy by the frequency convolution   #
!======================================================================#
        subroutine sub_calc_freqconvl(m0,m1)
            implicit none 
            integer, intent(in):: m0, m1
 
            integer :: imu,inu,ie2,iom,jom,icg,iat,ic
            real(8)    :: omg
            real(8)    :: enk2
            complex(8) :: scmn  ! acumlated sum over ie2
            complex(8) :: wint(nomeg)

            do inu = ibgw, nbgw       !* Loop over bands ie1
              do imu = ibgw, nbgw 
                do iom = 1, nomeg    !* Loop over frequencies for analytical continuation
                  omg=omega(iom)
                  do ie2 = m0,m1
                    if(ie2.gt.nbmax) then 
                      icg=ie2-nbmax
                      iat=corind(1,icg)
                      ic=corind(3,icg)
                      enk2=eigcore(ic,iat,isp)
                    else 
                      enk2=bande(ie2,jrk,isp)
                    endif 
                    wint=mwm(ie2,imu,inu,1:nomeg)
                    call freq_convl(iom,nomeg,omg,enk2,wint,omega,womeg,scmn) 
                    sc(imu,inu,iom) = sc(imu,inu,iom) + scmn
                  enddo
                enddo ! iom 
              enddo ! inu
            enddo ! imu
        end subroutine 

!======================================================================#
!              calculate static COHSEX correlation self-energy         #
!======================================================================#

        subroutine sub_calc_cohsex(m0,m1)
            implicit none 
            integer,intent(in):: m0, m1
            integer::imu,inu,ie2
            complex(8):: scmn 

            do inu=ibgw,nbgw 
              do imu=ibgw,nbgw
                scmn=czero 
                do ie2=m0,m1 
                  if(ie2.le.nhomo.or.ie2.gt.nbmax) then 
                    scmn = scmn - 0.5d0*mwm(ie2,imu,inu,1)
                  else 
                    scmn = scmn + 0.5d0*mwm(ie2,imu,inu,1)
                  endif 
               enddo 
               sc(imu,inu,1)=sc(imu,inu,1)+scmn
             enddo
            enddo 
        end subroutine 
!
!     This subroutine is used as a generic interface to calculate M*W*M
! 
        subroutine sub_calcmwm(xnm,m0,m1,iom)
            implicit none 
            integer,intent(in)::m0,m1,iom
            complex(8),intent(out)::xnm(m0:m1,ibgw:nbgw,ibgw:nbgw)
            integer::nmdim,imu,inu,ie2
            real(8) :: t1,t2,coefs1,coefs2
            complex(8):: xs
            complex(8), allocatable:: wm(:,:,:) 

            coefs2=singc2*vi4pi
            coefs1=singc1*sqvi4pi

            nmdim=nbandsgw*(m1-m0+1)
            if(nmdim.le.0) then 
              write(6,*) trim(sname)//":WARNING -- nmdim <=0 in sub_calcmwm"
              return 
            endif 
            allocate(wm(matsiz,m0:m1,ibgw:nbgw))
            call cpu_time(t1)
            call zhemm('l','u',matsiz,nmdim,cone,eps(:,:,iom),matsiz,  &
     &                 minm,matsiz,czero,wm,matsiz)

            do inu=ibgw,nbgw
              do imu=ibgw,nbgw 
                do ie2=m0,m1
                  xnm(ie2,imu,inu)=wkq*zdotc(matsiz,minm(:,ie2,imu),1,    &
     &                               wm(:,ie2,inu),1)
                  if(iq.eq.1.and.ie2.ge.ibgw.and.ie2.le.nbgw) then   !! add singular contributions 
                    xs = coefs2*head(iom)*qpwf_coef(imu,ie2,irk,isp)      &
     &                    *conjg(qpwf_coef(inu,ie2,irk,isp))              &
     &                 + coefs1*qpwf_coef(imu,ie2,irk,isp)                &
     &                    *zdotu(matsiz,minm(:,ie2,inu),1,epsw2,1)        &
     &                 + coefs1*conjg(qpwf_coef(inu,ie2,irk,isp))         &
     &                    *zdotc(matsiz,minm(:,ie2,imu),1,epsw1,1)
                    xnm(ie2,imu,inu)=xnm(ie2,imu,inu)+xs
                  endif  
                enddo
              enddo
            enddo
            deallocate(wm)
            call cpu_time(t2)
            time_lapack=time_lapack+t2-t1

        end subroutine 
              
      end subroutine calcselfcmn
!EOC        
