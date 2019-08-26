!BOP
!
! !ROUTINE: calcmwm
!
! !INTERFACE: 
      subroutine calcmwm_2d(iminm,isp,iq,irk,ie2_f,ie2_l,iom_f,iom_l,&
     &                   mwm)
      
! !DESCRIPTION:
!
!     this subroutine Setup the full M*W*M matrix 
!

! !USES:
      use bands,      only: ibgw,nbgw,nbandsgw
      use barcoul,    only: iop_coul
      use bzinteg,    only: singc1co,singc2co
      use constants,  only: cone,czero,pi,fourpi,sqrt4pi
      use core,       only: ncg_c
      use dielmat,    only: eps,head,epsw1,epsw2
      use kpoints,    only: nqp,kqid,idikp,kpirind,weightq
      use mixbasis,   only: matsiz
      use struk,      only: vi 
      use task,       only: time_lapack,lrestart,time_mwm, fid_outdbg,&
     &                      fid_outgw

! !INPUT PARAMETERS:

      implicit none
    
      integer, intent(in) :: iminm      ! control how to handle Minm 
      integer, intent(in) :: isp        ! index for spin  
      integer, intent(in) :: iq         ! index for q-point 
      integer, intent(in) :: irk        ! index for irreducible k-points 
      integer, intent(in) :: ie2_f, ie2_l  ! the range for m-index 
      integer, intent(in) :: iom_f,iom_l  ! the range of frequency points 
      complex(8),intent(out):: mwm(ie2_f:ie2_l,ibgw:nbgw,iom_f:iom_l)
! !REVISION HISTORY:
!
! Created Aug. 8, 2011 by JH
!
!EOP

!BOC
! !LOCAL VARIABLES:            

      integer:: ik     ! (Counter) Runs over k-points
      integer:: ierr
      logical:: lprt=.false.
      logical:: ltest_disable_sing=.false.
      real(8):: tstart,tend
      character(len=10)::sname='calcmwm'
      complex(8), allocatable :: minm(:,:,:)

      mwm = 0.d0 

      !write(fid_outgw,*) "Calculate M*W*M by ", sname
      call cpu_time(tstart)
      ik=idikp(irk)
      if(lrestart) then 
        call io_mwm('r',mwm,ie2_f,ie2_l,iom_f,iom_l,isp,iq,irk,ierr)
        if(ierr.eq.0) then 
!          write(6,*) "Using existing mwm data"
          return  
        endif 
      endif
      
      if(ie2_f.gt.ncg_c .or. ie2_l.le.ncg_c) then  
        allocate(minm(matsiz,ie2_f:ie2_l,ibgw:nbgw),stat=ierr)
        call errmsg(ierr.ne.0,sname,"fail to allocate minm")
        if(ie2_f.gt.ncg_c) then 
          call get_minm(iminm,'nm',minm,ibgw,nbgw,                      &
     &                ie2_f-ncg_c,ie2_l-ncg_c,ik,iq,isp)
        else 
          call get_minm(iminm,'nc',minm,ibgw,nbgw,ie2_f,ie2_l,ik,iq,isp)
        endif 
        call sub_setmwm(ie2_f, ie2_l) 
        deallocate(minm)
      else 
        allocate(minm(matsiz,ie2_f:ncg_c,ibgw:nbgw),stat=ierr)
        call errmsg(ierr.ne.0,sname,"fail to allocate minm")
        call get_minm(iminm,'nc',minm,ibgw,nbgw,ie2_f,ncg_c,ik,iq,isp) 
        call sub_setmwm(ie2_f,ncg_c) 
        deallocate(minm)
       
        allocate(minm(matsiz,ncg_c+1:ie2_l,ibgw:nbgw),stat=ierr)
        call errmsg(ierr.ne.0,sname,"fail to allocate minm")
        call get_minm(iminm,'nm',minm,ibgw,nbgw,1,ie2_l-ncg_c,ik,iq,isp)
        call sub_setmwm(ncg_c+1,ie2_l)
        deallocate(minm) 
      endif 

      call io_mwm('w',mwm,ie2_f,ie2_l,iom_f,iom_l,isp,iq,irk,ierr)

      call cpu_time(tend)
      time_mwm = time_mwm + tend-tstart

      contains
!
!     This subroutine is used as a generic interface to calculate M*W*M
! 
        subroutine sub_setmwm(mst,mend)
        use anisotropy, only: iop_aniso, aniso_calc_sing_q0_1, aniten
        use bzinteg,    only: iop_q0
        implicit none 
        integer,intent(in):: mst, mend

        integer:: iom,ie1,ie2,nmdim
        real(8):: wkq     ! Weight of the k-q combinatiion
        real(8):: coefs1,coefs2
        real(8):: t1,t2
        real(8) :: thres=1.0D-10
        complex(8) :: term_singular_h, term_singular_w
        complex(8) :: accum_term_singular_h,accum_term_singular_w
        complex(8), allocatable :: wm(:,:,:)    ! (matsiz,nmdim) 
        complex(8), external :: zdotc,zdotu
        external zhemm

        coefs2=singc2co*fourpi*vi
        coefs1=singc1co*sqrt(fourpi*vi)
        ! N_c = nqp
        wkq = dble(weightq(iq))/dble(nqp)
        if(lprt) write(*,*) irk, iq, dble(weightq(iq)), wkq
  
        nmdim = nbandsgw*(mend-mst+1) 
        allocate(wm(matsiz,mst:mend,ibgw:nbgw),stat=ierr)
        call errmsg(ierr.ne.0,sname,"fail to allocate wm")

        !write(fid_outdbg, *) "iop_coul = ", iop_coul
        call cpu_time(t1)
        do iom=iom_f,iom_l
          !write(fid_outdbg,"(A27)")"MWM Singu. Contrib. "
          !write(fid_outdbg,"(A6,I5,A6,I5)")"Freq=", iom, " irk=",irk
          !write(fid_outdbg,"(A6,I5,A6,I5)")"mst= ",mst," mend=",mend
          call zhemm('l','u',matsiz,nmdim,cone,eps(:,:,iom),matsiz,  &
     &           minm,matsiz,czero,wm,matsiz)
          ! Note that ie1 are index of valence bands
          do ie1=ibgw,nbgw
            accum_term_singular_h = czero
            accum_term_singular_w = czero
            do ie2=mst,mend ! ie2=n'
              mwm(ie2,ie1,iom)=wkq*zdotc(matsiz,minm(:,ie2,ie1),1, &
     &          wm(:,ie2,ie1),1)
              if(iq.eq.1.and.iop_coul.eq.-1.and.ie1.eq.ie2-ncg_c) then
                if(iop_aniso.ne.-1.and.iop_q0.eq.1)then
                  term_singular_h = coefs2*head(iom)
                !write(fid_outdbg,"(A15,4I4,2E15.6)")"Sing.(H)(iso):",iq,irk,iom,ie1,term_singular_h
                  call aniso_calc_sing_q0_1(aniten,iom,minm(:,ie2,ie1), &
                                            term_singular_h, term_singular_w)
                !write(fid_outdbg,"(A15,4I4,2E15.6)")"Sing.(H)(ani):",iq,irk,iom,ie1,term_singular_h
                else
                  term_singular_h = coefs2*head(iom)
                  term_singular_w = &
     &             coefs1*( zdotu(matsiz,minm(:,ie2,ie1),1,epsw2,1)  &
     &                     +zdotc(matsiz,minm(:,ie2,ie1),1,epsw1,1) )    
                endif
                ! test the gap without head and wing contrib
                if(ltest_disable_sing) then
                  term_singular_h = czero
                  term_singular_w = czero
                endif

                !write(fid_outdbg,"(A10,4I4,2E15.6)") "Sing.(H):",iq,irk,iom,ie1,term_singular_h
                accum_term_singular_h=accum_term_singular_h + &
                    term_singular_h
                accum_term_singular_w=accum_term_singular_w + &
                    term_singular_w
                mwm(ie2,ie1,iom)=mwm(ie2,ie1,iom) + term_singular_h + &
     &              term_singular_w
              endif  
!              write(*,"(A20,5I5,2E18.10)") "mwm iqirkome2e1 ", &
!     &              iq,irk,iom,ie2,ie1,mwm(ie2,ie1,iom)
            enddo ! ie2
            !if(abs(term_singular_h).ge.thres)then
            !  write(fid_outdbg,"(A6,I5,A3,2E15.6)") " Band ",ie1," H ",&
            !    term_singular_h
            !endif
            !if(abs(term_singular_w).ge.thres)then
            !  write(fid_outdbg,"(A14,2E15.6)")" W ",term_singular_w
            !endif
          enddo ! ie1
        enddo ! iom 
        call cpu_time(t2)
        time_lapack=time_lapack+t2-t1
        deallocate(wm) 
        end subroutine 
              
      end subroutine calcmwm_2d
!EOC        
