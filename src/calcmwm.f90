!BOP
!
! !ROUTINE: calcmwm
!
! !INTERFACE: 
      subroutine calcmwm(iminm,isp,iq,irk,ie2_f,ie2_l,iom_f,iom_l,&
     &                   mwm)
      
! !DESCRIPTION:
!
!     this subroutine Setup the full M*W*M matrix 
!

! !USES:

      use bands,      only: ibgw,nbgw,nbandsgw
      use barcoul,    only: iop_coul
      use bzinteg,    only: singc1,singc2
      use constants,  only: cone,czero,pi,fourpi,sqrt4pi
      use core,       only: ncg_c
      use dielmat,    only: eps,head,epsw1,epsw2
      use kpoints,    only: nqp,kqid,idikp,kpirind,weightq
      use mixbasis,   only: matsiz
      use struk,      only: vi 
      use task,       only: time_lapack,lrestart,time_mwm

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
      real(8):: tstart,tend
      character(len=10)::sname='calcmwm'
      complex(8), allocatable :: minm(:,:,:)

      mwm = 0.d0 

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
        implicit none 
        integer,intent(in):: mst, mend

        integer:: iom,ie1,ie2,nmdim
        real(8):: wkq     ! Weight of the k-q combinatiion
        real(8):: coefs1,coefs2
        real(8):: t1,t2
        complex(8), allocatable :: wm(:,:,:)    ! (matsiz,nmdim) 
        complex(8), external :: zdotc,zdotu
        external zhemm

        coefs2=singc2*fourpi*vi
        coefs1=singc1*sqrt(fourpi*vi)
        wkq = dble(weightq(iq))/dble(nqp)
  
        nmdim = nbandsgw*(mend-mst+1) 
        allocate(wm(matsiz,mst:mend,ibgw:nbgw),stat=ierr)
        call errmsg(ierr.ne.0,sname,"fail to allocate wm")

        call cpu_time(t1)
        do iom=iom_f,iom_l
          call zhemm('l','u',matsiz,nmdim,cone,eps(:,:,iom),matsiz,  &
     &           minm,matsiz,czero,wm,matsiz)
          do ie1=ibgw,nbgw
            do ie2=mst,mend
              mwm(ie2,ie1,iom)=wkq*zdotc(matsiz,minm(:,ie2,ie1),1, &
     &           wm(:,ie2,ie1),1)

              if(iq.eq.1.and.iop_coul.eq.0.and.ie1.eq.ie2-ncg_c) then 
                mwm(ie2,ie1,iom) = mwm(ie2,ie1,iom)                     &
     &           + coefs2*head(iom)                                     &
     &           + coefs1*( zdotu(matsiz,minm(:,ie2,ie1),1,epsw2,1)     &
     &                     +zdotc(matsiz,minm(:,ie2,ie1),1,epsw1,1) )    
              endif  
            enddo ! ie1
          enddo ! ie2
        enddo ! iom 
        call cpu_time(t2)
        time_lapack=time_lapack+t2-t1
        deallocate(wm) 
        end subroutine 
              
      end subroutine calcmwm
!EOC        
