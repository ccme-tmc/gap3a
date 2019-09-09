!BOP
!
! !ROUTINE: calcselfx
!
! !INTERFACE: 
      subroutine calcselfx_2d(iq,iop_minm)
      
! !DESCRIPTION:
!
! This subroutine calculate the diagonal exchange self energy for given q 
!

! !USES:

      use bands,      only: bande,efermi,nbmax,nspin,        &
     &                      nomaxs,numins,ibgw,nbgw,nbandsgw 
      use barcoul,    only: iop_coul,zcut_coul
      use bzinteg,    only: singc2ex,kiw
      use constants,  only: cone,czero,pi,fourpi,sqrt4pi
      use core,       only: corind, eigcore,ncg,iop_core,ncg_x 
      use kpoints,    only: nirkp,nqp,kqid,idikp,kpirind
      use minmmat,    only: mblksiz
      use mixbasis,   only: matsiz
      use selfenergy, only: sigx_q,sigx
      use struk,      only: vi
      use task,       only: time_selfx,lrestart 

      use modmpi,     only: nproc_row,myrank_row,mycomm_row, &
     &                      nproc_3rd,myrank_3rd,mycomm_3rd, & 
     &                      nproc_max
#ifdef MPI 
      use modmpi,     only: mpi_set_range, mpi_sum_array, &
     &                      mpi_gather_array
#endif 

! !INPUT PARAMETERS:

      implicit none
 
      integer, intent(in) :: iq                !! index for q-point 
      integer, intent(in) :: iop_minm          !! control how to handle Minm

! !REVISION HISTORY:
!
!  Created 23.06.05 by RGA.
! Modified 09 Sept 2019 by MYZ
!
!EOP

!BOC
! !LOCAL VARIABLES:            

      integer :: ie1,ie2     ! index for bands.
      integer :: ik,irk      ! index for k-points
      integer :: isp         ! index for spin 
      integer :: jk,jrk      ! index for k-q
      integer :: ie2_f,ie2_l    ! lower and upper bound for the summation over bands when parallelized
      integer :: nomx,nbands
      real(8) :: tstart,tend

      integer :: irk_f, irk_l, nirk_p, irk_cnts(0:nproc_max),irk_dspl(0:nproc_max)

      complex(8),allocatable:: minm(:,:,:)
      complex(8),allocatable:: sx_p(:,:)
      complex(8),allocatable:: sx_a(:,:) !! arrays of the x-selfenergy used in the parallelization 

      integer :: ierr
      character(len=10)::sname='calcselfx'
      logical :: ldbg = .true.

! !EXTERNAL ROUTINES: 
      complex(8), external :: zdotc
      external zhemm

      call cpu_time(tstart) 

#ifdef DEBUG
      ldbg = .true.
#endif 



      if(ldbg) call linmsg(6,'-','calcselfx')

      sigx_q = 0.d0 

      do isp=1,nspin 
        nomx = nomaxs(isp)
        nbands = ncg_x + nomx

!!      preparation for parallelization over irk and m 
#ifdef MPI
        call mpi_set_range(nproc_row,myrank_row,nirkp,1,irk_f,irk_l, &
     &                   irk_cnts,irk_dspl)
        call mpi_set_range(nproc_3rd,myrank_3rd,nbands,1,ie2_f,ie2_l)
#else
        irk_f = 1
        irk_l = nirkp 
        ie2_f = 1
        ie2_l = nbands 
#endif 

        if(ldbg) then 
          write(6,*) "nproc_row,irk_f,irk_l=",nproc_row, irk_f, irk_l
          write(6,*) "nproc_3rd,ie2_f,ie2_l=",nproc_3rd, ie2_f, ie2_l
        endif 

        allocate(minm(matsiz,ie2_f:ie2_l,ibgw:nbgw),  &
     &           sx_a(ibgw:nbgw,nirkp),                 &
     &           sx_p(ibgw:nbgw,irk_f:irk_l),stat=ierr)
        call errmsg(ierr.ne.0,sname,"fail to allocate minm")
        minm = czero
        sx_a = 0.d0 
        sx_p = 0.d0
        

        do irk=irk_f,irk_l 
          ik=idikp(irk)  
          jk=kqid(ik,iq)
          jrk=kpirind(jk)
           
          !* calculate minm and/or minc
          if(ie2_l.le.nomx) then
            call get_minm(iop_minm,'nm',minm,ibgw,nbgw,ie2_f,ie2_l,ik,iq,isp)
          elseif(ie2_f.gt.nomx) then
            call get_minm(iop_minm,'nc',minm,ibgw,nbgw,   &
     &                  ie2_f-nomx,ie2_l-nomx,ik,iq,isp) 
          else
            do ie1=ibgw,nbgw 
              call get_minm(iop_minm,'nm',minm(:,ie2_f:nomx,ie1),ie1,ie1,   &
     &                  ie2_f,nomx,ik,iq,isp)
              call get_minm(iop_minm,'nc',minm(:,nomx+1:ie2_l,ie1),ie1,ie1,&
     &                  1,ie2_l-nomx,ik,iq,isp) 
            enddo 
          endif

          call sub_setselfx(sx_p(:,irk),ie2_f,ie2_l)

#ifdef MPI          
          if(nproc_3rd.gt.1) then  
            call mpi_sum_array(0,sx_p(:,irk),nbandsgw,mycomm_3rd)
          endif 
#endif 
        enddo   ! irk

        !! gather data for different irk's calculated on different processes
#ifdef MPI 
        if(nproc_row.gt.1.and.myrank_3rd.eq.0) then 
          call mpi_gather_array(0,sx_p,irk_l-irk_f+1,nbandsgw, &
     &            sx_a,irk_cnts,irk_dspl,mycomm_row)
        else 
#endif 
          sx_a = sx_p 

#ifdef MPI
        endif  
#endif 
        if(myrank_3rd.eq.0.and.myrank_row.eq.0) then 
          sigx_q(:,:,isp) = sigx_q(:,:,isp) + sx_a
        endif ! myrank_3rd.eq.0)
        deallocate(minm,sx_p,sx_a) 
      enddo ! isp

      sigx = sigx + sigx_q 

      call cpu_time(tend) 
      time_selfx=time_selfx+tend-tstart 

      return
 101  format(2i5,3f12.4) 
    
      contains 

        subroutine sub_setselfx(sx,ie2_f,ie2_l)
        integer,intent(in)::ie2_f,ie2_l
        complex(8),intent(out)::sx(ibgw:nbgw)
        integer::ie1,ie2
        complex(8) :: mvm
        real(8) :: wt,sxs2,wts

        sxs2= fourpi*vi*zcut_coul
        do ie1=ibgw,nbgw

          do ie2=ie2_f,ie2_l

            if(ie2.le.nomx) then 
              wt=kiw(ie2,jrk,isp)
            else 
              wt=kiw(1,jrk,isp)
            endif 
            mvm = - zdotc(matsiz,minm(:,ie2,ie1),1,minm(:,ie2,ie1),1)
            sx(ie1)=sx(ie1) + mvm*wt
            ! contribution from singular term
            if(iop_coul.eq.-1.and.iq.eq.1.and.ie1.eq.ie2.and.ie1.le.nomx) then 
              wts = wt*nqp
              if(ldbg) write(6,*) "ie1=",ie1,"wts=",wts
              sx(ie1) = sx(ie1) - sxs2*singc2ex*wts
            endif 
          enddo
        enddo
        end subroutine ! internal subroutine sub_setselfx
              
      end subroutine calcselfx_2d
!EOC        
