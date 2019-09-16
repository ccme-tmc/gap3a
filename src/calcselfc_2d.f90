!BOP
!
! !ROUTINE: calcselfc
!
! !INTERFACE: 
      subroutine calcselfc_2d(iq,iminm)
      
! !DESCRIPTION:
!
!This subroutine calculate the correlation self-energy, this
!includes the weight for the BZ integration and the frequency integration
!is hardcoded inside
!

! !USES:

      use bands,      only: bande,nspin,nomaxs,numins,ibgw,nbgw, &
     &                      nbandsgw,nbands_c 
      use barcoul,    only: iop_coul
      use core,       only: corind,eigcore,ncg_c 
      use freq,       only: nomeg,omega,womeg
      use kpoints,    only: nirkp,nqp,kqid,idikp,kpirind,weightq
      use selfenergy, only: iop_sxc,sigc_q, sigc 
      use task,       only: time_selfc

      use modmpi,     only: nproc_row,myrank_row,mycomm_row, &
     &                      nproc_3rd,myrank_3rd,mycomm_3rd, &
     &                      nproc_max
#ifdef MPI 
      use modmpi,     only: mpi_set_range, mpi_sum_array, &
     &                      mpi_gather_array
#endif 




! !INPUT PARAMETERS:

      implicit none
      
      integer, intent(in) :: iq       ! index for q-point 
      integer, intent(in) :: iminm  ! control how to handle Minm 
! !REVISION HISTORY:
!
! Created 01.06.2011 by JH
!
!EOP

!BOC
! !LOCAL VARIABLES:            

      integer:: ie1     ! (Counter) Runs over bands.
      integer:: ik     ! (Counter) Runs over k-points
      integer:: irk
      integer:: iom     ! (Counter) Runs over frequencies.
      integer:: jk     ! (Counter) Runs over k-points
      integer:: jrk
      integer:: nomx
      integer:: ierr
      integer:: isp     ! Index for spin 
      integer:: irk_f, irk_l, irk_cnts(0:nproc_max), irk_dspl(0:nproc_max) 
      integer:: ie2_f, ie2_l
      integer:: nd_sc 

      logical:: ldbg=.false.
      real(8):: tstart,tend
      character(len=20)::sname='calcselfc'
      complex(8), allocatable:: sc_p(:,:,:), sc_a(:,:,:) !! local arrays needed for parallelized calculation of sigc
      complex(8), allocatable :: mwm(:,:,:) 
       
      call cpu_time(tstart)


      if(ldbg) call linmsg(6,'-',sname)

      sigc_q = 0.d0 

      if(iop_sxc.eq.0) then 
        nd_sc = nomeg
      elseif(iop_sxc.eq.2) then 
        nd_sc = 3
      elseif(iop_sxc.eq.3) then 
        nd_sc = 1
      endif 

!!      preparation for parallelization over irk and m
#ifdef MPI 
      call mpi_set_range(nproc_row,myrank_row,nirkp,1,irk_f,irk_l, &
     &                   irk_cnts,irk_dspl)
      call mpi_set_range(nproc_3rd,myrank_3rd,nbands_c,1,ie2_f,ie2_l)
      if(ldbg) write(6,*) "nproc_row,irk_f/l=", nproc_row,irk_f,irk_l 
      if(ldbg) write(6,*) "nproc_3rd,ie2_f/l=", nproc_3rd,ie2_f,ie2_l 
#else 
      irk_f = 1
      irk_l = nirkp
      ie2_f = 1
      ie2_l = nbands_c
#endif 

      allocate( sc_p(nd_sc,ibgw:nbgw,irk_f:irk_l), &
     &          sc_a(nd_sc,ibgw:nbgw,1:nirkp),     &
     &          mwm(ie2_f:ie2_l,ibgw:nbgw,nomeg),  &
     &          stat=ierr)
      call errmsg(ierr.ne.0,sname,"fail to allocate mwm,sc_a,sc_p")
      sc_a = 0.d0 


      do isp=1,nspin 
        nomx = nomaxs(isp)
  
        sc_p = 0.d0 
        do irk=irk_f, irk_l
          ik=idikp(irk)
          jk=kqid(ik,iq)
          jrk=kpirind(jk)
          call calcmwm_2d(iminm,isp,iq,irk,ie2_f,ie2_l,1,nomeg,mwm) 
          if(iop_sxc.eq.0) then 
            call sub_sigc_gw 
          elseif(iop_sxc.eq.2) then
            call sub_sigc_cohsex
          elseif(iop_sxc.eq.3) then
            call sub_sigc_sex
          endif  

#ifdef MPI         
          !! sum over the parallelization over ie2 
          if(nproc_3rd.gt.1) then
            call mpi_sum_array(0,sc_p(:,:,irk),nd_sc,nbandsgw,mycomm_3rd)
          endif
#endif 
        enddo ! irk

        !! gather data for different irk's calculated on different processes
#ifdef MPI 
        if(nproc_row.gt.1.and.myrank_3rd.eq.0) then
          call mpi_gather_array(0,sc_p,irk_l-irk_f+1,nd_sc,nbandsgw, &
     &            sc_a,irk_cnts,irk_dspl,mycomm_row)
        else
#endif 
          sc_a = sc_p

#ifdef MPI
        endif
#endif 

        if(myrank_3rd.eq.0.and.myrank_row.eq.0) then
          do irk=1,nirkp
            do ie1 = ibgw, nbgw
              sigc_q(:,ie1,irk,isp)=sigc_q(:,ie1,irk,isp) + sc_a(:,ie1,irk)
            enddo
          enddo
        endif ! myrank_3rd.eq.0
      enddo ! isp

      sigc = sigc + sigc_q 
      deallocate(mwm,sc_p,sc_a) 

      call cpu_time(tend)
      time_selfc = time_selfc + tend - tstart 

      return

      contains 

!======================================================================#
!     calculate correlation self-energy by the frequency convolution   #
!======================================================================#
        subroutine sub_sigc_gw 
        integer:: ie1,ie2,iom,iat,ic
        real(8):: omg,enk2
        complex(8):: sc, wint(nomeg), sc_n

        do ie1 = ibgw, nbgw       !* Loop over bands ie1
          do iom = 1, nomeg    !* Loop over frequencies for analytical continuation
            omg=omega(iom)
            do ie2=ie2_f,ie2_l
              if(ie2.le.ncg_c) then !! core states 
                iat=corind(1,ie2)
                ic =corind(3,ie2)
                enk2=eigcore(ic,iat,isp)
              else 
                enk2=bande(ie2-ncg_c,jrk,isp) 
              endif 
              wint=mwm(ie2,ie1,1:nomeg)
              call freq_convl(iom,nomeg,omg,enk2,wint,omega,womeg,sc) 
              sc_p(iom,ie1,irk)=sc_p(iom,ie1,irk) + sc
            enddo
          enddo ! iom
        enddo ! ie1
      end subroutine 

!======================================================================#
!              calculate static COHSEX correlation self-energy         #
!======================================================================#

      subroutine sub_sigc_cohsex
        implicit none 
        integer::ie1,ie2
        complex(8):: sc(3) 

        do ie1=ibgw,nbgw 
          sc = 0.0
          do ie2=ie2_f,ie2_l 
            if(ie2.le.nomx+ncg_c) then         !! occupied
              sc(1) = sc(1) - 0.5d0*mwm(ie2,ie1,1)
              sc(2) = sc(2) - mwm(ie2,ie1,1)
            else                               !! unoccupied 
              sc(1) = sc(1) + 0.5d0*mwm(ie2,ie1,1)
            endif 
            sc(3)=sc(3) + 0.5d0*mwm(ie2,ie1,1)
          enddo 
          sc_p(1:3,ie1,irk)=sc_p(1:3,ie1,irk)+sc
        enddo 
      end subroutine 

      subroutine sub_sigc_sex
        implicit none
        integer::ie1,ie2
        complex(8):: sc

        do ie1=ibgw,nbgw
          sc = 0.0
          do ie2=ie2_f,ie2_l
            if(ie2.le.nomx+ncg_c) then         !! occupied
              sc = sc - mwm(ie2,ie1,1)
            endif
          enddo
          sc_p(1,ie1,irk)=sc_p(1,ie1,irk)+sc 
        enddo
      end subroutine

              
      end subroutine calcselfc_2d
!EOC        
