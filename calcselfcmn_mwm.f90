!BOP
!
! !ROUTINE: calcselfcmn_mwm
!
! !INTERFACE: 
      subroutine calcselfcmn_mwm(iq,isxc,isc,isym)
      
! !DESCRIPTION:
!
! This subroutine calculate the full correlation self-energy matrix
! assuming that the mwm elements are already calculated and stored 
!

! !USES:

      use bands,      only: bande,nbmaxsc,nspin,ibgw,nbgw,nbands_c 
      use core,       only: corind,eigcore,ncg
      use freq,       only: nomeg,omega,womeg
      use kpoints,    only: nirkp,nkp,nqp,kqid,kpirind,get_kvec
      use selfenergy, only: sigm,sigm_f,qpwf_coef 
      use task,       only: time_selfc,savdir,casename
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
      integer:: ierr
      integer:: isp     ! Index for spin 
      integer:: iop_minm 
      integer:: fout=999

      logical:: ldbg = .true.
      real(8):: wkq     ! Weight of the k-q combinatiion
      real(8):: vi4pi,sqvi4pi
      real(8):: time1,time2,tstart,tend
      real(8):: kvec(3) 
      
      character(len=20)  :: sname='calcselfcmn_mwm'
      character(len=80)  :: fn_scq
!
!     Auxiliary arrays 
!
      complex(8), allocatable :: sc(:,:,:)             ! local scmn
      complex(8), allocatable :: minm(:,:,:)           ! 

! !EXTERNAL ROUTINES: 
      character(len=10),external::int2str
      complex(8), external :: zdotc,zdotu
      external zhemm

      call cpu_time(tstart)

      if(ldbg) then
        fn_scq=trim(savdir)//trim(casename)//".sc-q"//trim(int2str(iq)) 
        open(fout,file=fn_scq,action='write') 
      endif 

      if(isym.eq.0) then
        nktot = nkp
      else
        nktot = nirkp
      endif

      allocate(sc(ibgw:nbgw,ibgw:nbgw,nomeg),stat=ierr)
      call errmsg(ierr.ne.0,sname,"fail to allocate sc")

      do isp=1,nspin 
        do iik=1, nktot
          ! set ik, irk, jk and jrk 
          call get_kvec(isym,iik,ik,irk)
          jk=kqid(ik,iq)
          jrk=kpirind(jk)

          sc = 0.d0
          if(isc.eq.-1) then 
            call sub_calcselfc(1,nbands_c)
          else 
            if(isc.eq.0) then 
              call sub_calcselfc(1,ibgw-1)      ! low lying valence states 
              call sub_calcselfc(nbgw+1,nbands_c) ! higher states and core states 
              if(myrank_row.eq.0) then 
                sigm_f(:,:,1:nomeg,iik,isp) = sigm_f(:,:,1:nomeg,iik,isp) + sc
              endif 
              call sub_calcselfc(ibgw,nbgw)
            else
              call sub_calcselfc(ibgw,nbgw)
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
        subroutine sub_calcselfc(mst,mend) 
        implicit none 
        integer,intent(in):: mst, mend

        integer:: iom,ie1,im,imu,inu,cst,cend
        integer:: icg,iat,ic 
        
        real(8):: em,omg
        complex(8) :: scmn
        complex(8), allocatable :: mwm_m(:,:,:),xtmp(:) 

        allocate(mwm_m(ibgw:nbgw,ibgw:nbgw,nomeg),xtmp(nomeg),stat=ierr)
        call errmsg(ierr.ne.0,sname,'Fail to allocate mwm') 

        do im=mst,mend 
          call io_mwm3('r',mwm_m,im,im,1,nomeg,isp,iq,irk,ierr) 
          call errmsg(ierr.ne.0,sname,'Fail to read mwm3') 

          if(im.gt.nbmaxsc) then
            icg = im-nbmaxsc
            iat = corind(1,icg)
            ic  = corind(3,icg)
            em  = eigcore(ic,iat,isp)
          else
            em = bande(im,jrk,isp)
          endif

          do inu=ibgw,nbgw 
            do imu=ibgw,nbgw 
              xtmp = mwm_m(imu,inu,1:nomeg)

              do iom=1,nomeg
                omg=omega(iom)
                call freq_convl(iom,nomeg,omg,em,xtmp,omega,womeg,scmn) 
                sc(imu,inu,iom) = sc(imu,inu,iom) + scmn
              enddo  ! iom
            enddo ! imu
          enddo ! inu

        enddo ! im
        end subroutine 
              
      end subroutine calcselfcmn_mwm
!EOC        
