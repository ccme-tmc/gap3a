!BOP
!
! !ROUTINE: calceps
!
! !INTERFACE:
      subroutine calceps(iq,iom_f,iom_l,isym,isc,iop,lread)

! !DESCRIPTION:
!
! This subroutine calculates the dielectric matrix in the RPA approximation 
! according to eq. \ref{polmatdefom} used in self-consistent GW
!
! !USES:

      use bands,       only: bande,nbmaxpol,nspin,fspin,nomaxs,numins,  &
     &                       eminpol,emaxpol,metallic,nomax,numin 
      use barcoul,     only: iop_coul_c
      use bzinteg,     only: kcw
      use constants,   only: cone, czero,pi,hev
      use core,        only: corind, ncg_p, eigcore 
      use crpa,        only: mill
      use dielmat,     only: eps,epsw1,epsw2,q0_eps,mask_eps, &
     &                       c0_head,emac,head
      use freq,        only: omega,iop_freq
      use kpoints,     only: kqid,kpirind,idikp,wkir,nkp,nirkp
      use mixbasis,    only: matsiz
      use minmmat,     only: mblksiz
      use mommat,      only: mmatvv,mmatcv,init_mommat,end_mommat,&
     &                       iop_mommat
      use selfenergy,  only: qpwf_coef
      use struk,       only: nat,vi
      use task,        only: casename,time_lapack,time_eps,l_save_dielmat, &
     &                       fid_outgw, fid_outdbg,time_aniso
      use modmpi
      use anisotropy,  only: ten_a_ani,ten_p_ani,&
     &                       vec_a_ani,vec_b_ani,vec_t_ani,vec_u_ani,&
     &                       iop_aniso,proj_head_on_ylm,&
     &                       calc_h_w_inv_ang_grid,angint_eps_sph

! !INPUT PARAMETERS:
      
      implicit none

      integer, intent(in):: iq          ! Index of the q vector for which the polarization matrix is calculated
      integer, intent(in):: iom_f,iom_l ! the range of frequency points
      integer, intent(in):: isym        ! = 0/1 if sum over all/irreducible k-points 
      integer, intent(in):: isc         ! indicates whether this is the first iteration in QSGW
      integer, intent(in):: iop         ! determine what is really calculated 
                                        !  0/1/2 - eps / eps^{-1} / eps^{-1}-1 / read from the file 
      logical,intent(in) :: lread       ! read from the eps file if possible
                                 

! !LOCAL VARIABLES:

      integer :: iat,icg,ic         ! indices for core states of a given atom.
      integer :: ie1                ! index occupied states   
      integer :: ie2                ! index for unoccupied states 
      integer :: nktot,iik,ik,irk,jk,jrk   ! indices for the k-points.
      integer :: iom                ! Counter: Runs over frequencies.
      integer :: isp                ! Counter: runs over spin
      integer :: imats              ! Counter: runs over matsiz
      integer :: i,j                  ! Counter: Cartesian axis
      integer :: ie12
      integer :: ierr
      integer :: im  ! index for basis set 
      integer :: nvbm,ncbm,nmdim,cmdim
      integer :: ie2_f,ie2_l  
      integer :: ie1_f,ie1_l  
      integer :: ik_f, ik_l 
      
      logical :: ldbg=.true.
      logical :: l_eps_invd=.false.  ! label if eps,epsw and head inverted
      logical :: l_body_invd=.false. ! label if eps (body only) inverted

      integer :: nblk,iblk,m_f,m_l
      integer :: nomg
      integer :: iop_minm

      real(8) :: edif,coefwing,aux,ffact,omgsq,coefcoul
      complex(8):: coefks,cedif,caux,epsw1_tmp,epsw2_tmp,head_tmp,ccoefcoul
      complex(8):: ten_a_ani_tmp(3,3)

      real(8) :: time1,time2,time3,time4,tstart,tend
      character(len=15)::sname='calceps'

      real(8),   allocatable :: enk(:)        !! local array for eigen-energies 
      complex(8),allocatable :: tmat(:,:),pm(:),wtmp(:),minm(:,:,:) 
      complex(8), allocatable :: u_ani_iom(:,:),vec_u_tmp(:,:)

      character(len=10),external:: int2str 
      real(8),external :: coul_coef

      real(8),parameter :: sqr3=1.732050807568877294d0

!
!EOP
!BOC

      if(ldbg) call linmsg(6,'-',sname)
      call cpu_time(tstart)

      if(lread) then 
        call io_eps('r',iq,iom_f,iom_l,ierr)
        if(ierr.eq.0) then 
          write(fid_outgw,*) "Restart mode: read eps from files"
          return 
        else 
          write(fid_outgw,*) "WARNING: fail to read eps from files"
        endif 
      endif 
!
!     isc.eq.0 indicates either it is a G0W0 calculation or it is
!     the first iteration of QSGW calculations, one just need to reads
!     from Minm file directly, in the case of isc.gt.0, indicating
!     a QSGW calculation, Minm are read from the the file and transformed 
!     in terms of both "n" and "m" by qpwf_c,ccoefwoef 
!
      if(isc.lt.0) then 
        iop_minm = -1
      elseif(isc.eq.0) then 
        iop_minm = 0
      else
        iop_minm = 2
      endif 

      if(isym.eq.0) then 
        nktot = nkp 
      else
        nktot = nirkp
      endif 
      nomg = iom_l - iom_f + 1
     
      ! calculate q-dep. BZint weights
      call bz_calcqdepw(iq)  

      time_aniso = 0.0D0
      ! first calculate momentum matrix and the head
      if(iq.eq.1.and.iop_coul_c.eq.-1) then
        call init_mommat(1,nomax,numin,nbmaxpol,nirkp,nspin)
        if(iop_aniso.ne.-1) then
        !  call cpu_time(time1)
        !  call init_aniso(iom_f,iom_l)
        !  call cpu_time(time2)
        !  time_aniso = time_aniso + time2 - time1
          allocate(u_ani_iom(3,matsiz))
        endif

        if(iop_mommat.eq.0) then
          call calcmommat(0,1,nomax,numin,nbmaxpol,0)
        else
          call w2k_readmommat(1,nomax,numin,nbmaxpol)
        endif
        call calchead(1,nirkp,iom_f,iom_l)
      endif
        
      ! set the parallelization over k (within the row)
      call mpi_set_range(nproc_row,myrank_row,nktot,1,ik_f,ik_l)
      write(fid_outgw,*) "calceps: myrank_row,ik_f,ik_l =",myrank_row,ik_f,ik_l

      ! TODO choose q0_eps
      coefcoul = coul_coef(q0_eps,iop_coul_c)
      ccoefcoul = cmplx(coefcoul,0.0D0,8)
      coefwing = sqrt(vi)

      do isp=1,nspin  
        do iik=ik_f, ik_l 
          ! set ik, irk, jk and jrk 
          if(isym.eq.0) then 
            ik=iik 
            irk=kpirind(ik)
            coefks = cone*fspin 
          else 
            irk=iik
            ik=idikp(irk)
            coefks = wkir(irk)*cone*fspin
          endif 

          jk = kqid(ik,iq)
          jrk= kpirind(jk)
          nvbm = nomaxs(isp) 
          ncbm = numins(isp) 

          !! set the local array for band energies  
          allocate(enk(ncg_p + nbmaxpol))
          ! set core state energy part
          do icg=1,ncg_p       
            iat = corind(1,icg)
            ic  = corind(3,icg)
            enk(icg) = eigcore(ic,iat,isp) 
          enddo
          ! set valence band energy part
          enk(ncg_p+1:ncg_p+nbmaxpol) = bande(1:nbmaxpol,irk,isp) 

          !! find the index for the band whose energy is higher than eminpol  
          ie1_l = ncg_p + nvbm 
          do ie1=1,ie1_l
            ie1_f = ie1
            if(enk(ie1) > eminpol) exit 
          enddo
          if(ldbg) write(fid_outgw,*) " ie1_f=",ie1_f, "ie1_l=", ie1_l

          !! set mask_eps that defines the transitions to be skipped 
          call crpa_setmask(irk,jrk,isp) 

          !!  determine the number of blocks used in minm operation 
          nblk=max((nbmaxpol-ncbm+1)/mblksiz,1) 

          if(ldbg) write(fid_outgw,*) "nblk=",nblk 

          do iblk=1,nblk                     !! loop over m-blocks in Minm
            m_f = ncbm+(iblk-1)*mblksiz
            m_l = m_f + mblksiz-1
            if(iblk.eq.nblk) m_l = nbmaxpol

            call mpi_set_range(nproc_3rd,myrank_3rd,m_l-m_f+1,m_f,ie2_f,ie2_l)
  
            nmdim = (ie2_l-ie2_f+1)*(ie1_l-ie1_f+1) ! all possible n-m combinations

            allocate(tmat(1:matsiz,nmdim),& 
     &               pm(nmdim),& 
     &               wtmp(matsiz),&
     &               minm(matsiz,ie2_f:ie2_l,ie1_f:ie1_l), & 
     &               vec_u_tmp(3,nmdim), &
     &               stat=ierr)
            call errmsg(ierr.ne.0,sname,"Fail to allocate tmat,pm")

            if(ie1_f.le.ncg_p) then 
              call get_minm(iop_minm,'cm',minm(:,:,ie1_f:ncg_p),  &
     &                      ie1_f, ncg_p, ie2_f, ie2_l, ik, iq, isp)
              call get_minm(iop_minm,'nm',minm(:,:,ncg_p+1:ie1_l),  &
     &                      1, ie1_l-ncg_p, ie2_f, ie2_l, ik, iq,isp)
            else  
              call get_minm(iop_minm,'nm',minm,ie1_f-ncg_p,ie1_l-ncg_p,&
     &                      ie2_f,ie2_l,ik,iq,isp)
            endif 

            !! Calculate pm(ie12) = pm(ie1,ie2) = p_{ie1,ie2}/(e_2-e_1)
            !! Needed for the wings
            if(iq.eq.1.and.iop_coul_c.eq.-1) then
              ie12=0
              do ie1=ie1_f,ie1_l       ! the index for occ. states
                do ie2=ie2_f,ie2_l     ! the index for unocc. states. 
                  ie12=ie12+1
                  edif = enk(ie2+ncg_p) - enk(ie1)

                  if(abs(edif).gt.1.0d-10)then
                    aux=coefwing/edif
                  else 
                    aux=0.d0
                  endif 
                  
                  caux=cmplx(aux,0.0D0,8)
                  if(ie1.le.ncg_p) then
                    pm(ie12)=sqrt(ccoefcoul)*caux &
     &                  *sum(mmatcv(:,ie1,ie2,irk,isp)*q0_eps)
                    if(iop_aniso.ne.-1)then
                      vec_u_tmp(:,ie12)=caux*mmatcv(:,ie1,ie2,irk,isp)
                    endif ! iop_aniso.ne.-1
                  else 
                    pm(ie12)=sqrt(ccoefcoul)*caux &
     &                  *sum(mmatvv(:,ie1-ncg_p,ie2,irk,isp)*q0_eps)
                    if(iop_aniso.ne.-1)then
                      vec_u_tmp(:,ie12)=caux*mmatvv(:,ie1-ncg_p,ie2,irk,isp)
                    endif ! iop_aniso.ne.-1
                  endif 
                enddo 
              enddo 
            endif   ! iq.eq.1

            do iom=iom_f,iom_l
              ! TODO omgsq not used?
              omgsq=omega(iom)**2
              ie12=0
              do ie1=ie1_f, ie1_l 
                do ie2=ie2_f, ie2_l 
                  ie12=ie12+1
                  tmat(:,ie12)=minm(:,ie2,ie1)*kcw(ie1,ie2,ik,iom,isp) &
     &                         *mask_eps(ie2,ie1)
                enddo
              enddo  

              call cpu_time(time1)
              ! TODO is coulomb interacton included in minm?
              call zgemm('n','c',matsiz,matsiz,nmdim,-coefks,tmat,matsiz, &
     &            minm,matsiz,cone,eps(:,:,iom),matsiz)   

              if(iq.eq.1.and.iop_coul_c.eq.-1) then 
                call zgemv('n',matsiz,nmdim,-coefks,tmat,matsiz,pm,1,czero,wtmp,1)
                !call zgemv('n',matsiz,nmdim,coef,tmat,matsiz,pm,1,cone,epsw1(:,iom),1)
                if(iop_aniso.ne.-1)then
                  call cpu_time(time3)
                  call zgemm('n','t',3,matsiz,nmdim,coefks,&
     &                 vec_u_tmp,3,tmat,matsiz,czero,u_ani_iom,3)
                  vec_u_ani(:,:,iom)=vec_u_ani(:,:,iom)+u_ani_iom(:,:)
                  call cpu_time(time4)
                  time_aniso = time_aniso + time4 - time3
                endif
                epsw1(:,iom)=epsw1(:,iom)+wtmp

                if(iop_freq.eq.2) then !! real freq 
                  ie12=0
                  do ie1=ie1_f, ie1_l 
                    do ie2=ie2_f, ie2_l 
                      ie12=ie12+1
                      tmat(:,ie12)=minm(:,ie2,ie1) &
     &                  *conjg(kcw(ie1,ie2,ik,iom,isp)) &
     &                  *mask_eps(ie2,ie1) 
                    enddo
                  enddo
                  call zgemv('n',matsiz,nmdim,-coefks,tmat,matsiz,pm,1,czero,wtmp,1)
                  !call zgemv('n',matsiz,nmdim,coefks,tmat,matsiz,pm,1,cone,epsw2(:,iom),1)
                endif
                if(iop_aniso.ne.-1)then
                  call cpu_time(time3)
                  call zgemm('n','t',3,matsiz,nmdim,coefks,&
     &                 vec_u_tmp,3,tmat,matsiz,czero,u_ani_iom,3)
                  vec_t_ani(:,:,iom)=vec_t_ani(:,:,iom)+conjg(u_ani_iom)
                  call cpu_time(time4)
                  time_aniso = time_aniso + time4 - time3
                endif
                epsw2(:,iom)=epsw2(:,iom)+conjg(wtmp)
              endif

              call cpu_time(time2)
              time_lapack=time_lapack+time2-time1
            enddo ! iom
            deallocate(minm,tmat,pm,wtmp,vec_u_tmp)

          enddo ! iblk  
          deallocate(enk) 

        enddo ! iik
      enddo ! isp

#ifdef MPI 
      if(nproc_ra3.gt.1) then 
        write(fid_outgw,*) "Collect eps: myrank_ra3=",myrank_ra3,"comm=",mycomm_ra3
        call mpi_sum_array(0,eps,matsiz,matsiz,nomg,mycomm_ra3)
        if(iq.eq.1) then
          call mpi_sum_array(0,epsw1,matsiz,nomg,mycomm_ra3)
          call mpi_sum_array(0,epsw2,matsiz,nomg,mycomm_ra3)
          ! TODO check if the mpi_sum_array works with vec_u_ani etc
          if(iop_aniso.ne.-1)then
            call cpu_time(time1)
            call mpi_sum_array(0,vec_u_ani,3,matsiz,nomg,mycomm_ra3)
            call mpi_sum_array(0,vec_t_ani,3,matsiz,nomg,mycomm_ra3)
            call cpu_time(time2)
            time_aniso = time_aniso + time2 - time1
          endif
        endif
      endif
      write(fid_outgw,"(A28,I3,A10,I3)")"eps/w/h Sum done: myrank",&
     & myrank_ra3, " in comm ",mycomm_ra3
      call MPI_Barrier(mycomm_ra3, ierr)
      write(fid_outgw,"(A28,I3,A10,I3,A2,I3)")"MPI_Barrier of myrank",&
     & myrank_ra3, "in comm ",mycomm_ra3, " :", ierr
#endif

      ! for anisotropy, calculate wings here with vec_u and vec_t
      if(iq.eq.1.and.iop_coul_c.eq.-1.and.iop_aniso.ne.-1) then
        call cpu_time(time1)
        write(fid_outdbg, *) "diff epsw from iso and aniso"
        write(fid_outdbg, "(A3,A4,A2,4A13)") "iom","imats", "ST",&
     &                      "ReW1","ImW1","ReW2","ImW2"
        do iom=iom_f, iom_l
          do imats=1, matsiz
            epsw1_tmp=epsw1(imats,iom)
            epsw2_tmp=epsw2(imats,iom)
            epsw1(imats,iom) = - sqrt(ccoefcoul) * & 
     &          sum(vec_u_ani(:,imats,iom)*cmplx(q0_eps(:),0.0D0,8))
            epsw2(imats,iom) = - sqrt(ccoefcoul) * & 
     &          sum(vec_t_ani(:,imats,iom)*cmplx(q0_eps(:),0.0D0,8))
            write(fid_outdbg, "(I3,I4,A2,4e13.5)") iom, imats,'O',&
     &                            epsw1_tmp,epsw2_tmp
            write(fid_outdbg, "(I3,I4,A2,4e13.5)") iom, imats,'N',&
     &                            epsw1(imats,iom),epsw2(imats,iom)
            write(fid_outdbg, "(I3,I4,A2,2L26)") iom, imats,'D',&
     &               abs(epsw1(imats,iom)-epsw1_tmp).lt.1.0D-12, &
     &               abs(epsw2(imats,iom)-epsw2_tmp).lt.1.0D-12
          enddo
        enddo
        call cpu_time(time2)
        time_aniso = time_aniso + time2 - time1
      endif

      !! add "1" to the diagonal elements
      if(myrank_ra3.eq.0) then 
        do iom=iom_f,iom_l
          do im=1,matsiz
            eps(im,im,iom) = eps(im,im,iom) + cone
          enddo
         enddo
      endif

      if(iq.eq.1.and.iop_coul_c.eq.-1) then
        call end_mommat
      endif 

      if(iop.eq.0) then  !! calculate only eps 
        call sub_write_emac(0)  
      else               !!  calculate inveps 
        if(myrank_ra3.eq.0) then
          call sub_invbody
          call sub_calcinveps 
          call sub_write_emac(1) 
        endif 
#ifdef MPI 
        if(nproc_ra3.gt.1) then
          write(fid_outgw,*) "bcast eps: myrank_ra3 = ",myrank_ra3
          call mpi_bcast(eps, matsiz**2*nomg, mpi_double_complex,0, &
     &                  mycomm_ra3,ierr)

          if(iq.eq.1.and.iop_coul_c.eq.-1) then
            call mpi_bcast(head, nomg, mpi_double_complex,0, &
     &                   mycomm_ra3,ierr)
            call mpi_bcast(epsw1, matsiz*nomg, mpi_double_complex,0, &
     &                  mycomm_ra3,ierr)
            call mpi_bcast(epsw2, matsiz*nomg, mpi_double_complex,0, &
     &                  mycomm_ra3,ierr)
          endif
        endif
#endif
      endif 

      if(iq.eq.1.and.iop_aniso.ne.-1.and.iop_coul_c.eq.-1) then
!        call end_aniso
        deallocate(u_ani_iom)
      endif 

      !! save eps to files 
      if(myrank_ra3.eq.0.and.l_save_dielmat) then 
        call io_eps('W',iq,iom_f,iom_l,ierr)
      endif 

      call cpu_time(tend)
      time_eps = time_eps + tend - tstart

      contains 

      subroutine sub_write_emac(iop)
      integer:: iop

      if(iq.ne.1) return 
      if(metallic) then
        write(fid_outgw,'(a,f12.4)')"Plasmon frequency (eV):",sqrt(c0_head)*hev
      endif
      open(unit=99,file=trim(casename)//".emac",action='write')

      if(iop.ne.0) then 
        write(99,10)
        write(99,11)
        do iom=iom_f,iom_l
          write(99,12) iom,omega(iom)*hev,emac(1:2,iom)
        enddo
      else
        write(99,13)
        do iom=iom_f,iom_l
          write(99,14) iom,omega(iom)*hev,head(iom)
        enddo 
      endif
      close(99)
 10   format(/,"# emac with and without local field effects")
 11   format("# iom",6x,"\omega(eV)",10x," \eps_M     ",10x, &
     &                      10x," \eps_M(NLF)")
 12   format(i4,f16.6,5e16.6)
 13   format(/,"# emac without local field effects")
 14   format(i4,f16.6,2e16.6)
      end subroutine 


      subroutine sub_invbody
      
      implicit none 
      integer :: iom  ! index for mixed basis and freq.
      integer :: lwork
      integer, allocatable :: ipiv(:)
      complex(8), allocatable :: work(:)

      lwork=64*matsiz
      allocate(ipiv(matsiz),work(lwork),stat=ierr)
      call errmsg(ierr.ne.0,sname,"fail to allocate work,ipiv")

      do iom=iom_f,iom_l
        call cpu_time(time1)
        !if(iop_freq.eq.2) then  !! real freq
          call zgetrf(matsiz,matsiz,eps(:,:,iom),matsiz,ipiv,ierr)
          call errmsg0(ierr,sname,"calling zgetrf")
          call zgetri(matsiz,eps(:,:,iom),matsiz,ipiv,work,lwork,ierr)
          call errmsg0(ierr,sname,"calling zgetri")
        !else   !! imaginary freq
        !  call zhetrf('u',matsiz,eps(:,:,iom),matsiz,ipiv,work,lwork,ierr)
        !  call errmsg0(ierr,sname,"calling zhetrf")
        !  call zhetri('u',matsiz,eps(:,:,iom),matsiz,ipiv,work,ierr)
        !  call errmsg0(ierr,sname,"calling zhetri")
        !endif
        call cpu_time(time2)
        time_lapack=time_lapack+time2-time1
      enddo ! iom
      deallocate(ipiv,work)
      l_body_invd=.true.
      end subroutine sub_invbody


      subroutine sub_calcinveps
      implicit none 
      
      integer :: im,jm,iom  ! index for mixed basis and freq.
      complex(8), allocatable :: w2b(:), bw1(:)
      complex(8) :: ten_aob_tmp(3,3)
      complex(8) :: bw1_tmp, w2b_tmp ! for test use
      complex(8),external:: zdotu,ten_rvctrv
      external :: zgemv

      ! sub_invbody should be first called to make the body B `eps`
      ! actually B^{-1}
      call errmsg(.not.l_body_invd,sname,"body eps not inverted.")
      allocate(bw1(matsiz),w2b(matsiz),stat=ierr)
      call errmsg(ierr.ne.0,sname,"fail to allocate bw1,w2b")

      do iom=iom_f,iom_l
        if(iq.eq.1.and.iop_coul_c.eq.-1) then  !!  Gamma point 
          call cpu_time(time1)
          if(iop_aniso.ne.-1) then
            call cpu_time(time3)
            ! calculate vector a and vector b
            call zgemm('n','t',3,matsiz,matsiz,-cone, &
                       vec_u_ani(:,:,iom),3,eps(:,:,iom),matsiz, &
                       czero, vec_a_ani(:,:,iom), 3, ierr)
            call zgemm('n','n',3,matsiz,matsiz,-cone, &
                       vec_t_ani(:,:,iom),3,eps(:,:,iom),matsiz, &
                       czero, vec_b_ani(:,:,iom), 3, ierr)
            call cpu_time(time4)
            time_aniso = time_aniso + time4 - time3
            call zgemv('n',matsiz,matsiz,cone,eps(:,:,iom),matsiz, &
     &             epsw1(:,iom),1,czero,bw1,1,ierr)
            call zgemv('t',matsiz,matsiz,cone,eps(:,:,iom),matsiz, &
     &             epsw2(:,iom),1,czero,w2b,1,ierr)
            write(fid_outgw, "(I3,A10,6f12.3)") iom,"A before:",ten_a_ani(1,:,iom)
            write(fid_outgw, "(   A13,6f12.3)") " ",ten_a_ani(2,:,iom)
            write(fid_outgw, "(   A13,6f12.3)") " ",ten_a_ani(3,:,iom)
            write(fid_outgw, "(I3,A10,6f12.3)") iom,"P       :",ten_p_ani(1,:,iom)
            write(fid_outgw, "(   A13,6f12.3)") " ",ten_p_ani(2,:,iom)
            write(fid_outgw, "(   A13,6f12.3)") " ",ten_p_ani(3,:,iom)
            ! calculate tensor A
            ! following two zgemm should give identical result
            call cpu_time(time3)
            call zgemm('n','t',3,3,matsiz,cone, &
                       vec_t_ani(:,:,iom),3,vec_a_ani(:,:,iom),3, &
                       cone, ten_a_ani(:,:,iom), 3, ierr)
            !call zgemm('n','t',3,3,matsiz,cone, &
            !           vec_b_ani(:,:,iom),3,vec_u_ani(:,:,iom),3, &
            !           cone, ten_a_ani(:,:,iom), 3, ierr)
            !    ! results from above two ways should be identical
            !write(*,*) "ten_a_ani ZGEMM ierr = ",ierr
            call cpu_time(time4)
            time_aniso = time_aniso + time4 - time3
            write(fid_outgw,"(I3,A10,6f12.3)") iom,"A after :",ten_a_ani(1,:,iom)
            write(fid_outgw,"(   A13,6f12.3)") " ",ten_a_ani(2,:,iom)
            write(fid_outgw,"(   A13,6f12.3)") " ",ten_a_ani(3,:,iom)
          else
            call zgemv('n',matsiz,matsiz,cone,eps(:,:,iom),matsiz, &
     &               epsw1(:,iom),1,czero,bw1,1,ierr)
            call zgemv('t',matsiz,matsiz,cone,eps(:,:,iom),matsiz, &
     &               epsw2(:,iom),1,czero,w2b,1,ierr)
          endif ! iop_aniso.ne.-1

          call cpu_time(time2)
          time_lapack=time_lapack+time2-time1

          ! test for bw1 and w2b
!          if(iop_aniso.ne.-1)then
!            do imats=1,matsiz
!              bw1_tmp = bw1(imats)
!              w2b_tmp = w2b(imats)
!              bw1(imats)=sqrt(ccoefcoul)*sum(vec_a_ani(:,imats,iom)*cmplx(q0_eps(:),0.0D0,8))
!              w2b(imats)=sqrt(ccoefcoul)*sum(vec_b_ani(:,imats,iom)*cmplx(q0_eps(:),0.0D0,8))
!              if(ldbg) then
!                  write(fid_outdbg,"(i3,i4,A2,4e13.5)") &
!     &                  iom,imats,"O",bw1_tmp,w2b_tmp      
!                  write(fid_outdbg,"(i3,i4,A2,4e13.5)") &
!     &                  iom,imats,"N", bw1(imats),w2b(imats)
!                  write(fid_outdbg,"(i3,i4,A2,2L26)") &
!     &                  iom,imats,"D", &
!     &                  abs(bw1(imats)-bw1_tmp).le.1.0D-12,&
!     &                  abs(w2b(imats)-w2b_tmp).le.1.0D-12
!              endif
!            enddo
!          endif ! iop_aniso.ne.-1

          emac(2,iom)=head(iom)
          ! head calculated from the following two ways should be equivalent
          !head(iom)=1.d0/(head(iom)-zdotu(matsiz,epsw2(:,iom),1,bw1,1))
          head(iom)=1.d0/(head(iom)-zdotu(matsiz,epsw1(:,iom),1,w2b,1))

          ! test for head for anisotropy
!          head_tmp = head(iom)
!          if(iop_aniso.ne.-1)then
!            head(iom) = cone / &
!     &         (cone+ccoefcoul*ten_rvctrv(3,ten_a_ani(:,:,iom),q0_eps))
!!            head(iom)=1.d0/ & 
!!     &        (1.0D0-ccoefcoul*ten_rvctrv(3,ten_p_ani(:,:,iom),q0_eps)&
!!     &            -zdotu(matsiz,epsw1(:,iom),1,w2b,1)) ! tested correct
!            write(fid_outgw,"(A6,I3,A2,2e13.4)") "e-100",iom,"O",head_tmp
!            write(fid_outgw,"(A6,I3,A2,2e13.4)") "e-100",iom,"N",head(iom)
!            write(fid_outgw,"(A6,I3,A2,2e13.4)") "e-100",iom,"D",head(iom)-head_tmp
!          endif
          emac(1,iom)=1.d0/head(iom)

          if(iop_aniso.ne.-1)then
!            do im=1,matsiz
!              epsw1(im,iom)=-head(iom)*sqrt(ccoefcoul)*&
!     &            sum(vec_a_ani(:,im,iom)*q0_eps)
!              epsw2(im,iom)=-head(iom)*sqrt(ccoefcoul)*&
!     &            sum(vec_b_ani(:,im,iom)*q0_eps)
!            enddo
            call zgemv('T',3,matsiz,-head(iom)*sqrt(ccoefcoul),&
     &          vec_a_ani(:,:,iom),3,cmplx(q0_eps,0.0D0,8),1,&
     &          czero,epsw1(:,iom),1,ierr)
            call zgemv('T',3,matsiz,-head(iom)*sqrt(ccoefcoul),&
     &          vec_b_ani(:,:,iom),3,cmplx(q0_eps,0.0D0,8),1,&
     &          czero,epsw2(:,iom),1,ierr)
            ! TODO body of the eps^-1 could be optimized
!            do im=1,matsiz
!              do jm=1,matsiz
!!                do i=1,3
!!                  do j=1,3
!!                    ten_aob_tmp(i,j)=vec_a_ani(i,im,iom)*vec_b_ani(j,jm,iom)
!!                  enddo
!!                enddo
!!                eps(im,jm,iom)=eps(im,jm,iom) + &
!!     &              ccoefcoul*head(iom)*ten_rvctrv(3,ten_aob_tmp,q0_eps)
!                eps(im,jm,iom)=eps(im,jm,iom)+epsw1(im,iom)*epsw2(jm,iom)/head(iom)
!              enddo
!            enddo
            ! Use eps averaged over smallq region for self energy calculation
            ! Only use this in 3D.
            ! TODO different treatment for 2D/1D cases
!            if(ldbg) then
!              write(fid_outdbg,"(A30)") "V/H wing before average"
!              do im=1,matsiz
!                write(fid_outdbg,"(I3,I4,A1,2e13.4)") iom,im,"V",epsw1(im,iom)
!                write(fid_outdbg,"(I3,I4,A1,2e13.4)") iom,im,"H",epsw2(im,iom)
!              enddo
!            endif
!            call cpu_time(time3)
            call calc_h_w_inv_ang_grid(iom)
            call proj_head_on_ylm(iom)
!            call cpu_time(time4)
!            time_aniso = time_aniso + time4 - time3
!            if(ldbg) then
!              write(fid_outdbg,"(A30)") "V/H wing after average"
!              do im=1,matsiz
!                write(fid_outdbg,"(I3,I4,A1,2e13.4)") iom,im,"V",epsw1(im,iom)
!                write(fid_outdbg,"(I3,I4,A1,2e13.4)") iom,im,"H",epsw2(im,iom)
!              enddo
!            endif
          else
            epsw1(:,iom)=-head(iom)*bw1(:)
            epsw2(:,iom)=-head(iom)*w2b(:)
          endif ! iop_aniso.ne.-1

          do jm=1,matsiz
            do im=1,matsiz
              eps(im,jm,iom)=eps(im,jm,iom)+epsw1(im,iom)*epsw2(jm,iom)/head(iom)
            enddo
          enddo

          if(iop_aniso.ne.-1)then
            call angint_eps_sph(iom, eps(:,:,iom), .FALSE.)
          endif

          if(ldbg) then
            write(fid_outdbg,"(A30)") "V/H wing of eps-1"
            do im=1,matsiz
              write(fid_outdbg,"(I3,I4,A1,2e13.4)") iom,im,"V",epsw1(im,iom)
              write(fid_outdbg,"(I3,I4,A1,2e13.4)") iom,im,"H",epsw2(im,iom)
            enddo
          endif

        endif ! iq.eq.1.and.iop_coul_c.eq.-1

       ! iop == 2, inveps - 1 is calculated  
        if(iop.eq.2) then
          if(iq.eq.1.and.iop_coul_c.eq.-1) head(iom)=head(iom)-cone
          do im=1,matsiz
            eps(im,im,iom)=eps(im,im,iom)-cone
          enddo
        endif ! iop.eq.2
      enddo ! iom
      deallocate(w2b,bw1)

      l_eps_invd=.true.

      end subroutine ! sub_calcinveps

      end subroutine calceps
!EOC

