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
      use task,        only:casename,time_lapack,time_eps,l_save_dielmat
      use modmpi

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

      integer :: iat,icg,ic          ! indices for core states of a given atom.
      integer :: ie1                 ! index occupied states   
      integer :: ie2                 ! index for unoccupied states 
      integer :: nktot,iik,ik,irk,jk,jrk   ! indices for the k-points.
      integer :: iom                 ! Counter: Runs over frequencies.
      integer :: isp                 ! Counter: runs over spin
      integer :: ie12
      integer :: ierr
      integer :: im  ! index for basis set 
      integer :: nvbm,ncbm,nmdim,cmdim
      integer :: ie2_f,ie2_l  
      integer :: ie1_f,ie1_l  
      integer :: ik_f, ik_l 
      
      logical :: ldbg=.false.

      integer :: nblk,iblk,m_f,m_l
      integer :: nomg
      integer :: iop_minm

      real(8) :: edif,coefw,aux,ffact,omgsq

      complex(8):: coef
      real(8) :: time1,time2,tstart,tend
      character(len=15)::sname='calceps'

      real(8),   allocatable :: enk(:)        !! local array for eigen-energies 
      complex(8),allocatable :: tmat(:,:),pm(:),wtmp(:),minm(:,:,:) 
      character(len=10),external:: int2str 

      real(8),parameter :: sqr3=1.732050807568877294d0

 
!
!EOP
!BOC

      if(ldbg) call linmsg(6,'-',sname)
      call cpu_time(tstart)

      if(lread) then 
        call io_eps('r',iq,iom_f,iom_l,ierr)
        if(ierr.eq.0) then 
          write(6,*) "Restart mode: read eps from files"
          return 
        else 
          write(6,*) "WARNING: fail to read eps from files"
        endif 
      endif 
!
!     isc.eq.0 indicates either it is a G0W0 calculation or it is
!     the first iteration of QSGW calculations, one just need to reads
!     from Minm file directly, in the case of isc.gt.0, indicating
!     a QSGW calculation, Minm are read from the the file and transformed 
!     in terms of both "n" and "m" by qpwf_coef 
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

      ! first calculate momentum matrix and the head
      if(iq.eq.1.and.iop_coul_c.eq.0) then
        call init_mommat(1,nomax,numin,nbmaxpol,nirkp,nspin)

        if(iop_mommat.eq.0) then
          call calcmommat(0,1,nomax,numin,nbmaxpol,0)
        else
          call w2k_readmommat(1,nomax,numin,nbmaxpol)
        endif
        call calchead(1,nirkp,iom_f,iom_l)
      endif

      ! set the parallelization over k (within the row)
      call mpi_set_range(nproc_row,myrank_row,nktot,1,ik_f,ik_l)
      write(6,*) "calceps: myrank_row,ik_f,ik_l =",myrank_row,ik_f,ik_l

      coefw = 2.0d0*sqrt(pi*vi)
      do isp=1,nspin  
        do iik=ik_f, ik_l 

          ! set ik, irk, jk and jrk 
          if(isym.eq.0) then 
            ik=iik 
            irk=kpirind(ik)
            coef= - cone*fspin 
          else 
            irk=iik
            ik=idikp(irk)
            coef= - wkir(irk)*cone*fspin
          endif 

          jk = kqid(ik,iq)
          jrk= kpirind(jk)
          nvbm = nomaxs(isp) 
          ncbm = numins(isp) 

          !! set the local array for band energies  
          allocate(enk(ncg_p + nbmaxpol))
          do icg=1,ncg_p       
            iat = corind(1,icg)
            ic  = corind(3,icg)
            enk(icg) = eigcore(ic,iat,isp) 
          enddo 
          enk(ncg_p+1:ncg_p+nbmaxpol) = bande(1:nbmaxpol,irk,isp) 

          !! find the index for the band whose energy is higher than eminpol  
          ie1_l = nvbm + ncg_p
          do ie1=1,ie1_l
            ie1_f = ie1
            if(enk(ie1) > eminpol) exit 
          enddo  
          if(ldbg) write(6,*) " ie1_f=",ie1_f 

          !! set mask_eps that defines the transitions to be skipped 
          call crpa_setmask(irk,jrk,isp) 

          !!  determine the number of blocks used in minm operation 
          nblk=max((nbmaxpol-ncbm+1)/mblksiz,1) 

          if(ldbg) write(6,*) "nblk=",nblk 

          do iblk=1,nblk                     !! loop over m-blocks in Minm
            m_f = ncbm+(iblk-1)*mblksiz
            m_l = m_f + mblksiz-1
            if(iblk.eq.nblk) m_l = nbmaxpol

            call mpi_set_range(nproc_3rd,myrank_3rd,m_l-m_f+1,m_f,ie2_f,ie2_l)
  
            nmdim = (ie2_l-ie2_f+1)*(ie1_l-ie1_f+1)  

            allocate(tmat(1:matsiz,nmdim),& 
     &               pm(nmdim),& 
     &               wtmp(matsiz),&
     &               minm(matsiz,ie2_f:ie2_l,ie1_f:ie1_l), & 
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
            if(iq.eq.1.and.iop_coul_c.eq.0) then
              ie12=0
              do ie1=ie1_f,ie1_l       ! the index for occ. states
                do ie2=ie2_f,ie2_l     ! the index for unocc. states. 
                  ie12=ie12+1
                  edif = enk(ie2+ncg_p) - enk(ie1)

                  if(abs(edif).gt.1.0d-10)then
                    aux=coefw/edif
                  else 
                    aux=0.d0
                  endif 

                  if(ie1.le.ncg_p) then 
                    pm(ie12)=aux*sum(mmatcv(:,ie1,ie2,irk,isp)*q0_eps)
                  else 
                    pm(ie12)=aux*sum(mmatvv(:,ie1-ncg_p,ie2,irk,isp)*q0_eps)
                  endif 
                enddo 
              enddo 
            endif   ! iq.eq.1

            do iom=iom_f,iom_l
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
              call zgemm('n','c',matsiz,matsiz,nmdim,coef,tmat,matsiz, &
     &            minm,matsiz,cone,eps(:,:,iom),matsiz)   

              if(iq.eq.1.and.iop_coul_c.eq.0) then 
                call zgemv('n',matsiz,nmdim,coef,tmat,matsiz,pm,1,czero,wtmp,1)
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
                  call zgemv('n',matsiz,nmdim,coef,tmat,matsiz,pm,1,czero,wtmp,1)
                endif
                epsw2(:,iom)=epsw2(:,iom)+conjg(wtmp)
              endif

              call cpu_time(time2)
              time_lapack=time_lapack+time2-time1
            enddo ! iom
            deallocate(minm,tmat,pm,wtmp) 

          enddo ! iblk  
          deallocate(enk) 

        enddo ! iik
      enddo ! isp

#ifdef MPI 
      if(nproc_ra3.gt.1) then 
        write(6,*) "Collect eps: myrank_ra3=",myrank_ra3,"comm=",mycomm_ra3
        call mpi_sum_array(0,eps,matsiz,matsiz,nomg,mycomm_ra3)
        if(iq.eq.1) then
          call mpi_sum_array(0,epsw1,matsiz,nomg,mycomm_ra3)
          call mpi_sum_array(0,epsw2,matsiz,nomg,mycomm_ra3)
        endif
      endif
#endif

      !! add "1" to the diagonal elements
      if(myrank_ra3.eq.0) then 
        do iom=iom_f,iom_l
          do im=1,matsiz
            eps(im,im,iom) = eps(im,im,iom) + cone
          enddo
         enddo
      endif

      if(iq.eq.1.and.iop_coul_c.eq.0) then
        call end_mommat
      endif 

      if(iop.eq.0) then  !! calculate only eps 
        call sub_write_emac(0)  
      else               !!  calculate inveps 
        if(myrank_ra3.eq.0) then
          call sub_calcinveps 
          call sub_write_emac(1) 
        endif 

#ifdef MPI 
        if(nproc_ra3.gt.1) then
          write(6,*) "bcast eps: myrank_ra3=",myrank_ra3
          call mpi_bcast(eps, matsiz**2*nomg, mpi_double_complex,0, &
     &                  mycomm_ra3,ierr)

          if(iq.eq.1.and.iop_coul_c.eq.0) then
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
        write(6,'(a,f12.4)')"Plasmon frequency (eV):",sqrt(c0_head)*hev
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

      subroutine sub_calcinveps
      implicit none 
      
      integer :: im,jm,iom  ! index for mixed basis and freq.
      integer :: lwork

      complex(8), allocatable :: w2b(:), bw1(:)
      integer, allocatable :: ipiv(:)
      complex(8), allocatable :: work(:)

      external zhemm
      complex(8),external::zdotc,zdotu

      lwork=64*matsiz
      allocate(ipiv(matsiz),work(lwork),bw1(matsiz),w2b(matsiz),stat=ierr)
      call errmsg(ierr.ne.0,sname,"fail to allocate work,ipiv")

      do iom=iom_f,iom_l

        call cpu_time(time1)
        if(iop_freq.eq.2) then  !! real freq
          call zgetrf(matsiz,matsiz,eps(:,:,iom),matsiz,ipiv,ierr)
          call errmsg0(ierr,sname,"calling zgetrf")
          call zgetri(matsiz,eps(:,:,iom),matsiz,ipiv,work,lwork,ierr)
          call errmsg0(ierr,sname,"calling zgetri")
        else   !! imaginary freq
          call zhetrf('u',matsiz,eps(:,:,iom),matsiz,ipiv,work,lwork,ierr)
          call errmsg0(ierr,sname,"calling zhetrf")
          call zhetri('u',matsiz,eps(:,:,iom),matsiz,ipiv,work,ierr)
          call errmsg0(ierr,sname,"calling zhetri")
        endif
        call cpu_time(time2)
        time_lapack=time_lapack+time2-time1

        if(iq.eq.1.and.iop_coul_c.eq.0) then  !!  Gamma point 
          call cpu_time(time1)

          if(iop_freq.eq.2) then
            call zgemv('n',matsiz,matsiz,cone,eps(:,:,iom),matsiz, &
     &              epsw1(:,iom),1,czero,bw1,1)
            call zgemv('t',matsiz,matsiz,cone,eps(:,:,iom),matsiz, &
     &              epsw2(:,iom),1,czero,w2b,1)
          else
            call zhemv('u',matsiz,cone,eps(:,:,iom),matsiz, &
     &              epsw1(:,iom),1,czero,bw1,1)
            call zhemv('u',matsiz,cone,eps(:,:,iom),matsiz, &
     &              epsw2(:,iom),1,czero,w2b,1)
            w2b=conjg(w2b)
          endif

          call cpu_time(time2)
          time_lapack=time_lapack+time2-time1

          emac(2,iom)=head(iom)
          head(iom)=1.d0/(head(iom)-zdotu(matsiz,epsw2(:,iom),1,bw1,1))
          emac(1,iom)=1.d0/head(iom)

          epsw1(:,iom)= - head(iom)*bw1(:)
          epsw2(:,iom)= - head(iom)*w2b(:)
          do jm=1,matsiz
            do im=1,matsiz
              eps(im,jm,iom)=eps(im,jm,iom)+epsw1(im,iom)*epsw2(jm,iom)/head(iom)
            enddo
          enddo
        endif

       ! iop == 2, inveps-1 is calculated  
        if(iop.eq.2) then
          if(iq.eq.1.and.iop_coul_c.eq.0) head(iom)=head(iom)-cone
          do im=1,matsiz
            eps(im,im,iom)=eps(im,im,iom)-cone
          enddo
        endif
      enddo ! iom
      deallocate(ipiv,work,w2b,bw1)

      end subroutine ! sub_calcinveps


      end subroutine calceps
!EOC

