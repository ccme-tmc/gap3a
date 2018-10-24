!-----------------------------------------------------------------------
!BOP
!
! !MODULE: consrpa
!  a module for the constrained RPA interaction 
!
      module crpa
      implicit none 

      integer:: nsp_crpa       ! the number of spin-components of CRPA matrices
                               !   1 -- for nspin = 1 
                               !   3 -- for nspin = 2, i.e. (u,u),(d,d),(u,d)
                                
      integer:: iop_crpa       ! indicate which cRPA scheme 
                               ! -1 -- calculate the matrix of the bare Coulomb interaction (v)
                               !  0 -- calculate the matrix of the screened Coulomb interaction (W)
                               !  1 -- cRPA U with the mask approach 
                               !  2 -- cRPA U with the weighted mask approach 
                               !        (as in Sasioglu et al. PRB 83, 121101(R)(2011))
                               !  3 -- the projection approach  (to do)
                               !  4 -- the impurity U approach  (to do) 
                               !  5 -- the disentanglement approach 
                               !       by Miyake et al.(PRB 80, 155134 (2009) (not implemented) 
      integer:: iop_wf         ! option about the Wannier function 
                               !  0 -- use the internally generated WF (not implemented yet!)
                               !  1 -- use WF generated by dmftproj 
                               !  2 -- use WF generated by wien2wannier + wannier90 (not implemented yet!) 
      integer:: iop_intercell  ! the option to control how inter-cell interactions are considered
                               !  0 -- no intercell interactions
                               ! -1 -- the neareast neighbouring(NN): (1,0,0),(0,1,0),(0,0,1) 
                               ! -2 -- including the second NN: ! (1,1,0),(
                               ! >1 -- read the cell information from gw.inp 

      integer::iop_pln_phase=0     ! if equal to 1, an additional phase factor is added to pln(ilm,ie,ik,isp)  

      integer:: nlorb_wf       ! number of correlated orbitals (counting a,l only [ a -- index for atoms]
      integer:: nlmorb         ! total number of correlated states (countering a,l,m)
                               !      nlmorb = \sum_i nmorb(i)
      integer,allocatable:: info_orb(:,:)   ! the information about each set of (a,l) orbitals 
                               ! (1,:) index of the atom that (a,l) belong to
                               ! (2,:) l
                               ! (3,:) number of different "m"-state for each (a,l)
                               ! (4,:) the type of (a,l)-orbitals 
      integer,allocatable :: info_wf(:,:) 

      integer:: nsrt           ! number of sort of atoms 
      integer,allocatable :: nmult(:)          ! multiplicity of atoms 
      integer,allocatable :: natcorr(:)        ! number of correlated atoms

      integer:: nbmax_wf,nbmin_wf   ! minimal/maximal band index included in the construction of Wannier functions  
      integer:: nbands_wf           ! the number of bands falling withiin the Wannier window

      real:: tol_shift_centers = 0.01 
      real(8),allocatable :: wf_centers(:,:)       ! the centers of Wannier functions in the unit cell 
      real(8),allocatable :: wf_centers_new(:,:)   ! the new centers of Wannier functions to be used  
      real(8),allocatable :: wf_centers_lm(:,:)  ! the centers of Wannier functions in the unit cell 
      real(8),allocatable :: weight_corr(:,:,:) ! the weight of correlated states in each KS band
      complex(8),allocatable:: cproj_corr(:,:,:,:) ! the weight of correlated states in each KS band
      integer,allocatable :: nbk_wf(:,:,:)   ! windows for bands to be included in cRPA calculations 
      complex(8),allocatable:: pln(:,:,:,:)    ! the projection of a Kohn-Sham vector to a Wannier function (wf)

      integer              :: ncell=1        ! number of unit cells for calc U_{0R0R} (=1 if only on-site interaction is considered)
      real(8), allocatable :: rcell(:,:)     ! cell coordinate in the internal basis (in the wien2k convention)
      real(8), allocatable :: rcell_cart(:,:)! cell coordinate in cartesian basis

      complex(8),allocatable:: mill(:,:,:)     
      complex(8),allocatable:: umat(:,:,:,:,:)                ! matrix elements of W_{r} in terms of Wannier functions
      complex(8),allocatable:: vmat(:,:,:,:)                  ! matrix elements of W_{r} in terms of Wannier functions

      character(len=20),private:: sname="module:crpa"

      contains

        subroutine init_crpa(iq)
        use bands,    only: nspin 
        use freq,     only: nomeg
        use kpoints,  only: nkp
        use mixbasis, only: matsiz
        implicit none 
        integer,intent(in):: iq
        integer:: ierr
        integer:: isp,ik,ib,jb,ilm,jlm
        integer:: nbot,ntop
        integer:: lldim
        logical:: ldbg = .false.

        nsp_crpa = nspin*(nspin+1)/2 
        lldim = nlmorb*nlmorb

        if (iq.eq.0) then 
          allocate(umat(lldim,lldim,nsp_crpa,nomeg,ncell),   &
     &             vmat(lldim,lldim,nsp_crpa,ncell), &
     &             stat=ierr)
          umat = 0.d0 
          vmat = 0.d0
          call errmsg(ierr.ne.0,sname,"Fail to allocate umat ... ")

          allocate(weight_corr(nbmin_wf:nbmax_wf,nkp,nspin),  &
     &     cproj_corr(nbmin_wf:nbmax_wf,nbmin_wf:nbmax_wf,nkp,nspin),  &
     &     stat=ierr)
          call errmsg(ierr.ne.0,sname,"Fail to allocate cproj")

          weight_corr = 0.d0
          cproj_corr = 0.d0

          do isp=1,nspin
            do ik=1,nkp
              nbot = nbk_wf(1,ik,isp)
              ntop = nbk_wf(2,ik,isp)
              do jb = nbot,ntop
                do ilm=1,nlmorb
                  weight_corr(jb,ik,isp) = weight_corr(jb,ik,isp) &
     &             +abs(pln(ilm,jb,ik,isp))**2
                enddo

                do ib = nbot,ntop
                  do ilm=1,nlmorb
                    cproj_corr(ib,jb,ik,isp) = cproj_corr(ib,jb,ik,isp) &
     &               + conjg(pln(ilm,ib,ik,isp))*pln(ilm,jb,ik,isp)
                  enddo
                enddo
              enddo
             
              if(ldbg) then 
                write(6,*)
                write(6,101) "k=",ik
                write(6,100) "Total d-weights in each state"
                write(6,201) "ib","weight"
                do ib = nbot,ntop
                  write(6,202) ib,weight_corr(ib,ik,isp)
                enddo

                write(6,100) "The d-projector in the Kohn-Sham space"
                write(6,203) "n","n'","C_{n,n'}"
                do jb=nbot,ntop
                  do ib=nbot,ntop
                    write(6,204) ib,jb,cproj_corr(ib,jb,ik,isp)
                  enddo
                enddo
                write(6,*)
              endif 
            enddo  ! ik
          enddo ! isp 
        else  ! iq > 0 
          allocate(mill(matsiz,lldim,nspin),stat=ierr)
          call errmsg(ierr.ne.0,sname,"Fail to allocate mill ")
        endif 
 100  format(a)
 101  format(A,i5)
 102  format(A,2i5)
 105  format(A,5i5)
 201  format(a5,a8)
 202  format(i5,f8.3)
 203  format(a5,a5,a16)
 204  format(i5,i5,2f8.3)

      end subroutine 

      subroutine end_crpa(iq)
        implicit none 
        integer,intent(in):: iq
        if (iq.eq.0) then 
          deallocate(umat,vmat,weight_corr,cproj_corr)
        else
          deallocate(mill)
        endif 
      end subroutine 

      end module 
!EOP
