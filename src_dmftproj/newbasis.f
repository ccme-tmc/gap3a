      SUBROUTINE newbasis(tetr)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine computes the new basis which diagonalizes the   %%
C %% density matrices for all the correlated orbitals in the         %%
C %% the energy window for which the projectors were previously      %%
C %% defined.                                                        %%
C %%                                                                 %%
C %% If tetr = .TRUE., the tetrahedron weights are used.             %%
C %%         = .FALSE., a simple point integration is used.          %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
      USE almblm_data
      USE common_data
      USE prnt
      USE projections
      USE reps
      USE symm
      IMPLICIT NONE
      LOGICAL :: tetr, ifparabis
      COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: D, cmat 
      COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zeros,tmat
      COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: conj_mat, mat
      COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:) :: work     
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: eigen, rwork
      INTEGER :: is, is1, ik, iorb, l, nbnd, iss 
      INTEGER :: lm, i, icrorb, m1, m2, m, ib 
      INTEGER :: iatom, jatom, jorb, imu, isrt, isym
      INTEGER :: istart, istart1
      INTEGER :: nbbot, nbtop, info
      REAL(KIND=8) :: q, qtot
C
        WRITE(buf,*)'Evaluation of the basis that diagonalizes' 
     &                   ,' the density matrix...'
        CALL printout(0)
        CALL printout(0)

        OPEN(25,file='Diagonal_density.dat',action='write')
C
C Initialization of the variable nsp and of the table crdensmat
C nsp describes the number of block necessary to describe the complete density matrix.
        nsp=ns
        IF(ifSO) nsp=4
C The only possible cases are therefore :
C nsp=1 for a computation without spin-polarized input files, without SO   [block up/up]
C nsp=2 for a computation with spin-polarized input files (but without SO) [blocks up/up and dn/dn]
C nsp=4 for a computation (with spin-polarized input files) with SO        [all the blocks]
C
C
C =============================================================
C Computation of the density matrices for correlated orbitals :
C =============================================================
C
C These computations are performed only if there are correlated orbitals in the system.
        IF(ncrorb.NE.0) THEN
C crdensmat is a table of size nsp*ncrorb. 
C For each correlated orbital icrorb, crdensmat(:,icrorb) is the corresponding density matrix.
        IF(.NOT.ALLOCATED(crdensmat)) THEN
         ALLOCATE(crdensmat(nsp,ncrorb))
        ENDIF
        DO icrorb=1,ncrorb
          DO is=1,nsp
            IF(ALLOCATED(crdensmat(is,icrorb)%mat)) 
     &       DEALLOCATE(crdensmat(is,icrorb)%mat)
            ALLOCATE(crdensmat(is,icrorb)%mat(1,1))
            crdensmat(is,icrorb)%mat(1,1)=0.d0
          ENDDO
        ENDDO
C 
C Loop on the correlated orbitals icrorb
C
        DO icrorb=1,ncrorb
          isrt=crorb(icrorb)%sort
          l=crorb(icrorb)%l
C
C -----------------------------------------------------------------------------------
C The s-orbitals are a particular case of a "non-mixing" basis and are treated here :
C -----------------------------------------------------------------------------------
          IF (l==0) THEN
C The field mat of crdensmat has already the good size (1 scalar element). 
C There's no need of a basis change since the representation of an s-orbital is Identity whatever the basis is.
           DO ik=1,nk
             DO iss=1,nsp
C
C Determination of the block indices :
C ------------------------------------
               IF(iss.LE.2) THEN
                is=iss
                is1=iss 
C If iss=1 (up), is=1 and is1=1 -> Description of the up/up block
C If iss=2 (down), is=2 and is1=2 -> Description of the dn/dn block
               ELSE 
                is=iss-2 
                is1=3-is
C If iss=3, is=1 (up) and is1=2 (down) -> Description of the up/dn block
C If iss=4, is=2 (down) and is1=1 (up) -> Description of the dn/up block
               ENDIF
C Only the kpoints with inlcuded bands are considered for the projectors.
               IF(.NOT.kp(ik,is)%included) CYCLE
               IF(.NOT.kp(ik,is1)%included) CYCLE
               nbnd=kp(ik,is)%nb_top-kp(ik,is)%nb_bot+1
               nbbot=kp(ik,is)%nb_bot
               nbtop=kp(ik,is)%nb_top
C for the diagonal blocks, is=is1 ; thus nbtop and nbtop are the same for is and is1. 
C for the off-diagonal blocks (calculated only when SO is considered), 
C it was checked that nbtop(up)=nbtop(dn) and nbbot(up)=nbbot(dn) [in set_projections.f]
C thus nbtop and nbtop fit again for is and is1.
C
C Computation of the density matrix using the tetrahedron weights for the integration :
C -------------------------------------------------------------------------------------
               IF(tetr) THEN
C The field pr_crorb%mat_rep is used to perform the computation (well defined for s-orbitals)  
                ALLOCATE(mat(1,nbbot:nbtop),conj_mat(1,nbbot:nbtop))
                mat(1,nbbot:nbtop)=
     &            pr_crorb(icrorb,ik,is)%mat_rep(1,nbbot:nbtop)*
     &            SQRT(kp(ik,is)%tetrweight(nbbot:nbtop))
                conj_mat(1,nbbot:nbtop)=CONJG( 
     &            pr_crorb(icrorb,ik,is1)%mat_rep(1,nbbot:nbtop))*
     &            SQRT(kp(ik,is1)%tetrweight(nbbot:nbtop))
C mat = P(icrorb,ik,is)*sqrt(tetrahedron-weight(ik,is))
C conj_mat = conjugate[ P(icrorb,ik,is1)) ]*sqrt(tetrahedron-weight(ik,is1))
                DO ib = nbbot,nbtop
                  crdensmat(iss,icrorb)%mat(1,1)=
     =              crdensmat(iss,icrorb)%mat(1,1)+
     +              mat(1,ib)*conj_mat(1,ib)
                ENDDO
C crdensmat = mat*transpose(conj_mat) which is a matrix of size 1
C The summation over the k points is done with the do loop ; 
C crdensmat(iss,icrorb) is therefore the block "iss" of the density matrix for the orbital icrorb.
                DEALLOCATE(conj_mat,mat)
C
C Computation of the density matrix using a simple point integration :
C --------------------------------------------------------------------
               ELSE
                DO ib = nbbot, nbtop
                  crdensmat(iss,icrorb)%mat(1,1)=
     =              crdensmat(iss,icrorb)%mat(1,1)+
     +              pr_crorb(icrorb,ik,is)%mat_rep(1,ib)*
     *              CONJG(pr_crorb(icrorb,ik,is1)%mat_rep(1,ib))*
     *              kp(ik,is)%weight
                ENDDO
C crdensmat = P(icrorb,ik,is)*transpose(conjugate(P(icrorb,ik,is1))) which is a matrix of size 1.
C The weight used is a geometric factor associated to k and does not depend on the variable is. 
C That's why we merely multiply by the "weight" each term while summing over the k points.
               ENDIF
             ENDDO    ! End of the iss loop  
           ENDDO      ! End of the ik loop
C
C -----------------------------------------------------------------------------------------------------
C If the basis representation needs a complete spinor rotation approach (matrices of size 2*(2*l+1) ) :
C -----------------------------------------------------------------------------------------------------
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C If this option is used, then ifSO=.TRUE. (restriction in set_ang_trans.f)
C Moreover ifSP=.TRUE. (since ifSO => ifSP) 
C As a result, we know that nb_bot(up)=nb_bot(dn) and nb_top(up)=nb_top(dn)
C
C The field mat of crdensmat must be resized.
C As the complete spinor rotation approach is used, only one matrix is necessary (with is=1).
           DEALLOCATE(crdensmat(1,icrorb)%mat)
           ALLOCATE(crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1)))
           crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))=0.d0
           ALLOCATE(D(1:2*(2*l+1),1:2*(2*l+1)))
C
C Computation of the density matrix using the tetrahedron weights for the integration :
C -------------------------------------------------------------------------------------
           IF(tetr) THEN
            DO ik=1,nk
C Only the kpoints with inlcuded bands are considered for the projectors.
              IF(.NOT.kp(ik,1)%included) CYCLE
              nbnd=kp(ik,1)%nb_top-kp(ik,1)%nb_bot+1
              nbbot=kp(ik,1)%nb_bot
              nbtop=kp(ik,1)%nb_top
C A distinction between up and dn states is necessary in order to use the tetrahedron weight.
C As a result, we will transform the projectors back in spherical harmonics basis 
C to perform the calculation and then put the resulting density matrix in the desired basis.
C
              ALLOCATE(mat(1:2*(2*l+1),nbbot:nbtop))
              ALLOCATE(conj_mat(1:2*(2*l+1),nbbot:nbtop))
C The representation of the projectors is put back in the spherical harmonics basis
              mat(1:2*(2*l+1),nbbot:nbtop)=MATMUL(TRANSPOSE(CONJG( 
     =          reptrans(l,isrt)%transmat
     &          (1:2*(2*l+1),1:2*(2*l+1)))),pr_crorb(icrorb,ik,1)%
     &          mat_rep(1:2*(2*l+1),nbbot:nbtop))
C mat = inverse(reptrans)*pr_crorb%mat_rep = <lm|new_i>*proj_{new_i} = proj_{lm} [temporarly]
              conj_mat(1:2*(2*l+1),nbbot:nbtop)=MATMUL( 
     &          TRANSPOSE(CONJG(reptrans(l,isrt)%
     &          transmat(1:2*(2*l+1),1:2*(2*l+1)) )),
     &          pr_crorb(icrorb,ik,1)%mat_rep
     &          (1:2*(2*l+1),nbbot:nbtop))
C conj_mat = inverse(reptrans)*pr_crorb%mat_rep = <lm|new_i>*proj_{new_i} = proj_{lm} [temporarly]
              DO m=1,2*l+1
                mat(m,nbbot:nbtop)=mat(m,nbbot:nbtop)*
     &            SQRT(kp(ik,1)%tetrweight(nbbot:nbtop))
C The first (2*l+1) lines are associated to up states. Hence a multiplication by tetrweight(up)
                mat((2*l+1)+m,nbbot:nbtop)=
     =            mat((2*l+1)+m,nbbot:nbtop)*
     &            SQRT(kp(ik,2)%tetrweight(nbbot:nbtop))
C The last (2*l+1) lines are associated to dn states. Hence a multiplication by tetrweight(dn)
C mat = P(icrorb,ik,is)*sqrt(tetrahedron-weight(ik,is))
                conj_mat(m,nbbot:nbtop)=
     =            CONJG(conj_mat(m,nbbot:nbtop))*
     &            SQRT(kp(ik,1)%tetrweight(nbbot:nbtop))
C The first (2*l+1) lines are associated to up states. Hence a multiplication by tetrweight(up)
                conj_mat((2*l+1)+m,nbbot:nbtop)=
     &            CONJG(conj_mat((2*l+1)+m,nbbot:nbtop))*
     &            SQRT(kp(ik,2)%tetrweight(nbbot:nbtop))
C The last (2*l+1) lines are associated to dn states. Hence a multiplication by tetrweight(dn)
C conj_mat = conjugate[ P(icrorb,ik,is1))*sqrt(tetrahedron-weight(ik,is1)) ]
              ENDDO
              CALL zgemm('N','T',2*(2*l+1),2*(2*l+1),nbnd,
     &          DCMPLX(1.D0,0.D0),mat,2*(2*l+1),conj_mat,2*(2*l+1),
     &          DCMPLX(0.D0,0.D0),D,2*(2*l+1))
              DEALLOCATE(conj_mat,mat)
C D = mat*transpose(conj_mat) is a matrix of size 2*(2*l+1)* 2*(2*l+1)
C
              crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     &          crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))
     &          +D(1:2*(2*l+1),1:2*(2*l+1))
            END DO       ! End of the ik loop 
C The summation over the k points is done.
C crdensmat(icrorb) is therefore the complete density matrix in spherical harmonic basis.
C
C The density matrix is then put into the desired basis, using reptrans(l,isrt)%transmat 
            crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     =        MATMUL(reptrans(l,isrt)%transmat
     &        (1:2*(2*l+1),1:2*(2*l+1)),crdensmat(1,icrorb)%mat
     &        (1:2*(2*l+1),1:2*(2*l+1)) )
            crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     =        MATMUL(crdensmat(1,icrorb)%mat
     &        (1:2*(2*l+1),1:2*(2*l+1)),TRANSPOSE( CONJG(
     &        reptrans(l,isrt)%transmat(1:2*(2*l+1),1:2*(2*l+1)))))
C crdensmat = (reptrans)*crdensmat_l*inverse(reptrans) 
C or crdensmat = <new_i|lm> crdensmat_{lm} <lm|new_i>
C crdensmat(icrorb) is now the complete density matrix in the desired basis.
C
C Computation of the density matrix using a simple point integration :
C --------------------------------------------------------------------
           ELSE
C No distinction between up and dn states is necessary because we use a
C geometric factor which depends only of the k-point.
C We can then use directly the projectors in the desired basis.
            DO ik=1,nk
C Only the kpoints with inlcuded bands are considered for the projectors.
              IF(.NOT.kp(ik,1)%included) CYCLE
              nbnd=kp(ik,1)%nb_top-kp(ik,1)%nb_bot+1
              nbbot=kp(ik,1)%nb_bot
              nbtop=kp(ik,1)%nb_top
              ALLOCATE(mat(1:2*(2*l+1),nbbot:nbtop))
              mat(1:2*(2*l+1),nbbot:nbtop)=
     =          pr_crorb(icrorb,ik,1)%mat_rep(1:2*(2*l+1),
     &          nbbot:nbtop)
              CALL zgemm('N','C',2*(2*l+1),2*(2*l+1),nbnd,
     &          DCMPLX(1.D0,0.D0),mat,2*(2*l+1),mat,2*(2*l+1),
     &          DCMPLX(0.D0,0.D0),D,2*(2*l+1))
              DEALLOCATE(mat)
C D = P(icrorb,ik,is)*transpose(conjugate(P(icrorb,ik,is))) is a matrix of size 2*(2*l+1) * 2*(2*l+1)
              crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))= 
     =          crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))
     +          +D(1:2*(2*l+1),1:2*(2*l+1))*kp(ik,1)%weight
            ENDDO      ! End of the ik loop
C The summation over the k  points is done in the do loop. 
C The weight used is a geometric factor associated to k and does not depend on the variable is. 
C That's why we merely multiply by the "weight" each term while summing over the k points.
           ENDIF
           DEALLOCATE(D)
C
C ----------------------------------------------------------------------------------------------
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only) :
C ----------------------------------------------------------------------------------------------
          ELSE
C The field mat of crdensmat must be resized.
C nsp matrices of size (2*l+1) are necessary to represent the whole density matrix.
           DO is=1,nsp
             DEALLOCATE(crdensmat(is,icrorb)%mat)
             ALLOCATE(crdensmat(is,icrorb)%mat(-l:l,-l:l))
             crdensmat(is,icrorb)%mat(-l:l,-l:l)=0.d0
           ENDDO
           ALLOCATE(D(-l:l,-l:l))
C All the computations can be performed in the new basis (using the field pr_crorb%mat_rep)
           DO ik=1,nk
             DO iss=1,nsp
C
C Determination of the block indices :
C ------------------------------------
               IF(iss.LE.2) THEN
                is=iss
                is1=iss 
C If iss=1 (up), is=1 and is1=1 -> Description of the up/up block
C If iss=2 (down), is=2 and is1=2 -> Description of the dn/dn block
               ELSE 
                is=iss-2 
                is1=3-is
C If iss=3, is=1 (up) and is1=2 (down) -> Description of the up/dn block
C If iss=4, is=2 (down) and is1=1 (up) -> Description of the dn/up block
               ENDIF
C Only the kpoints with inlcuded bands are considered for the projectors.
               IF(.NOT.kp(ik,is)%included) CYCLE
               IF(.NOT.kp(ik,is1)%included) CYCLE
               nbnd=kp(ik,is)%nb_top-kp(ik,is)%nb_bot+1
               nbbot=kp(ik,is)%nb_bot
               nbtop=kp(ik,is)%nb_top
C for the diagonal blocks, is=is1 ; thus nbtop and nbtop are the same for is and is1. 
C for the off-diagonal blocks (calculated only when SO is considered), 
C it was checked that nbtop(up)=nbtop(dn) and nbbot(up)=nbbot(dn) [in set_projections.f]
C thus nbtop and nbtop fit again for is and is1.
C
C Computation of the density matrix using the tetrahedron weights for the integration :
C -------------------------------------------------------------------------------------
               IF(tetr) THEN
C The representation of the projectors in the desired basis is used (field pr_crorb%mat_rep)
                ALLOCATE(mat(-l:l,nbbot:nbtop))
                ALLOCATE(conj_mat(-l:l,nbbot:nbtop))
                DO m=-l,l
                  mat(m,nbbot:nbtop)=
     =              pr_crorb(icrorb,ik,is)%mat_rep(m,nbbot:nbtop)*
     &              SQRT(kp(ik,is)%tetrweight(nbbot:nbtop))
                  conj_mat(m,nbbot:nbtop)=CONJG(
     &              pr_crorb(icrorb,ik,is1)%mat_rep(m,nbbot:nbtop))
     &              *SQRT(kp(ik,is1)%tetrweight(nbbot:nbtop))
C mat = P(icrorb,ik,is)*sqrt(tetrahedron-weight(ik,is))
C conj_mat = conjugate[ P(icrorb,ik,is1))*sqrt(tetrahedron-weight(ik,is1)) ]
                ENDDO
                CALL zgemm('N','T',(2*l+1),(2*l+1),nbnd,
     &            DCMPLX(1.D0,0.D0),mat,(2*l+1),conj_mat,(2*l+1),
     &            DCMPLX(0.D0,0.D0),D(-l:l,-l:l),(2*l+1))
                DEALLOCATE(conj_mat,mat)
C D = mat*transpose(conj_mat) is a matrix of size (2*l+1)*(2*l+1)
                crdensmat(iss,icrorb)%mat(-l:l,-l:l)=
     =            crdensmat(iss,icrorb)%mat(-l:l,-l:l)+D(-l:l,-l:l)
C The summation over the k points is done.
C crdensmat(icrorb) is therefore the complete density matrix of the orbital icrorb.
C
C Computation of the density matrix using a simple point integration :
C --------------------------------------------------------------------
               ELSE
                ALLOCATE(mat(-l:l,nbbot:nbtop))
                ALLOCATE(conj_mat(-l:l,nbbot:nbtop))
                mat(-l:l,nbbot:nbtop)=pr_crorb(icrorb,ik,is)
     &            %mat_rep(-l:l,nbbot:nbtop)
                conj_mat(-l:l,nbbot:nbtop)=pr_crorb(icrorb,ik,is1)
     &            %mat_rep(-l:l,nbbot:nbtop)
                CALL zgemm('N','C',(2*l+1),(2*l+1),nbnd,
     &            DCMPLX(1.D0,0.D0),mat,(2*l+1),conj_mat,(2*l+1),
     &            DCMPLX(0.D0,0.D0),D(-l:l,-l:l),(2*l+1))
                DEALLOCATE(mat,conj_mat)
C D = P(icrorb,ik,is)*transpose(conjugate(P(icrorb,ik,is))) is a matrix of size (2*l+1)*(2*l+1)
                crdensmat(iss,icrorb)%mat(-l:l,-l:l)=
     =           crdensmat(iss,icrorb)%mat(-l:l,-l:l)
     +           +D(-l:l,-l:l)*kp(ik,is)%weight
C The weight used is a geometric factor asoociated to k and does not depend on the variable is. 
C That's why we merely multiply by the "weight" each term while summing over the k points.
               ENDIF
             ENDDO    ! End of the iss loop  
           ENDDO      ! End of the ik loop
           DEALLOCATE(D)
          ENDIF       ! End of the basis if-then-else
C
        ENDDO         ! End of the icrorb loop
C
C
C ===============================
C Symmetrization to the full BZ : 
C ===============================
         CALL symmetrize_mat(crdensmat,crorb,ncrorb)
C
C
C ============================================================================
C Application of the Rloc transformation to go back to the local coordinates :
C ============================================================================
         CALL rotdens_mat(crdensmat,crorb,ncrorb)
C
C
C =========================================================
C Finding the basis which diagonalizes the density matrix :
C =========================================================
C
        DO icrorb=1,ncrorb
          isrt=crorb(icrorb)%sort
          l=crorb(icrorb)%l
C
C -----------------------------------------------------------------------------------
C The s-orbitals are a particular case of a "non-mixing" basis and are treated here :
C -----------------------------------------------------------------------------------
          IF (l==0) THEN
C In this case, nothing is done.
           WRITE(buf,'(a)')'The correlated orbitals are s-orbitals.'
           CALL printout(0)
           WRITE(buf,'(a)')'END OF THE PRGM'
           CALL printout(0)
C Modified 20/09/10 - deleting STOP
C           STOP
C
C -----------------------------------------------------------------------------------------------------
C If the basis representation needs a complete spinor rotation approach (matrices of size 2*(2*l+1) ) :
C -----------------------------------------------------------------------------------------------------
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
           ALLOCATE(cmat(1:2*(2*l+1),1:2*(2*l+1)))
           cmat(1:2*(2*l+1),1:2*(2*l+1))=
     &       crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))
C Calculations of eigenvalues and eigenvectors
           ALLOCATE(work(4*(2*l+1)-1))
           ALLOCATE(rwork(6*(2*l+1)-2))
           ALLOCATE(eigen(1:2*(2*l+1)))
           CALL ZHEEV('V', 'U', 2*(2*l+1), cmat, 2*(2*l+1), 
     &       eigen, work,4*(2*l+1)-1,rwork,info)
           DEALLOCATE(work,rwork)
C If there is a pbm in zheev during execution, the prgm stops.
           IF (info.ne.0) THEN
            WRITE(buf,'(a)')
     &       'The subroutine zheev ends with info = ',info
            CALL printout(0)
            WRITE(buf,'(a)')'In newbasis, a pbm occurs in zheev.'
            CALL printout(0)
            WRITE(buf,'(a)')'END OF THE PRGM'
            CALL printout(0)
            STOP
           ENDIF       
C If info=0, eigen contains the eigenvalues 
C and cmat the desccription of the eigenvectors in the |new_i> basis: cmat=<new_i|eig_j>
           cmat(1:2*(2*l+1),1:2*(2*l+1))=
     =       MATMUL(TRANSPOSE(CONJG(reptrans(l,isrt)%transmat
     &       (1:2*(2*l+1),1:2*(2*l+1)))),
     &       cmat(1:2*(2*l+1),1:2*(2*l+1)))
C the eigenvectors are transformed in the complex basis.
C cmat = transpose(reptrans).cmat= <lm|new_i>.<new_i|eig_j>
C
C Writing the results in the file fort.25
           DO m=1,2*(2*l+1)
             WRITE(25,'(28(F10.7,x))') cmat(:,m)        
           ENDDO
           WRITE(25,'()')
           WRITE(25,'()')
           WRITE(25,'(a,i2,a,i2,a)') 'Eigenvectors for isrt= ',isrt,
     &       ' and l= ',l , ' calculated from'
           WRITE(25,'(a,f10.5,a,f10.5,a)') 
     &        'the density for the energy window from',
     &        e_bot,' Ry  to  ',e_top,' Ry'
           WRITE(25,'()')
           WRITE(25,'(a)') 'The corresponding eigenvalues are :'
           WRITE(25,'(14(F12.9,2x))') eigen(:)
           WRITE(25,'()')
           WRITE(25,'()')
C
           DEALLOCATE(eigen,cmat) 
C
C ----------------------------------------------------------------------------------------------
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only) :
C ----------------------------------------------------------------------------------------------
          ELSE
           IF(.not.ifSP) THEN
C The calculation is not spinpolarized (and then without SO), nsp=1
            ALLOCATE(cmat(1:2*l+1,1:2*l+1))
            ALLOCATE(zeros(1:2*l+1,1:2*l+1))
            zeros(1:2*l+1,1:2*l+1)=0.d0
            cmat(1:2*l+1,1:2*l+1)=
     &        crdensmat(1,icrorb)%mat(-l:l,-l:l)
C Calculations of eigenvalues and eigenvectors
            ALLOCATE(work(2*(2*l+1)-1))
            ALLOCATE(rwork(3*(2*l+1)-2))
            ALLOCATE(eigen(1:2*l+1))
            CALL ZHEEV('V', 'U', 2*l+1, cmat, 2*l+1, 
     &        eigen, work,2*(2*l+1)-1,rwork,info)
            DEALLOCATE(work,rwork)
C If there is a pbm in zheev during execution, the prgm stops.
            IF (info.ne.0) THEN
             WRITE(buf,'(a)')
     &        'The subroutine zheev ends with info = ',info
             CALL printout(0)
             WRITE(buf,'(a)')'In newbasis, a pbm occurs in zheev.'
             CALL printout(0)
             WRITE(buf,'(a)')'END OF THE PRGM'
             CALL printout(0)
             STOP
            ENDIF       
C If info=0, eigen contains the eigenvalues 
C and cmat the desccription of the eigenvectors in the |new_i> basis: cmat=<new_i|eig_j>
            cmat(1:2*l+1,1:2*l+1)=
     =        MATMUL(TRANSPOSE(CONJG(reptrans(l,isrt)%transmat
     &        (-l:l,-l:l) )),cmat(1:2*l+1,1:2*l+1))
C the eigenvectors are transformed in the complex basis.
C cmat = transpose(reptrans).cmat= <lm|new_i>.<new_i|eig_j>
C
C Writing the results in the file fort.25
            DO m=1,2*l+1
              WRITE(25,'(28(F10.7,x))') cmat(:,m),zeros(:,m)        
            ENDDO
            DO m=1,2*l+1
              WRITE(25,'(28(F10.7,x))') zeros(:,m),cmat(:,m)        
            ENDDO
            WRITE(25,'()')
            WRITE(25,'()')
            WRITE(25,'(a,i2,a,i2,a)') 'Eigenvectors for isrt= ',isrt,
     &        ' and l= ',l , ' calculated from'
            WRITE(25,'(a,f10.5,a,f10.5,a)') 
     &        'the density for the energy window from',
     &        e_bot,' Ry  to  ',e_top,' Ry'
            WRITE(25,'()')
            WRITE(25,'(a)') 'The corresponding eigenvalues are :'
            WRITE(25,'(14(F12.9,2x))') eigen(:)
            WRITE(25,'()')
            WRITE(25,'()')
C
            DEALLOCATE(eigen,cmat,zeros) 
C
           ELSE
C The calculation is spinpolarized (SO may be taken into account)
            ALLOCATE(cmat(1:2*(2*l+1),1:2*(2*l+1)))
            cmat(1:2*(2*l+1),1:2*(2*l+1))=0.d0
            DO iss=1,nsp
C Determination of the block indices :
              IF(iss.LE.2) THEN
               is=iss
               is1=iss 
C If iss=1 (up), is=1 and is1=1 -> Description of the up/up block
C If iss=2 (down), is=2 and is1=2 -> Description of the dn/dn block
              ELSE 
               is=iss-2 
               is1=3-is
C If iss=3, is=1 (up) and is1=2 (down) -> Description of the up/dn block
C If iss=4, is=2 (down) and is1=1 (up) -> Description of the dn/up block
              ENDIF
              cmat(1+(is-1)*(2*l+1):is*(2*l+1),
     &          1+(is1-1)*(2*l+1):is1*(2*l+1))=
     &          crdensmat(iss,icrorb)%mat(-l:l,-l:l)
            ENDDO
            WRITE(26,*) icrorb
            DO m=1,2*(2*l+1)
              WRITE(26,'(28(F15.10,x))') cmat(:,m)        
            ENDDO
            WRITE(26,'()')
C Calculations of eigenvalues and eigenvectors
            ALLOCATE(work(4*(2*l+1)-1))
            ALLOCATE(rwork(6*(2*l+1)-2))
            ALLOCATE(eigen(1:2*(2*l+1)))
            CALL ZHEEV('V', 'U', 2*(2*l+1), cmat, 2*(2*l+1), 
     &        eigen, work,4*(2*l+1)-1,rwork,info)
            DEALLOCATE(work,rwork)
C If there is a pbm in zheev during execution, the prgm stops.
            IF (info.ne.0) THEN
             WRITE(buf,'(a)')
     &        'The subroutine zheev ends with info = ',info
             CALL printout(0)
             WRITE(buf,'(a)')'In newbasis, a pbm occurs in zheev.'
             CALL printout(0)
             WRITE(buf,'(a)')'END OF THE PRGM'
             CALL printout(0)
             STOP
            ENDIF
            DO m=1,2*(2*l+1)
              WRITE(26,'(28(F10.7,x))') cmat(:,m)        
            ENDDO
            WRITE(26,'()')      
C If info=0, eigen contains the eigenvalues 
C and cmat the desccription of the eigenvectors in the |new_i> basis: cmat=<new_i|eig_j>
            ALLOCATE(tmat(1:2*(2*l+1),1:2*(2*l+1)))
            tmat(1:2*(2*l+1),1:2*(2*l+1))=0.d0
            tmat(1:(2*l+1),1:(2*l+1))=
     &        reptrans(l,isrt)%transmat(-l:l,-l:l)
            tmat((2*l+2):2*(2*l+1),(2*l+2):2*(2*l+1))=
     &        reptrans(l,isrt)%transmat(-l:l,-l:l)
            cmat(1:2*(2*l+1),1:2*(2*l+1))=
     =        MATMUL(TRANSPOSE(CONJG(tmat
     &        (1:2*(2*l+1),1:2*(2*l+1)))),
     &        cmat(1:2*(2*l+1),1:2*(2*l+1)))
            DO m=1,2*(2*l+1)
              WRITE(26,'(28(F10.7,x))') tmat(:,m)        
            ENDDO
            WRITE(26,'()')
            DO m=1,2*(2*l+1)
              WRITE(26,'(28(F10.7,x))') cmat(:,m)        
            ENDDO
            WRITE(26,'()')
            WRITE(26,'()')
C the eigenvectors are transformed in the complex basis.
C cmat = transpose(reptrans).cmat= <lm|new_i>.<new_i|eig_j>
C
C Writing the results in the file fort.25
            DO m=1,2*(2*l+1)
              WRITE(25,'(28(F10.7,x))') cmat(:,m)        
            ENDDO
            WRITE(25,'()')
            WRITE(25,'()')
            WRITE(25,'(a,i2,a,i2,a)') 'Eigenvectors for isrt= ',isrt,
     &        ' and l= ',l , ' calculated from'
            WRITE(25,'(a,f10.5,a,f10.5,a)') 
     &        'the density for the energy window from',
     &        e_bot,' Ry  to  ',e_top,' Ry'
            WRITE(25,'()')
            WRITE(25,'(a)') 'The corresponding eigenvalues are :'
            WRITE(25,'(14(F12.9,2x))') eigen(:)
            WRITE(25,'()')
            WRITE(25,'()')
C
            DEALLOCATE(eigen,cmat,tmat) 
C
           ENDIF      ! End of the ifSP if-then-else
          ENDIF       ! End of the basis if-then-else
        ENDDO         ! End of the icrorb loop
C
        ENDIF  ! End of the if ncrorb=0 if-then-else
      RETURN
      END

