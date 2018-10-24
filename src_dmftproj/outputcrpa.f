        SUBROUTINE outcrpa(elecn,qbbot,emin,emax)
C      
C Output data (projectors) for cRPA
C
        USE almblm_data
        USE common_data
        USE file_names
        USE prnt
        USE reps
        USE symm
        USE projections
        IMPLICIT NONE
C
C
        INTEGER :: iorb, icrorb, irep, isrt
        INTEGER :: icrorb2,iat,iatcorr,iatcorr2,icr
        INTEGER :: l,m,is,i1,i2
        INTEGER :: ik,il,ib,ir,n
        INTEGER :: ibmin,ibmax,tmp1,tmp2
        INTEGER, ALLOCATABLE :: ibmins(:),ibmaxs(:)  
        INTEGER :: ind1,ind2,iatom
        INTEGER :: ii,jj,nsymetr,is2,ispin,m2,iiden,iss
        INTEGER :: ikbz,ikbz0,ikbz1,ikbz2,ikbz3
        INTEGER :: nksym,nksym2
        REAL(KIND=8) ::  qbbot, elecn,emin,emax
        REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: kpoint
        REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: kpBZ
        REAL(KIND=8) :: dum1,dum2,dum3,det
        REAL(KIND=8) :: test1,test2,test3,diff1,diff2,diff3,tol
        REAL(KIND=8),DIMENSION(1:3) :: kvec,kpBZvec        
        logical symineq

!***********************************************************************************
c****************** Read in case.klist the cartesian coordinates of the K in BZ ****
!***********************************************************************************
c 5 corresponds to the number of columns defined in case.klist
c kpoint(ik,ik:kx:ky:kz:M)

        ALLOCATE(kpoint(1:nk,1:5))
c        write(oucRPA,'(a)')'K-points in IBZ'
        DO ik=1,nk
           read(iuklist2,*)(kpoint(ik,il),il=1,5)
        ENDDO

!        DO ikbz=1,nk*nsym
!           write(oucRPA,*)(kpBZ(ikbz,il),il=1,5)
!        ENDDO

!**************************************************************************
!**************************************************************************
!********************** WRITE THE OUTPUT FILE *****************************
!**************************************************************************
!**************************************************************************

        WRITE(buf,'(a)')'Writing the file case.projcrpa...'
        CALL printout(0)

        WRITE(oucRPA,'(f8.4,1x,f8.4,5x,a)') emin/rev,emax/rev
     &     ,'% unit of energy eV (Wannier window)'
        WRITE(oucRPA,'(i4,18x,a)')nk,'% number of k-points'

! Number of correlated orbitals counting only orbital moment l 
        WRITE(oucRPA,'(a,x,a)')
     &'% --------- num. of correlated orbitals,total num. of atom sorts'
        WRITE(oucRPA,'(20(i4,2x))') ncrorb,nsort
        WRITE(oucRPA,'(a)')'% --------- multiplicity of each atom sort'
        WRITE(oucRPA,'(20(i4,2x))') nmult(1:nsort)

!***********************************************************************
C********* Description of each correlated orbital "icrorb" *************
!***********************************************************************
       WRITE(oucRPA,'(a)')'% --------- atom| sort| lmax| l| m| rep| sym'
       DO icrorb=1,ncrorb 
           l=crorb(icrorb)%l
           isrt=crorb(icrorb)%sort
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
           IF (l==0) THEN
C For the s-orbitals, the only irep possible is the matrix itself.
              IF (ifSO) THEN
C If SO is taken into account, spinor rotation matrix is considered. 
C The spin is not a good quantum number, so the whole representation is used. 
                 WRITE(oucRPA,'(6(i6,x),a)') crorb(icrorb)%atom,isrt,0,
     &                2,crorb(icrorb)%ifSOat,1,defbasis(isrt)%typebasis
              ELSE
C Without SO, only the rotation matrix in orbital space is necessary.
                 WRITE(oucRPA,'(5(i6,x),a)') crorb(icrorb)%atom,isrt,
     &        orbitcorr(isrt)%orbmax,0,
C     &        1,crorb(icrorb)%ifSOat,1
C     &        1
     &        1,defbasis(isrt)%typebasis
              ENDIF
C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
           ELSEIF(reptrans(l,isrt)%ifmixing) THEN
C In this case, the SO is necessary considered, spinor rotation matrices are used.
              IF (crorb(icrorb)%ifsplit) THEN
C If only an irep is correlated
                 DO irep=1,reptrans(l,isrt)%nreps
                    IF (crorb(icrorb)%correp(irep)) THEN
                       WRITE(oucRPA,'(6(i6,x))') crorb(icrorb)%atom, 
     &                      isrt,l,reptrans(l,isrt)%dreps(irep),
     &                      crorb(icrorb)%ifSOat,irep
                    ENDIF
                 ENDDO
              ELSE
C If no particular irep is correlated
                 WRITE(oucRPA,'(7(i6,x))') crorb(icrorb)%atom,isrt,
     &                orbitcorr(isrt)%orbmax,l,
     &                2*(2*l+1),crorb(icrorb)%ifSOat,1
              ENDIF
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
           ELSE
              IF (ifSO) THEN 
C If SO is taken into account, spinor rotation matrices are considered.
C The spin is not a good quantum number, so the whole representation is used. 
                 WRITE(oucRPA,'(6(i6,x))') crorb(icrorb)%atom,isrt,l,
     &        2*(2*l+1),crorb(icrorb)%ifSOat,1
              ELSE
C Without SO, only the rotation matrix in orbital space is necessary.
                 IF (crorb(icrorb)%ifsplit) THEN
C If only an irep is correlated
                    DO irep=1,reptrans(l,isrt)%nreps
                       IF (crorb(icrorb)%correp(irep)) THEN
                         WRITE(oucRPA,'(5(i6,x),a)')crorb(icrorb)%atom,
     &                                   isrt,orbitcorr(isrt)%orbmax,l,
     &                                    reptrans(l,isrt)%dreps(irep),
     &                                    defbasis(isrt)%typebasis
                       ENDIF
                    ENDDO
                 ELSE
C If no particular irep is correlated

                   WRITE(oucRPA,'(5(i6,x),a)')crorb(icrorb)%atom,isrt,
     &                   orbitcorr(isrt)%orbmax,l,
     &                   (2*l+1),defbasis(isrt)%typebasis
                 END IF         ! End of the ifsplit if-then-else
              END IF            ! End of the ifSO if-then-else
           END IF               ! End of the ifmixing if-then-else
        END DO                  ! End of the icrorb loop
C an orbital "iorb" is described by :
C  - the associated atom
C  - the corresponding atomic sort
C  - the max orbital number corresponded to an atom 
C  - the considered orbital number l
C  - the size of the "correlated" submatrix (can be the whole matrix) 
C  - the flag ifSOat which states that SO is considered or not for this orbital
C  - the number of the irep


!**************************************************************************
!****** Write the k-points in full BZ (in cartesian coordinates)***********
!****** In the units of (2pi/a,2pi/b,2pi/c) *******************************
!**************************************************************************
!        DO ispin=1,ns
!           DO ikbz=1,nk
!              write(oucRPA,906)(NINT(kpoint(ikbz,il)),il=1,5)
!           ENDDO
!        ENDDO

!**************************************************************************
!********************** IBMIN AND IBMAX for all K *************************
!**************************************************************************
        ALLOCATE(ibmins(ns),ibmaxs(ns))
        ibmins = 0
        ibmaxs = 0

        DO ispin=1,ns
c The eigenvalues have no spin quatum number for Spin-orbit case
           IF (ifSO.and.is.eq.2) cycle
           DO ikbz=1,nk 
              ibmins(ispin)=kp(ikbz,ispin)%nb_bot
              ibmaxs(ispin)=kp(ikbz,ispin)%nb_top
              DO ikbz2=1,nk
               if(ikbz2.lt.ikbz) then
                 ibmins(ispin)=MIN(ibmins(ispin),kp(ikbz2,ispin)%nb_bot)
                 ibmaxs(ispin)=MAX(ibmaxs(ispin),kp(ikbz2,ispin)%nb_top)
               endif
              enddo !!ikbz2        
           ENDDO !!ikbz
        ENDDO !! ispin

        tmp1 = ibmins(1)
        tmp2 = ibmaxs(1)
        DO ispin = 1,ns
           if(tmp1.ge.ibmins(ispin)) then
             tmp1 = ibmins(ispin)
           endif
           if(tmp2.le.ibmaxs(ispin)) then
             tmp2 = ibmaxs(ispin)
           endif
        ENDDO !!ispin 
        ibmin = tmp1 
        ibmax = tmp2
        WRITE(oucRPA,'(a)')
     &     '% --------- min.KS band | max.KS band (over all k)'
        WRITE(oucRPA,*)ibmin,ibmax 
 
!**************************************************************************
!**************************** PROJECTORS **********************************
!**************************************************************************
        WRITE(oucRPA,'(a)')
     &    '% --------- h |k |l |M |weight |min.KS |max.KS ' 
        WRITE(oucRPA,'(a)')'% --------- Wannier projectors' 

        DO ispin=1,ns
          DO ikbz=1,nk
              write(oucRPA,904)(NINT(kpoint(ikbz,1+il)),il=1,4)
     &             ,kp(ikbz,ispin)%weight,kp(ikbz,ispin)%nb_bot,
     &             kp(ikbz,ispin)%nb_top

              DO ib=kp(ikbz,ispin)%nb_bot,kp(ikbz,ispin)%nb_top
                 DO icrorb=1,ncrorb
                    l=crorb(icrorb)%l
                    isrt=crorb(icrorb)%sort
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
                    IF (l==0) THEN
C For the s-orbitals, the only irep possible is the matrix itself.

                       WRITE(oucRPA,902)ib,isrt,icrorb,0, 
     &                      REAL(pr_crorb(icrorb,ikbz,ispin)%mat_rep(1,
     &                      ib)),
     &                      AIMAG(pr_crorb(icrorb,ikbz,ispin)%mat_rep(1,
     &                      ib))

C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
                    ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the SO is necessary considered, spinor rotation matrices are used.
                       IF(crorb(icrorb)%ifsplit) THEN
C If only 1 irep is correlated
                          ind1=1
                          DO irep=1,reptrans(l,isrt)%nreps
                             IF(crorb(icrorb)%correp(irep)) THEN
                                ind2=ind1+reptrans(l,isrt)%dreps(irep)-1
                                DO m=ind1,ind2
                                   WRITE(oucRPA,902) ib,isrt,icrorb,m,
     &                                  REAL(pr_crorb(icrorb,ikbz,ispin)
     &                                  %mat_rep(m,ib))
     &                                ,AIMAG(pr_crorb(icrorb,ikbz,ispin)
     &                                  %mat_rep(m,ib))
     
                                ENDDO

                             ENDIF
                             ind1=ind1+reptrans(l,isrt)%dreps(irep)
                          ENDDO

                       ELSE
C If no particular irep is correlated
                          DO m=1,2*(2*l+1)
                             WRITE(oucRPA,902)ib,isrt,icrorb,m, 
     &                          REAL(pr_crorb(icrorb,ikbz,1)%mat_rep(m,
     &                            ib)) 
     &                         ,AIMAG(pr_crorb(icrorb,ikbz,1)%mat_rep(m,
     &                            ib))
                          ENDDO
                          
                       ENDIF
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
                    ELSE
                       IF ((.not.ifSO).AND.crorb(icrorb)%ifsplit) THEN
C If only 1 irep is correlated (case without SO)
                          ind1=-l
                          DO irep=1,reptrans(l,isrt)%nreps
                             IF(crorb(icrorb)%correp(irep)) THEN
                                ind2=ind1+reptrans(l,isrt)%dreps(irep)-1
                                DO m=ind1,ind2
                                   WRITE(oucRPA,902)ib,isrt,icrorb,m, 
     &                                  REAL(pr_crorb(icrorb,ikbz,ispin)
     &                                  %mat_rep(m,ib)) 
     &                                ,AIMAG(pr_crorb(icrorb,ikbz,ispin)
     &                                  %mat_rep(m,ib))
                                ENDDO
                             ENDIF
                             ind1=ind1+reptrans(l,isrt)%dreps(irep)
                          ENDDO
                       ELSE
C If not only particular irep is correlated (case with and without SO)
                          DO m=-l,l
                             WRITE(oucRPA,902)ib,isrt,icrorb,m, 
     &                            REAL(pr_crorb(icrorb,ikbz,ispin)
     &                            %mat_rep(m,ib))
     &                            ,AIMAG(pr_crorb(icrorb,ikbz,ispin)
     &                            %mat_rep(m,ib))
                          ENDDO

                       END IF   ! End of the ifsplit if-then-else
                    END IF      ! End of the ifmixing if-then-else
                 END DO         ! End of the icrorb loop
              END DO            ! End of the ib loop
          ENDDO                 ! End of the ikbz loop
        ENDDO                   ! End of the ispin loop
C for each k-point and each correlated orbital, the corresponding projector is described by :
C  - the real part of the "correlated" submatrix
C  - the imaginary part of the "correlated" submatrix 
c*************************************************************************************

        DEALLOCATE(kpoint,ibmins,ibmaxs)

 900    format (f15.8,4(2x,f15.8))
 901    format (3(2x,f8.3))
!902    format (4(i3,2x),2(f15.10,2x))
 902    format (i6,2x,3(i3,2x),2(f15.10,2x))
 903    format (2(i3,2x),2(3f15.10,4x))
 904    format (4(i6),f15.8,2(i4,2x))
 905    format (5(i4,2x))
 906    format (5(9x,i4))
        RETURN
        END


