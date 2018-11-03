!BOP
! !ROUTINE: bz_setkmesh
!
! !INTERFACE:
       subroutine bz_setkmesh

! !DESCRIPTION:
!
! This subroutine generates the k- and q-point grids   
! and initilizes all arrays used in the tetrahedral integration. 
! The crystal structure is read from the \verb"case.struct" file of the
! wien2k code. The result for the k-mesh has to be identical to that of the
! WIEN2k code. The q-mesh is the same but without k0shift. The tetrahedra 
! are different since the q-vector breaks their
! symmetry all of them have to be taken into account separatedly. The
! geometrical weights of the set of q points depends on the corresponding
! k-point...
!
! !USES:a
      use struk,   only: alat,ortho,lattic,imat,izmat,nsym,rbas,gbas,br2
      use kpoints, only: nkdivs,idvk,idvkir,idvq,kirlist, klist,    &
     &                   kpirind,kqid, nirkp, nkp,qlist,k0shift,&
     &                   weightk,weightq, wkir, nqp,   &
     &                   kirvecs,nkp_ene,kvecs
      use eigenvec, only: lsymvector
      use task,     only: fid_outkpt

! !LOCAL VARIABLES:

      implicit none
      character(len=20):: sname="bz_setkmesh"
      


! !EXTERNAL ROUTINES: 


      external cartezian
      external kgen
      external kqgen
      external sym2int
      external writeklist
      external writeqgen
      external writesym

!EOP
!BOC

!
!     Calculate the maximum number of k-points
!
      call linmsg(6,'-',sname) 
      nkp= nkdivs(1)*nkdivs(2)*nkdivs(3)
      nqp=nkp
      nirkp=nkp

!     Allocate the arrays for the kpoints
!
      allocate(kpirind(nkp))
      allocate(kqid(nkp,nqp))
      allocate(weightk(nkp))
      allocate(weightq(nqp))
      allocate(klist(3,nkp))
      allocate(qlist(3,nqp))

      call kgen(gbas,nsym,izmat,nkdivs,k0shift,2,nirkp,                 &
     &       klist,idvkir,kpirind,weightk)

      allocate(kirlist(3,nirkp),wkir(nirkp))

      kirlist(1:3,1:nirkp)=klist(1:3,1:nirkp)
      wkir(1:nirkp)=weightk(1:nirkp)

!
!     Generate the k- and q-points meshes
!
      write(6,*) " Number of symmetry operations: nsym=",nsym
      write(6,*) " Total number of k-points: nkp=",nkp
      write(6,*) " Number of irred. k-points:  nirkp=",nirkp

      call kqgen(gbas,nkdivs,k0shift,2,nkp,klist,qlist,idvk,idvq,kqid)

      call bz_setiksym

      weightq = 1    !! q is currently over full BZ

      write(6,*) " idvk=",idvk 
      write(6,*) " idvq=",idvq 

      write(6,*) ""
      write(fid_outkpt,*) "#k-mesh in internal coordinates"

      write(fid_outkpt,*) "#Irreducible k-points "
      call w2k_writeklist(fid_outkpt,nirkp,nkdivs,idvkir,wkir,kirlist)

      write(fid_outkpt,*) "# All k-points "
      call w2k_writeklist(fid_outkpt,nkp,nkdivs,idvk,weightk,klist)

      write(fid_outkpt,*) "# All q-points "
      call w2k_writeklist(fid_outkpt,nkp,nkdivs,idvq,weightq,qlist)
!
!     Transform the k- and q-points to cartesian coords. if necesary
!
      if(ortho.or.(lattic(1:3).eq.'CXZ'))then
        call cartezian(nirkp,alat,br2,kirlist,idvkir)
        call cartezian(nkp,alat,br2,klist,idvk)
        call cartezian(nqp,alat,br2,qlist,idvq)
      endif

      allocate(kirvecs(3,nirkp))
      kirvecs=dble(kirlist)/idvkir
      allocate(kvecs(3,nkp)) 
      kvecs=dble(klist)/idvk
!
!     Write k-points info in WIEN2k format
!
      write(fid_outkpt,*) "#k-mesh in WIEN2k format"
      write(fid_outkpt,*) "#Irreducible k-points "
      call w2k_writeklist(fid_outkpt,nirkp,nkdivs,idvkir,wkir,kirlist)

      write(fid_outkpt,*) "# All k-points "
      call w2k_writeklist(fid_outkpt,nkp,nkdivs,idvk,weightk,klist)

      write(fid_outkpt,*) "# All q-points "
      call w2k_writeklist(fid_outkpt,nkp,nkdivs,idvq,weightq,qlist)

!      write(fid_outkpt,*) "# the tetrahedra data "
!      call writeqgen(fid_outkpt)

      ! Check the consistency with the k-mesh used in WIEN2k 
      if(lsymvector) then 
        call w2k_readenergy(-1) 
        if(nirkp.ne.nkp_ene) then 
          write(6,'(a)') "ERROR: the irreducible k-mesh generated here &
     &is inconsistent with that used in WIEN2k!!!"
          write(6,'(a)') " - there may be some bugs in libbzint "
          write(6,'(a)') " - use KS-vectors in the full BZ (check the & 
     &flag '-sv' in gap2_init) to avoid such problems"
          stop 
        endif 
      endif 
      end subroutine bz_setkmesh
!EOm
