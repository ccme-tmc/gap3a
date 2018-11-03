!BOP
!
! !ROUTINE: writeqgen
!
! !INTERFACE: 
      subroutine writeqgen(fid)

! !DESCRIPTION:
!
! This subroutine writes the identification number of the nodal points of
! each tetrahedron to file case.qgen
!
! !USES:
      use kpoints, only: kqid, link, nkp, ntet, tnodes, tvol

! !LOCAL VARIABLES:

      implicit none
      integer(4),intent(in) :: fid

      integer(4) :: i,j,k,itet,kmax
!EOP
!BOC 
      write(fid,*) 
      write(fid,*) "#writeeqgen: Nodal points of tetrahedron"
  
      write(fid,101) nkp,ntet,tvol
      do itet=1,ntet
        write(fid,100)itet,(tnodes(i,itet),i=1,4),(link(itet,j),j=1,nkp)
      enddo 
!
!     Write the k-dependent q and k' weights
!
      write(fid,*)' k,k-q links'
      do i=1,nkp,8
        kmax=7
        if((nkp-i).lt.kmax)kmax=nkp-i
        write(fid,103)(i+k,k=0,kmax)
        do j=1,nkp
          write(fid,102)j,(kqid(j,i+k),k=0,kmax)
        enddo
      enddo

  100 format(i6,4i4,4i6)
  101 format(2i6,e16.8)
  102 format(i4,'   |',8i4)      
  103 format(' ik/iq |',8i4)      
      
      end subroutine writeqgen
!EOC
