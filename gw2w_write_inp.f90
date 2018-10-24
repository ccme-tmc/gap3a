!BOP
!
! !ROUTINE: io_sxcmn
!
! !INTERFACE: 
      subroutine gw2w_write_inp(isym,ierr)

! !DESCRIPTION:
!
!  This subroutine writes the full selfenergy matrix 
!  represented by Wannier functions  
!
! !USES:
      use bands,      only: nspin,eferks
      use crpa,       only: nlmorb,nbmin_wf, nbmax_wf 
      use constants,  only: hev
      use freq,       only: nomeg,omega
      use kpoints,    only: nkp,nirkp,get_kvec,nkdivs
      use selfenergy, only: sigm,sxc_wann,vxc_wann,hks_wann,beta_mats, &
     &                      omega_mats,nomeg_mats,sxc_wann_mats
      use struk,      only: alat,rbas,gbas
      use task,       only: casename,savdir
      
      implicit none
     
      integer,intent(in)   :: isym ! 0/1 - full/irreducible BZ for k-mesh 
      integer,intent(out)  :: ierr
      
      
! !LOCAL VARIABLES:
      integer :: fid 
      integer :: i,j,isp,iik,ik,irk,nktot
      integer :: iom
  
      real(8) :: kvec(3)
      character(80)::fname
      character(20):: sname="gw2w_write_inp "

! !REVISION HISTORY:
!
!EOP
!BOC
!

!
!     Set file name 
!

      if(isym.eq.0) then 
        nktot = nkp
      else
        nktot = nirkp
      endif 
      fid=999

      !! 
      !! write out the important parameters 
      !!
      fname=trim(casename)//".gw2wann_inp"
      open(unit=fid,file=fname,action='write',iostat=ierr)
      if(ierr.ne.0) then 
        write(6,*) "WARNING: fail to open "//trim(fname)
        return 
      endif  
      write(fid,11) nlmorb
      write(fid,12) nspin 
      write(fid,13) nomeg,nomeg_mats 
      write(fid,14) nkp,nkdivs(1:3)
      write(fid,15) beta_mats
      write(fid,16) eferks
      write(fid,10)"# Unit cell vectors in cartesian coordinates:"
      do i=1,3
        write(fid,17) rbas(1:3,i),i
      enddo 
      write(fid,10) "# Reciprocal lattice vectors "
      do i=1,3
        write(fid,18) gbas(1:3,i),i
      enddo 
  
      write(fid,10) "# k-vectors" 
      do ik=1,nkp
        call get_kvec(0,ik,iik,irk,kvec)
        write(fid,20) ik,kvec(1:3)  
      enddo 
 10   format(A) 
 11   format(I6,4X,'# Number of Wannier orbitals') 
 12   format(I6,4X,'# Number of spin channels')
 13   format(2I6,4X,'# Number of imag. freq.')
 14   format(4I6,4X,'# nr. of k-points, nkx, nky, nkz')   
 15   format(E16.8,4X,'# Beta (inv.Temp.) in 1/Ha')
 16   format(E16.8,4X,'# Chemical potential (in Ha)') 
 17   format(3F12.6,4X,'# a',i1)
 18   format(3F12.6,4X,'# b',i1) 
 20   format(i6,3F12.6,4X,'# ik, k(1:3)') 
      close(fid) 

      end subroutine gw2w_write_inp 
!EOC          
            
