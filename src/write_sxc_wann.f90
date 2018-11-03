!BOP
!
! !ROUTINE: io_sxcmn
!
! !INTERFACE: 
      subroutine write_sxc_wann(isym,iop_hgw,ierr)

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
      use selfenergy, only: sigm,sxc_wann,vxc_wann,hks_wann,hgw_wann,&
     &                 beta_mats,omega_mats,nomeg_mats,sxc_wann_mats
      use struk,      only: alat,rbas,gbas
      use task,       only: casename,savdir
      
      implicit none
     
      integer,intent(in)   :: isym ! 0/1 - full/irreducible BZ for k-mesh 
      integer,intent(in)   :: iop_hgw ! 0/1 - indicating whether output the effect GW Hamiltonian  
      integer,intent(out)  :: ierr
      
      
! !LOCAL VARIABLES:
      integer :: fid 
      integer :: i,j,isp,iik,ik,irk,nktot
      integer :: iom
  
      real(8) :: kvec(3)
      character(80)::fname
      character(20):: sname="write_sxc_wann"

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
 15   format(E16.8,4X,'# Beta (inv.Temp.)')
 16   format(E16.8,4X,'# Chemical potential') 
 17   format(3F12.6,4X,'# a',i1)
 18   format(3F12.6,4X,'# b',i1) 
 20   format(i6,3F12.6,4X,'# ik, k(1:3)') 
      close(fid) 
      
      !! 
      !! write out the LDA Kohn-Sham Hamiltonian in the Wannier basis 
      !!
      fname=trim(casename)//".hks_wann"
      open(unit=fid,file=fname,action='write',iostat=ierr)
      if(ierr.ne.0) then 
        write(6,*) "WARNING: fail to open "//trim(fname)
        return 
      endif  

      write(fid,100) "# Kohn-Sham Hamiltonian in the Wannier basis"
      do isp=1,nspin
        do ik=1,nkp 
          call get_kvec(0,ik,iik,irk,kvec)
          write(fid,110,advance='no') ik,kvec(1:3) 
          do j=1,nlmorb
            do i=1,nlmorb
              write(fid,120,advance='no') hks_wann(i,j,ik,isp)
            enddo 
          enddo  
          write(fid,100) ' '
        enddo
      enddo 
      close(fid) 

      !! 
      !! write out the effective GW Hamiltonian in the Wannier basis 
      !!
      if(iop_hgw.gt.0) then 
        fname=trim(casename)//".hgw_wann"
        open(unit=fid,file=fname,action='write',iostat=ierr)
        if(ierr.ne.0) then 
          write(6,*) "WARNING: fail to open "//trim(fname)
          return 
        endif 
        write(fid,100) "# Effective GW Hamiltonian in the Wannier basis"
        do isp=1,nspin
          do ik=1,nkp
            call get_kvec(0,ik,iik,irk,kvec)
            write(fid,110,advance='no') ik,kvec(1:3)
            do j=1,nlmorb
              do i=1,nlmorb
                write(fid,120,advance='no') hgw_wann(i,j,ik,isp)
              enddo
            enddo
            write(fid,100) ' '
          enddo
        enddo
        close(fid)
      endif 

      !! 
      !! write out the LDA Vxc in the Wannier basis 
      !!
      fname=trim(casename)//".vxc_wann"
      open(unit=fid,file=fname,action='write',iostat=ierr)
      if(ierr.ne.0) then 
        write(6,*) "WARNING: fail to open "//trim(fname)
        return 
      endif  
      write(fid,100) "# Vxc in the Wannier representation "
      do isp=1,nspin
        do ik=1,nkp
          call get_kvec(0,ik,iik,irk,kvec)
          write(fid,110,advance='no') ik, kvec(1:3)
          do j=1,nlmorb
            do i=1,nlmorb
              write(fid,120,advance='no') vxc_wann(i,j,ik,isp)
            enddo
          enddo
          write(fid,100) ' '
        enddo
      enddo
      close(fid)

      !! 
      !! write out the GW selfenergy (Sxc) in the Wannier basis 
      !!
      fname=trim(casename)//".sxc_wann"
      open(unit=fid,file=fname,action='write',iostat=ierr)
      if(ierr.ne.0) then
        write(6,*) "WARNING: fail to open "//trim(fname)
        return
      endif
      write(fid,100) "# Sxc in the Wannier representation "
      do isp=1,nspin
        do ik=1,nkp
          call get_kvec(0,ik,iik,irk,kvec)
          write(fid,130) ik,kvec(1:3),irk
          do iom=1,nomeg
            write(fid,140,advance='no') omega(iom)
            do j=1,nlmorb
              do i=1,nlmorb
                write(fid,120,advance='no') sxc_wann(i,j,iom,ik,isp)  
              enddo
            enddo
            write(fid,100) ' '
          enddo
          write(fid,100) ' '
        enddo
      enddo
      close(fid)

      !! 
      !! write out the GW selfenergy (Sxc) in the KS basis 
      !!
      fname=trim(casename)//".sxc_ks"
      open(unit=fid,file=fname,action='write',iostat=ierr)
      if(ierr.ne.0) then
        write(6,*) "WARNING: fail to open "//trim(fname)
        return
      endif
      write(fid,100) "# Sxc in the Kohn-Sham representation "
      do isp=1,nspin
        do iik=1,nktot
          call get_kvec(isym,iik,ik,irk,kvec)
          write(fid,131) iik,kvec(1:3),irk,ik
          do iom=1,nomeg
            write(fid,140,advance='no') omega(iom)
            do j= nbmin_wf,nbmax_wf 
              do i= nbmin_wf, nbmax_wf 
                write(fid,120,advance='no') sigm(i,j,0,iik,isp) & 
     &                                     +sigm(i,j,iom,iik,isp)
              enddo
            enddo
            write(fid,100) ' '   !! to end the output line 
          enddo
          write(fid,100) ' '  !! to output a blank line 
        enddo
      enddo
      close(fid)

      if(nomeg_mats.eq.0) return  

      !! 
      !! write out the GW selfenergy (Sxc) in the Wannier basis along
      !  the Matsubara frequency
      !!
      fname=trim(casename)//".sxc_wann_mats"
      open(unit=fid,file=fname,action='write',iostat=ierr)
      if(ierr.ne.0) then
        write(6,*) "WARNING: fail to open "//trim(fname)
        return
      endif
      write(fid,100) "# Sxc_wann interpolated to Matsubara freq."
      do isp=1,nspin
        do ik=1,nkp
          call get_kvec(0,ik,iik,irk,kvec)
          write(fid,130) ik, kvec(1:3),irk
          do iom=1,nomeg_mats 
            write(fid,140,advance='no') omega_mats(iom)
            do j=1,nlmorb
              do i=1,nlmorb
                write(fid,120,advance='no') sxc_wann_mats(i,j,iom,ik,isp)
              enddo
            enddo
            write(fid,100) ' '
          enddo
          write(fid,100) ' '
        enddo
      enddo
      close(fid)
      return
 100  format(a) 
 110  format(I6,3F12.6)  
 120  format(2g24.16) 
 130  format('#',I6,3F12.6,I6,4X,'# ik, k(1:3),irk') 
 131  format('#',I6,3F12.6,2I6,4X,'# iik, k(1:3),irk,ik') 
 140  format(F16.6)

      end subroutine write_sxc_wann
!EOC          
            
