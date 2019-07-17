!BOP
!
! !ROUTINE: io_sxcmn
!
! !INTERFACE: 
      subroutine gw2w_write_h(isym,ierr)

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
      character(20):: sname="gw2w_write_wann "

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
      !! write out the LDA Kohn-Sham Hamiltonian in the Wannier basis 
      !!
      fname=trim(casename)//".hks_wann"
      open(unit=fid,file=fname,action='write',iostat=ierr)
      if(ierr.ne.0) then 
        write(6,*) "WARNING: fail to open "//trim(fname)
        return 
      endif  

      ! Write the DFT Hamitlonian
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
      !! write out the LDA Vxc in the Wannier basis 
      !!
      fname=trim(casename)//".vxc_wann"
      open(unit=fid,file=fname,action='write',iostat=ierr)
      if(ierr.ne.0) then 
        write(6,*) "WARNING: fail to open "//trim(fname)
        return 
      endif  
      write(fid,100) "# Vxc in the Wannier representation "
      ! Write the DFT Hamitlonian
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
      ! Write the DFT Hamitlonian
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
      ! Write the DFT Hamitlonian
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
      ! Write the DFT Hamitlonian
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

      end subroutine gw2w_write_wann 
!EOC          
            
