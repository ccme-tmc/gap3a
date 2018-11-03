!BOP
!
! !ROUTINE: calc_eqp_diag
!
! !INTERFACE:
      subroutine calc_eqp_diag(nsp,nk,nb,omg,nom,eks,vxc,sx,sc,eqp)
      use selfenergy, only: iop_ac,npar_ac,sacpar
      use constants,  only: hev 
      implicit none 

      integer,intent(in):: nsp ! the number of spins
      integer,intent(in):: nk  ! the number of k-points 
      integer,intent(in):: nb  ! the number of bands
      integer,intent(in):: nom ! the number of imaginary freq. 
      real(8),intent(in):: omg ! the imaginary freq. 
      real(8),intent(in):: eks(nb,nk,nsp) ! KS band energies 
      real(8),intent(in):: vxc(nb,nk,nsp) ! KS Vxc diagonal elements 
      real(8),intent(in):: sx(nb,nk,nsp) ! exchange selfenergy 
      complex(8),intent(in):: sc(nom,nb,nk,nsp) ! correlation selfenergies along imaginary freq.  
      real(8),intent(out):: eqp(nb,nk,nsp) ! QP energies 

! !DESCRIPTION:
! Given the matrix elements $\langle
!   \Psi_{n\vec{k}}|\Sigma(\vec{k},\omega)|\Psi_{n\vec{k}}\rangle$, $\langle
!   \Psi_{n\vec{k}}|V^{xc}|\Psi_{n\vec{k}}\rangle$ and
!   $\varepsilon^{DFT}_{n\vec{k}}$
! this subroutine calculates the quasi-particle energies $\varepsilon^{qp}_{n\vec{k}}$
! 

      integer:: ie,isp,ik   !(Counter) index for band,spin,k-points
      integer:: ierr
      real(8):: enk,vxcnk,znk,sxnk,scnk,snk,enk0
      complex(8):: sig,ein
      complex(8):: dsig

      real(8):: omg_ac(nom) ! Parameters of the selfenergya
      complex(8) ::sc_ac(nom)

      do isp=1,nspin      ! Loop over spin
        do ik=1,nk
          do ie=1,nb
            enk   = eks(ie,ik,isp)  
            vxcnk = vxc(ie,ik,isp)   
            sxnk  = sx(ie,ik,isp)

            ein=cmplx(enk,0.d0) 
            if(enk.ge.0.d0) then 
              omg_ac = omg 
              sc_ac = sc(1:nom,ie,ik,isp)
            else 
              omg_ac = -omg
              sc_ac = conjg(sc(1:nom,ie,ik,isp))
            endif 
        
            !! Calculate the new quasi-particle energy 
            call calcacfreq(0,iop_ac,nom,omg_ac,sc_ac,npar_ac, &
     &                sacpar(:,ie,ik,isp),ein,sig,dsig)
            znk=1.0d0/(1.0d0-real(dsig))
            scnk=real(sig)
            snk=scnk+sxnk

            if(znk.gt.1.d0 .or.znk.lt.0.d0) then
              write(6,*) "WARNING: nonphysical Znk found!"
              write(6,*) " -- check case.outqp for details "
              write(fout,100) ik,ie,enk*hev,dsig,znk
              znk=1.0
              dsig = 0.d0
            endif
            eqp(ie,ik,isp) = enk + znk*(snk-vxcnk)
          enddo ! ie
        enddo ! ik
      enddo ! isp

 100  format('WARNING: nonphysical Znk',/,'ik=',i4,'  ie=',i4, &
     &'  enk=',f8.3,' eV',/,' dsig=',2g12.4, '  Znk=',f8.3,/, &
     &' - reset to 1.0') 
      end subroutine calc_eqp_diag

!EOC
          
