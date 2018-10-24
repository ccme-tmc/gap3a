      subroutine calc_spectfun_nk(spf,oms,nom,exnk,sigc,oms_im,nom_im,fout)
!     This subroutine calculate the spectral function for a given state
!     (nk) based on the diagonal correlation selfenergy along the imaginary
!     frequency
      use constants,  only: hev
      use freq,       only: freq_eta 
      use selfenergy, only: iop_ac,npar_ac
      implicit none 
      integer,intent(in):: nom             ! the number of real freq.  
      integer,intent(in):: nom_im          ! the number of imag. freq. 
      real(8),intent(in):: oms(nom)        ! real freq. grid 
      real(8),intent(in):: oms_im(nom_im)  ! imag. freq. grid 
      real(8),intent(in):: exnk            ! the exchange part of qp energy 
      complex(8),intent(in):: sigc(nom_im) ! correlation self. energy
      real(8),intent(out):: spf(nom)       ! the spectral function 
      integer,intent(in)::fout 

      integer:: iom
      real(8):: omg
      complex(8),allocatable:: acpar(:,:) ! AC parameters for positive and negative frequencies
      real(8):: sc_re, sc_im 

      complex(8):: ein, sc, dsc
      real(8),parameter:: coef=0.318309886184 ! 1/Pi

      ! set AC parameters for positive frequencies 
      allocate(acpar(npar_ac,2)) 
      call setsac(iop_ac,nom_im,npar_ac, 1.0,sigc,oms_im,acpar(:,1))

      ! set AC parameters for negative frequencies 
      call setsac(iop_ac,nom_im,npar_ac,-1.0,sigc,oms_im,acpar(:,2))

      spf=0.d0 
      do iom=1,nom 
        omg = oms(iom) 
        ein = cmplx(omg,0.d0)
        if(omg.gt.0.0d0) then 
          call getsac(iop_ac,nom_im,npar_ac,1.0,ein,oms_im,acpar(:,1),sc,dsc) 
        else 
          call getsac(iop_ac,nom_im,npar_ac,-1.0,ein,oms_im,acpar(:,2),sc,dsc) 
        endif
        sc_re = real(sc) 
        sc_im = max(abs(dimag(sc)),freq_eta) 
        
        spf(iom) = coef*sc_im/((omg-exnk-sc_re)**2 + sc_im**2)
        if(fout.gt.0) then 
          write(fout,'(6F12.6)') omg*hev,spf(iom),sc_re,dimag(sc)  
        endif 
      enddo 
      deallocate(acpar) 
      end subroutine 

