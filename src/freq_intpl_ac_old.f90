      subroutine freq_intpl_ac(iac,f1,w1,nw1,f2,w2,nw2,iop) 
!
!     This subroutine interpolate a complex function to a discete set of frequency 
!     (real or imaginary, indicated by iop ), represented by "w2", by
!     using the data on a discrete set of positive imaginary frequency,
!     represented by w1 and f1 
! 
      implicit none 
      integer, intent(in):: iac     !! the option to choose the AC scheme 
      integer, intent(in):: nw1     !! the number of imag. freq. 
      integer, intent(in):: nw2     !! the number of freq. to be interpolated 
      real(8), intent(in):: w1(nw1) !!  
      real(8), intent(in):: w2(nw2) 
      complex(8), intent(in):: f1(nw1) 
      complex(8), intent(out)::f2(nw2) 
      integer, intent(in):: iop   ! 0/1/2 - positive imag./ positive real/ negative real   

      integer:: iw 
      integer:: npol  ! the number of poles 
      integer:: npar  ! Number of parameters (== 2*npol)
      real(8) :: w1t(nw1)
      complex(8):: z,fz,dfz 
      complex(8):: f1t(nw1),w2c(nw2)  
      complex(8), allocatable :: apar(:)   ! fitted paramters 
      
      if(iac.le.0) then 
        npol = nw1/2
      else 
        npol = iac + 1
      endif 
      npar=2*npol
      allocate(apar(npar))

      if(iop.eq.0) then 
        z = cmplx(0.d0,w2(1)) 
        w2c = cmplx(0.d0,w2) 
        w1t = w1 
        f1t = f1
        call calcacfreq(0,iac,nw1,w1t,f1t,npar,apar,z,f2(1),dfz)
      else if (iop.eq.1) then 
        z =   cmplx(w2(1),0.d0) 
        w2c = cmplx(w2,0.d0) 
        w1t = w1 
        f1t = f1
        call calcacfreq(0,iac,nw1,w1t,f1t,npar,apar,z,f2(1),dfz)
      else 
        z = cmplx(w2(1),0.d0)  
        w2c = cmplx(w2,0.d0) 
        w1t = - w1
        f1t = conjg(f1) 
        call calcacfreq(0,iac,nw1,w1t,f1t,npar,apar,z,f2(1),dfz)
      endif 

      do iw=2,nw2 
        call calcacfreq(1,iac,nw1,w1t,f1t,npar,apar,w2c(iw),f2(iw),dfz) 
      enddo 
      deallocate(apar) 

      end subroutine 
          
