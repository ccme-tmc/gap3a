      real(8) function coul_coef(qvec, iop_cut)
      ! Return the coefficient used when dealing with head and wing.
      ! Essentially the options of iop_cut are same to iop_coul(_c/x)
      ! in the barcoul module. 
      ! Note that q_para_vec should be the component of q parallel to 
      ! the BZ plane(2D) or direction(1D)
      ! Returns
      !     4\pi   
      !            for iop_cut = -1,0,3,4
      !     4\pi(1-\exp(-q_{\parallel}*zcut_coul*\cos(zcut_coul*q_{\perp}))
      !            for iop_cut = 2
      ! TODO: support for 1D, i.e. iop_cut = 1

          use constants, only: pi
          use barcoul,   only: zcut_coul, acut_coul, bcut_coul
          implicit none 
          real(8),intent(in) :: qvec(3)
          integer,intent(in) :: iop_cut
          real(8) :: t
      
          real(8),external :: veclen
      
          select case (iop_cut)
            case (2)
              t = 1.0D0 - exp(-veclen(qvec)*zcut_coul)
            case default
              t = 1.0D0
          end select 
      
          coul_coef = 4.0D0 * pi * t
      
      end function coul_coef
  
