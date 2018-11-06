      real(8) function coul_coef(q_para_vec,iop_cut)
      ! Return the coefficient used when dealing with head and wing.
      ! Essentially the options of iop_cut are same to iop_coul(_c/x)
      ! in the barcoul module. 
      ! Note that q_para_vec should be the component of q parallel to 
      ! the 2D BZ
      ! Returns
      !     4\pi                                  iop_cut = -1,0,3,4
      !     4\pi(1-\exp(-abs{q_{||}}*zcut_coul)   iop_cut = 2
      ! TODO: support for 1D, i.e. iop_cut = 1

          use constants, only: pi
          use barcoul,   only: zcut_coul, acut_coul, bcut_coul
          implicit none 
          real(8),intent(in) :: q_para_vec(3)
          integer,intent(in) :: iop_cut
          real(8) :: t
      
          real(8),external :: veclen
      
          select case (iop_cut)
            case (2)
              t = 1.0D0 - exp(-veclen(q_para_vec)*zcut_coul)
            case default
              t = 1.0D0
          end select 
      
          coul_coef = 4.0D0 * pi * t
      
      end function coul_coef
  
