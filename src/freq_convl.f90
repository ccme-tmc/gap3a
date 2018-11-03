        subroutine freq_convl(iom,nomg,om,enk,mwm,omg,womg,sc)
!!
!! subroutine for the following frequency convolution 
!!            S(i\omega; \epsilon ) 
!!             = 1/\pi \int_0^\infity \frac{ \epsilon - i \omega ) W(i\omega')} 
!!               { (\epsilon - i \omega)^2 + \omega'^2} d\omega'
!!    
        implicit none
        integer,intent(in) :: iom   !! control how to do convolution 
                                    !  iom.eq.0 -- direct integraion 
                                    !  iom.gt.0 -- indicating that om =
                                    !  omg(iom) so that special treatment as explained in
                                    !  Developer Guide (H.4) is needed 
        integer,   intent(in) :: nomg  !! number of frequency points used for the convolution
        real(8),   intent(in) :: om    !! the frequency \omega 
        real(8),   intent(in) :: enk   !! \epsilon 
        complex(8),intent(in) :: mwm(nomg)  !! W(i\omega')
        real(8),   intent(in) :: omg(nomg)  !! the array for \omega'
        real(8),   intent(in) :: womg(nomg) !! integraion weight 
        complex(8),intent(out):: sc

        integer :: jom
        complex(8) :: ffac ! the frequecy dependent factor in each term
        complex(8) :: comsq
        complex(8) :: ommek   ! i*omega-e_nk
        real(8),parameter:: pi=3.14159265358979d0
        ommek=cmplx(enk,-1.0d0*om,8)
        sc=0.d0
        if(iom.gt.0) then
          do jom = 1, nomg    ! integration frequencies
            if(jom.eq.iom) cycle
            comsq=cmplx(omg(jom)*omg(jom),0.0d0,8)
            ffac =cmplx(womg(jom),0.0d0,8)/(ommek*ommek+comsq)
            sc=sc+ffac*(mwm(jom)-mwm(iom))
          enddo
          sc=ommek*sc/pi+mwm(iom)*sign(0.5d0,enk)
        else
          do jom = 1, nomg
            comsq=cmplx(omg(jom)*omg(jom),0.0d0,8)
            ffac =cmplx(womg(jom),0.0d0,8)/(ommek*ommek+comsq)
            sc=sc+ffac*mwm(jom)
          enddo
          sc=ommek*sc/pi
        endif
        end subroutine 
