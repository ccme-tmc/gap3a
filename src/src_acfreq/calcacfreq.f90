!BOP
!
! !ROUTINE: calcacfreq
!
! !INTERFACE: 
      subroutine calcacfreq(init,iac,nomeg,omega,sc,npar,apar,z,fz,dfz)

! !DESCRIPTION:
!
! This subroutine sets up the analytical continuation (AC) parameters, 
! based on a set of data along the imaginary frequency and then
! calculate the value of the function at z  
! !USES:

      implicit none
      integer,   intent(in):: init   ! (0/1) -- control whether to calculate AC parameters (0 -> yes, 1 -> no)
      integer,   intent(in):: iac    ! Option for which method to use 
                                  ! 0  -- Thiele's reciprocal difference method as described in
                                  !       H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977)
                                  ! 1  -- Rojas, Godby and Needs (PRL 74, 1827 (1996)

      integer,   intent(in)   :: nomeg        ! Number of frequency points along the imaginary axis
      integer,   intent(in)   :: npar         ! Number of parameters (== 2*npol (npol= the number of poles)  
      real(8),   intent(in)   :: omega(nomeg) ! frequency points along the imaginary frequency 
      complex(8),intent(in)   :: sc(nomeg)    ! correlation selfenergy along the imaginary frequency 
      complex(8),intent(inout):: apar(npar)   ! fitted paramters to calculate selfenergy 
      complex(8),intent(in)   :: z
      complex(8),intent(out)  :: fz
      complex(8),intent(out)  :: dfz
      
     
! !LOCAL VARIABLES:   
      integer :: iw,step,ip,ierr
      integer :: iwpa(npar)
      real(8) :: varsq    ! Square root of the mean square error
      real(8)    :: xin(nomeg),anl(2*npar)      ! Values of the function
      complex(8) :: yin(nomeg)             ! Values of the function
      complex(8) :: cx(npar),cy(npar)        ! input for Pade's approximation 
      complex(8) :: coefs(0:npar/2)
      complex(8) :: comega(npar)
      logical::lpolish


! !INTRINSIC ROUTINES: 
      
      intrinsic dble

! !EXTERNAL ROUTINES: 
      external setrgn
      external setpatrd

! !REVISION HISTORY:
! Created: Nov. 19, 2007 by  JH
!EOP
!BOC
       call setwpa(nomeg,npar,omega,iwpa)  !! Choose the frequency points for the Pade's approximation     

       if(init.eq.0) then 
         xin(1:nomeg)=omega(1:nomeg)
         yin(1:nomeg)=sc(1:nomeg)

         do ip=1,npar
           iw=iwpa(ip)
           cx(ip)=cmplx(0.d0,xin(iw))
           cy(ip)=yin(iw)
         enddo

         if(iac.eq.0 .or. npar .eq. nomeg)  then 
           call setpatrd(npar,cx,cy,apar)
         else 
           call init_c(cx,cy,apar,npar)

           do ip=1,npar
             anl(ip)=real(apar(ip))
             anl(ip+npar)=aimag(apar(ip))
           enddo
           call nllsq(xin,yin,nomeg,anl,2*npar,varsq)
           do ip=1,npar
             apar(ip)=cmplx(anl(ip),anl(ip+npar),8)
           enddo
         endif 
       endif 

       if(iac.eq.1) then
         call acrgn(npar,z,apar,fz,dfz)
       elseif(iac.eq.0) then
         do ip=1,npar
           iw=iwpa(ip)
           comega(ip)=cmplx(0.d0,omega(iw))
         enddo
         call acpatrd(npar,z,comega,apar,fz,dfz)
       endif 
       return
          
      endsubroutine 
!EOC        
