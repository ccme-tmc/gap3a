!BOP
!
! !ROUTINE: edifwt
!
! !INTERFACE:

      subroutine stweight_imag(deltae_vert,omeg,equiv_flag,weight_vert)
!
! !DESCRIPTION:
!
! This subroutine calculates the weight on the whole small tetrahedron
! in which the $k$ states are fully occupied and $k-q$ states are fully 
! unoccupied. This is for the $sigfreq=3$ case when we consider the 
! imaginary frequency. 
!                                             
! !USES:
      use tetra_internal, only: weighttol, weightwarn,fout, &
     &                          vol_small_tetra

! !INPUT PARAMETERS:
      implicit none
      
      real(8), intent(in) :: deltae_vert(4)    ! difference of the energy 
                                               ! in k-mesh tetrahedron vertices 
                                               ! and k-q mesh tetrahedron vertices.

      real(8), intent(in) :: omeg              ! the frequency omega to be calculated

      integer(4), intent(in) :: equiv_flag     ! == 4, none is equal 
                                               ! == 6, deltae_vert(1)=deltae_vert(2).
                                               ! ==8,  deltae_vert(1)=deltae_vert(2) and deltae_vert(3)=deltae_vert(4).
                                               ! ==10, deltae_vert(1)=deltae_vert(2)=deltae_vert(3).
                                               ! ==16, deltae_vert(1)=deltae_vert(2)=deltae_vert(3)=deltae_vert(4). 
 
! !OUTPUT PARAMETERS:            
      real(8), intent(out) :: weight_vert(4)   ! the weight on the whole tetrahedron.

! !LOCAL VARIABLES:

      integer(4) :: j,k
      integer(4) :: ivert

      real(8)    :: aa, bb, cc, dd, bb1, bb3,bb4,omg2
      real(8)    :: vp


      integer    :: nmax_taylor = 4     ! maximum expansion order to be tried 
      real(8)    :: tol_taylor = 1.e-3  ! when the expansion variable is smaller than this number, 
                                        ! a Taylor expansion will be used to calculate weights
      real(8)    :: conv_taylor = 1.e-5 ! the expansion stops when the change is smaller than this 
       

      real(8), dimension(4) :: ev,tans, logs
      real(8), dimension(4,4) :: vdif
      character(20):: sname="stweight_imag"
      logical:: lwarn =.false.

! 
! !SYSTEM ROUTINES:

      intrinsic datan
      intrinsic dlog
      intrinsic dsign
 
! !REVISION HISTORY:
!
! Created 05.11.2004 by XZL.
!

!EOP
!BOC
      weight_vert(1:4) = 0.0d0

      select case(equiv_flag)

      case(4)                   ! for the case none of them are equal

        do ivert=1,4
          do j=1,4
            k=mod(j+ivert-2,4)+1       
            ev(k)=deltae_vert(j)
          enddo
          call sub_set_vdif_tans_logs  
          
          aa=2.0d0*(omeg**2-ev(1)**2)*vdif(1,2)*vdif(1,3)*vdif(1,4)*    &
     &      vdif(2,3)*vdif(2,4)*vdif(3,4) 

          !! A0
          omg2=omeg**2
     &    dd= 3.0d0*ev(1)**4                                            &
     &      -omg2*(ev(3)*ev(4)+ev(2)*(ev(3)+ev(4)))                     &
     &      -3.0*ev(1)**2*(omg2+ev(3)*ev(4)+ev(2)*(ev(3)+ev(4)))        &
     &      +2.0*ev(1)*( omg2*(ev(3)+ev(4))+ev(2)*(omg2+3.0*ev(3)*ev(4)))
          dd = dd * 2.0d0*omeg
     
          aa=aa+dd*vdif(2,3)*vdif(2,4)*vdif(3,4)*tans(1)
       
          !! A1
          aa=aa+2.0d0*omeg*(omeg**2-3.0d0*ev(2)**2)*vdif(1,3)**2*          &
     &       vdif(1,4)**2*vdif(3,4)*tans(2) 
        
          !! A2
          aa=aa-2.0d0*omeg*(omeg**2-3.0d0*ev(3)**2)*vdif(1,2)**2*          &
     &       vdif(1,4)**2*vdif(2,4)*tans(3) 
     
          !! A3
          aa=aa+2.0d0*omeg*(omeg**2-3.0d0*ev(4)**2)*vdif(1,2)**2*          &
     &       vdif(1,3)**2*vdif(2,3)*tans(4) 
     
          !! B0
          dd= -3.0d0*omeg**2*ev(2)*ev(3)*ev(4)                          &
     &      + ev(1)**4*(ev(2)+ev(3)+ev(4))                              &
     &      -2.0*ev(1)**3*(3.0*omg2+ev(3)*ev(4)+ev(2)*(ev(3)+ev(4)))    &
     &      +3.0*ev(1)**2*(omg2*(ev(2)+ev(3)+ev(4))+ev(2)*ev(3)*ev(4))

          bb=dd*vdif(2,3)*vdif(2,4)*vdif(3,4)*logs(1) 
        
          !! B1
          bb=bb+ev(2)*(3.0d0*omg2-ev(2)**2)*vdif(1,3)**2*vdif(1,4)**2*  &
     &     vdif(3,4)*logs(2)
     
          !! B2
          bb=bb+ev(3)*(-3.0d0*omg2+ev(3)**2)*vdif(1,2)**2*vdif(1,4)**2*  &
     &      vdif(2,4)*logs(3) 
     
          !! B3
          bb=bb-ev(4)*(-3.0d0*omg2+ev(4)**2)*vdif(1,2)**2*vdif(1,3)**2*  &
     &      vdif(2,3)*logs(4)
     
          !! commond denominator 
          cc=6.0d0*vdif(1,2)**2*vdif(1,3)**2*vdif(1,4)**2*vdif(2,3)*      &
     &       vdif(2,4)*vdif(3,4)
     
          call sub_set_weight(ivert)

        enddo ! ivert  
     
      case(6)                ! for the case when ev(1)=ev(1)
        do ivert=1,2
          ev(1)=deltae_vert(ivert)
          ev(2)=ev(1)
          ev(3:4)=deltae_vert(3:4)
          call sub_set_vdif_tans_logs  
          
          dd=ev(1)**3-2.0d0*omeg**2*(ev(3)+ev(4))-3.0d0*ev(1)**2*(ev(3)+ev(4))+ &
     &       ev(1)*(4.0d0*omeg**2+5.0d0*ev(3)*ev(4))   
        
          aa=vdif(3,1)*vdif(1,4)*vdif(3,4)*dd
        
          dd=-omeg**2*(3.0d0*ev(1)**2+ev(3)**2+ev(3)*ev(4)+ev(4)**2-3.0d0*ev(1)*&
     &       (ev(3)+ev(4)))+3.0d0*(-3.0d0*ev(1)**2*ev(3)*ev(4)+ev(3)**2*ev(4)**2+&
     &       ev(1)**3*(ev(3)+ev(4)))   
          aa=aa-2.0d0*omeg*vdif(3,4)*dd*tans(1)
          
          aa=aa-2.0d0*omeg*(omeg**2-3.0d0*ev(3)**2)*vdif(1,4)**3*tans(3)

          aa=aa+2.0d0*omeg*vdif(1,3)**3*(omeg**2-3.0d0*ev(4)**2)*tans(4)

          dd=-3.0d0*omeg**2*ev(3)*ev(4)*(ev(3)+ev(4))-3.0d0*ev(1)**2*ev(3)*    &
     &     ev(4)*(ev(3)+ev(4))+3.0d0*ev(1)*ev(3)*ev(4)*(3.0d0*omeg**2+ev(3)*  &
     &     ev(4))+ev(1)**3*(-3.0d0*omeg**2+ev(3)**2+ev(3)*ev(4)+ev(4)**2)  

          bb1=-vdif(3,4)*dd
          
          bb3=ev(3)*vdif(1,4)**3*(ev(3)**2-3.0d0*omeg**2)
     
          bb4=-ev(4)*vdif(1,3)**3*(ev(4)**2-3.0d0*omeg**2)
 
          bb=bb1*logs(1)+bb3*logs(3)+bb4*logs(4)
      
          cc=6.0d0*vdif(1,3)**3*vdif(1,4)**3*vdif(3,4)
          call sub_set_weight(ivert)
        enddo ! ivert  

        do ivert=3,4
          ev(1)=(deltae_vert(1)+deltae_vert(2))*0.5d0
          ev(2)=ev(1)
          do j=3,4
            k=mod(j+ivert,2)+3
            ev(k)=deltae_vert(j)
          enddo
          call sub_set_vdif_tans_logs  

          dd=-ev(1)*ev(3)**2+omeg**2*(ev(1)+ev(3)-2.0d0*ev(4))+ev(3)**2* & 
     &      ev(4)+ev(1)**2*vdif(4,3)
          aa=2.0d0*dd*vdif(1,3)*vdif(1,4)*vdif(3,4)
          
          dd=3.0d0*ev(1)**3+3.0d0*ev(1)**2*ev(3)+omeg**2*(ev(3)+2.0d0*   &
     &      ev(4))-3.0d0*ev(1)*(omeg**2+2.0d0*ev(3)*ev(4))
          aa=aa+2.0d0*vdif(3,4)**2*omeg*dd*tans(1)

          dd=3.0d0*ev(3)*(ev(3)**2+ev(1)*(ev(3)-2.0d0*ev(4)))+omeg**2*   &
     &      (ev(1)-3.0d0*ev(3)+2.0d0*ev(4))    
          aa=aa-2.0d0*vdif(1,4)**2*omeg*dd*tans(3)

          aa=aa+2.0d0*vdif(1,3)**3*omeg*(omeg**2-3.0d0*ev(4)**2)*        &
     &      tans(4)

          dd=ev(1)**2*(-3.0d0*ev(3)*ev(4)+ev(1)*(2.0d0*ev(3)+ev(4)))
          dd=dd+omeg**2*(-6.0d0*ev(1)**2+3.0d0*ev(1)*ev(4)+3.0d0*ev(3)*  &
     &      ev(4))
          bb=vdif(3,4)**2*dd*logs(1)

          dd=ev(3)**2*(2.0d0*ev(1)*ev(3)-3.0d0*ev(1)*ev(4)+ev(3)*ev(4))
          dd=dd+3.0d0*omeg**2*(-2.0d0*ev(3)**2+(ev(1)+ev(3))*ev(4))
        
          bb=bb-vdif(1,4)**2*dd*logs(3)

          bb=bb-ev(4)*vdif(1,3)**3*(ev(4)**2-3.0d0*omeg**2)*             &
     &      logs(4)
         
          cc=6.0d0*vdif(1,3)**3*vdif(1,4)**2*vdif(3,4)**2

          call sub_set_weight(ivert)
        enddo ! ivert  

      case(8)          !for the case when ev(1)=ev(1) and ev(3)=ev(4)
        ev(1)=(deltae_vert(1)+deltae_vert(2))/2.d0
        ev(2)=ev(1)
        ev(3)=(deltae_vert(3)+deltae_vert(4))/2.d0
        ev(4)=ev(3) 
        call sub_set_vdif_tans_logs  

        dd=6.0d0*omeg**2+ev(1)**2-5.0d0*ev(1)*ev(3)-2.0d0*ev(3)**2
        aa=vdif(3,1)*dd
  
        dd=6.0d0*omeg*(omeg**2 - ev(3)**2 - 2.0d0*ev(1)*ev(3))
        aa=aa+dd*(tans(1)-tans(3))

        bb=3.0d0*((ev(1)+2.0d0*ev(3))*omeg**2 - ev(1)*ev(3)**2 )
        bb=bb*(logs(1)-logs(3))

        cc=6.0d0*vdif(1,3)**4

        x=vdif(3,1)/sqrt(ev(1)**2+omeg**2)
        if(abs(x) < tol_taylor ) then 
          call sub_taylor_8(ev(1),vdif(3,1),x,weight_vert(1))
        else 
          call sub_set_weight(1)  
        endif 
        weight_vert(2) = weight_vert(1)

        dd=6.0d0*omeg**2-2.0d0*ev(1)**2-5.0d0*ev(1)*ev(3)+ev(3)**2
        aa=vdif(1,3)*dd
          
        dd=6.0d0*omeg*(omeg**2 - ev(1)**2 - 2.0d0*ev(1)*ev(3))
        aa=aa+dd*(tans(3)-tans(1))

        dd=3.0d0*( (2.0d0*ev(1) + ev(3))*omeg**2 - ev(1)**2*ev(3) )
        bb=dd*(logs(3)-logs(1))
         
        cc=6.0d0*vdif(1,3)**4
         
        x=vdif(1,3)/sqrt(ev(3)**2+omeg**2)
        if(abs(x) < tol_taylor ) then
          call sub_taylor_8(ev(3),vdif(1,3),x,weight_vert(3))
        else
          call sub_set_weight(3)
        endif
        weight_vert(4) = weight_vert(3)

      case(10)             ! for the case when ev(1)=ev(1)=ev(3)
        ev(1:3)=sum(deltae_vert(1:3))/3.0d0
        ev(4)=deltae_vert(4) 
        call sub_set_vdif_tans_logs  

        aa=vdif(1,4)*(6.0d0*omeg**2-2.0d0*ev(1)**2+7.0d0*ev(1)*ev(4)- &
     &       11.0d0*ev(4)**2)     
          
        dd=6.0d0*omeg*(omeg**2-3.0d0*ev(4)**2)
        aa=aa-dd*(tans(1)-tans(4))
          
        dd=3.0d0*ev(4)*(ev(4)**2-3.0d0*omeg**2)
        bb=dd*(logs(1)-logs(4))
         
        cc=18.0d0*vdif(1,4)**4

        call sub_set_weight(1) 
        weight_vert(2:3)=weight_vert(1) 

        dd=6.0d0*omeg**2+ev(1)**2-5.0d0*ev(1)*ev(4)-2.0d0*ev(4)**2
        aa=vdif(4,1)*dd
          
        dd=6.0d0*omeg*(omeg**2-ev(4)**2-2.0d0*ev(1)*ev(4))
        aa=aa+dd*(tans(1)-tans(4))
        
        dd=-3.0d0*ev(1)*ev(4)**2+3.0d0*omeg**2*(ev(1)+2.0d0*ev(4))
        bb=dd*(logs(1)-logs(4))
         
        cc=6.0d0*vdif(1,4)**4

        call sub_set_weight(4)

      case(16)
        ev(1:4)=sum(deltae_vert(1:4))/4.0
        aa=-ev(1)
        bb=0.0d0
        cc=12.0d0*(omeg**2+ev(1)**2)
        call sub_set_weight(1) 
        weight_vert(2:4)=weight_vert(1)
        
       case default
         write(6,*)'ERROR in stweight_imat'
         write(6,*)'wrong equiv_flag: ',equiv_flag
         stop "ERROR in stweight_imat"

      end select

      contains

        subroutine sub_set_weight(iv) 
        integer,intent(in):: iv

        if((dabs(cc)*weighttol).gt.dabs(aa+bb)) then
          weight_vert(iv)=(aa+bb)/cc
        else  
          weight_vert(iv)=0.d0
          write(fout,1)iv,(aa+bb)/cc,equiv_flag,omeg, &
     &                  deltae_vert,vdif,aa,bb,cc,vol_small_tetra
        endif
    1   format('WARNING in stweight_imag: unexpected big weight!',/ &
     &    'vertix:',i4,&
     &    ' weightt =',g16.6,' case:',i4,/,'omeg =',g16.6,/, &
     &    'deltae_vert =',4g16.6,    &
     &    /,'vdif = ',/,4(4g16.6,/),' aa =',g16.6,' bb =',g16.6,      &
     &    ' cc =',g16.6,' vol_small_tetra=',g16.6)     
        end subroutine 
 
        subroutine sub_set_vdif_tans_logs 
        implicit none
        integer(4) :: ii, jj
        do ii=1,4
          do jj=1,4
            vdif(ii,jj)=ev(ii)-ev(jj)
          enddo
          if(omeg.gt.1.0e-20)then
            tans(ii)=datan(ev(ii)/omeg)
          else
            tans(ii)= dsign(1.0d0,ev(ii))*2.0d0*datan(1.0d0)
          endif  
          if((omeg**2+ev(ii)**2).gt.1.0d-10)then
            logs(ii)=dlog(ev(ii)**2+omeg**2)
          else 
            logs(ii)=0.0d0
          endif
        enddo ! ii  
        end subroutine sub_set_vdif_tans_logs 

        subroutine sub_taylor_8(d0,omg,x,wt) 
!
!       Calculate weight by Talyor expansion for case 8
!
        real(8),intent(in)  :: d0,omg,x
        real(8),intent(out) :: wt
 
        integer:: n       
        real(8):: sqd2w2,wt_old,coef

        sqd2w2= sqrt(d0**2+omg**2)
        wt = - d0/(12.d0*sqd2w2**2)
        do n=1,nmax_taylor 
          wt_old=wt 
          if(n.eq.1) then 
            coef = (d0**2-omg**2)/(30.d0*sqd2w2**3)
          elseif(n.eq.2) then  
            coef = (-d0**3+3.d0*d0*omg**2)/(60.d0*sqd2w2**4)
          elseif(n.eq.3) then 
            coef = (d0**4-6.0*d0**2*omg**2+omg**4)/(105.d0*sqd2w2**5) 
          elseif(n.eq.4) then 
            coef = (-d0**5+10.0*d0**3*omg**2-5.d0*d0*omg**4) &
     &               /(168.d0*sqd2w2**6)
          endif 
          wt = wt + coef*x**n
          if(abs(wt-wt_old).lt.conv_taylor) exit 
        enddo 
        if(n.gt.nmax_taylor) then 
          write(fout,80) d0,omg,x,wt,wt_old 
        endif 
 80     format("WARNING from sub_taylor_8: Convergence not reached",/, &
     &         '  d0=',g16.6,' omg=',omg,/, &
     &         '  final wt=',wt, ' last wt=',wt_old,/) 
        end subroutine 
      end subroutine stweight_imag
!EOC      
