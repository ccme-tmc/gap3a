!BOP
!
! !ROUTINE: w2k_vxccub
!
! !INTERFACE:
       subroutine w2k_vxccub(iat,isp)
!
! !DESCRIPTION:
      
! calculates cubic harmonics after kara \& kurki-suonio
! acta cryst a 1981 37 201-210
! gkhm 2/5-01

! !USES:

      use xcpot,    only: lmxc,lxcm,vxclm
      use struk,    only: iatnr,nrpt

! !INPUT PARAMETERS:

       implicit none
       
       integer(4), intent(in) :: iat
       integer(4), intent(in) :: isp
       

! !OUTPUT PARAMETERS:


! !LOCAL VARIABLES:
       
       integer(4) :: i
       integer(4) :: irp
       integer(4) :: npt
       integer(4) :: j
       integer(4) :: lxc
       
       real(8) :: c1,c2,c3
       real(8) :: sq1,sq2
       real(8) :: c_kub(0:10,0:10)
       character(20) :: sname="w2k_vxccub"
       character(200):: msg

! !REVISION HISTORY:
!
! Taken from w2k_vxccub.f (LAPW5, Wien2k)
! Last modified: 18.08.05 by RGA
!
!EOP
!BOC
      npt=nrpt(iat)
      if(iatnr(iat).gt.0)then                                                      
        do i=1,10
          do j=1,10
            c_kub(i,j)=0.0d0
         enddo
        enddo
        c_kub(0,0)=1.d0
        c_kub(3,2)=1.d0
        c_kub(4,0)=.5d0*sqrt(7.d0/3.d0)
        c_kub(4,4)=.5*sqrt(5.d0/3.d0)
        c_kub(6,0)=.5d0*sqrt(.5d0)
        c_kub(6,2)=.25d0*sqrt(11.d0)
        c_kub(6,4)=-.5d0*sqrt(7.d0/2.d0)
        c_kub(6,6)=-.25d0*sqrt(5.d0)
        c_kub(7,2)=.5d0*sqrt(13.d0/6.d0)
        c_kub(7,6)=.5d0*sqrt(11.d0/16.d0)
        c_kub(8,0)=.125d0*sqrt(33.d0)
        c_kub(8,4)=.25d0*sqrt(7.d0/3.d0)
        c_kub(8,8)=.125d0*sqrt(65.d0/3.d0)
        c_kub(9,2)=.25d0*sqrt(3.d0)
        c_kub(9,4)=.5d0*sqrt(17.d0/6.d0)
        c_kub(9,6)=-.25d0*sqrt(13.d0)
        c_kub(9,8)=-.5d0*sqrt(7.d0/6.d0)
        c_kub(10,0)=.125*sqrt(65.d0/6.d0)
        c_kub(10,2)=.125*sqrt(247.d0/6.d0)
        c_kub(10,4)=-.25*sqrt(11.d0/2.d0)
        c_kub(10,6)=0.0625d0*sqrt(19.d0/3.d0)
        c_kub(10,8)=-.125*sqrt(187.d0/6.d0)
        c_kub(10,10)=-.0625d0*sqrt(85.d0)
        
        sq2=sqrt(2.0d0)

        lxc=1
        do while (lxc.le.lxcm(iat))
          select case (lmxc(1,lxc,iat))

          case(0)
            if(lmxc(2,lxc,iat).eq.0) then
              lxc=lxc+1
            else
              goto 991
            endif

          case(-3)      
            if(lmxc(2,lxc,iat).eq.2)then  
              do irp=1,npt
                vxclm(irp,lxc,iat,isp)=-vxclm(irp,lxc,iat,isp)/sq2
              enddo  
              lxc=lxc+1
            else
              goto 991
            endif
          
          case(4,6,-7,-9)
            c1=c_kub(abs(lmxc(1,lxc,iat)),lmxc(2,lxc,iat))
            c2=c_kub(abs(lmxc(1,lxc,iat)),lmxc(2,lxc,iat)+4)
            sq1=sq2
            if(lmxc(2,lxc,iat).eq.0)sq1=1.0d0
            do irp=1,npt
              vxclm(irp,lxc,iat,isp)=  vxclm(irp,lxc,iat,isp)*c1        &
     &                               + vxclm(irp,lxc+1,iat,isp)*c2
              vxclm(irp,lxc+1,iat,isp)=vxclm(irp,lxc,iat,isp)*c2/sq2
              vxclm(irp,lxc,iat,isp)=vxclm(irp,lxc,iat,isp)*c1/sq1
            enddo  
            lxc=lxc+2
            
          case(8,10)  
            sq1=sq2
            if(lmxc(2,lxc,iat).eq.0)sq1=1.0d0
            c1=c_kub(abs(lmxc(1,lxc,iat)),lmxc(2,lxc,iat))
            c2=c_kub(abs(lmxc(1,lxc,iat)),lmxc(2,lxc,iat)+4)
            c3=c_kub(abs(lmxc(1,lxc,iat)),lmxc(2,lxc,iat)+8)
            do irp=1,npt
              vxclm(irp,lxc,iat,isp)=vxclm(irp,lxc,iat,isp)  *c1        &
     &                             + vxclm(irp,lxc+1,iat,isp)*c2        &
     &                             + vxclm(irp,lxc+2,iat,isp)*c3
              vxclm(irp,lxc+1,iat,isp)=vxclm(irp,lxc,iat,isp)*c2/sq2
              vxclm(irp,lxc+2,iat,isp)=vxclm(irp,lxc,iat,isp)*c3/sq2
              vxclm(irp,lxc,iat,isp)  =vxclm(irp,lxc,iat,isp)*c1/sq1
            enddo  
            
            lxc=lxc+3

          case default
            goto 991
            
          end select  

        enddo ! lxc
      endif
      
      return
      
  991 write(msg,*) 'uncorrect lmxc list for cubic structure',&
     &            'l=',lmxc(1,lxc,iat), ' m=',lmxc(2,lxc,iat)
      call outerr(sname,msg)

      end subroutine w2k_vxccub
!EOC      
  
