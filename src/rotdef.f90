!BOP
! !ROUTINE: rotdef
!
! !INTERFACE:
      subroutine rotdef
!
! !DESCRIPTION:
!
!        Finds the symmetry operations, which transform equivalent
!        atoms into each other
!        the operation (ROTIJ(3,3,idf), TAU(3,idf)) transforms
!        a position of a "not equivalent" atom to the position of
!        a corresponding equivalent atom (idf).
!

! !USES:

      use constants, only: pi
      use struk,     only: gbas, imat, nsym, lattic, mult, nat, ortho,  &
     &                     pos, rbas, tau,rotij

!
! !LOCAL VARIABLES:
!
      implicit none

      integer(4) :: isym 
      integer(4) :: idf 
      integer(4) :: idf0 
      integer(4) :: iat 
      integer(4) :: i 
      integer(4) :: ieq
      integer(4) :: ncount
      
      real(8)    :: sym_pos(1:3)
      real(8)    :: delta_pos(1:3) 
      real(8)    :: delta_cpos(1:3)
     
      character(20):: sname="rotdef" 
      character(len=120) :: msg
      
      logical :: found_sym

! !DEFINED PARAMETERS:
!
      real(8), parameter :: one = 1.0d+0
      real(8), parameter :: half = 0.500000001d0
      real(8), parameter :: toler = 1.0d-4

! !REVISION HISTORY:
!
! Last modified Jan. 27th. 2004. (RGA)
! 
!
!EOP
!BOC
!
!        external subroutines
!
      external outerr
      external locdef
!
!        intrinsic functions
!
      intrinsic abs
      intrinsic mod
      intrinsic sum
      
      idf = 0
      ncount = 0
      do iat = 1, nat
        idf0 = idf + 1
        do ieq = 1, mult(iat)
          idf = idf + 1

          found_sym = .false.
sym_loop: do isym = 1, nsym
            do i = 1, 3
              sym_pos(i) = sum(imat(:,i,isym)*pos(:,idf0))
              sym_pos(i) = sym_pos(i) + tau(i,isym) + 1.0d+0
              sym_pos(i) = mod(sym_pos(i),one)
              delta_pos(i) = abs(sym_pos(i)-pos(i,idf))
            enddo  
            found_sym = (delta_pos(1).lt.toler).and.(delta_pos(2).lt.   &
     &           toler) .and. (delta_pos(3) .lt. toler)
            if (found_sym) exit sym_loop

            !!....check positions for centered lattices
            select case (lattic(1:1))
            case('B')
              do i = 1, 3
                delta_cpos(i)=mod(delta_pos(i)+half,one)
              enddo  
              found_sym = (delta_cpos(1) .lt. toler) .and. (delta_cpos(2)  &
     &                    .lt. toler) .and. (delta_cpos(3) .lt. toler)
              if (found_sym) exit sym_loop
             
            case('F')
              do i = 1, 3
                delta_cpos(i)=mod(delta_pos(i)+half,one)
              enddo
              do i=1,3
                delta_cpos(i)=delta_pos(i)  
                found_sym = (delta_cpos(1) .lt. toler) .and. (delta_cpos(2)  &
     &                    .lt. toler) .and. (delta_cpos(3) .lt. toler)
                if (found_sym) exit sym_loop
                delta_cpos(i)=mod(delta_pos(i)+half,one)
              enddo
              
            case('C')

              select case (lattic(2:3))
              case('XY')
                delta_cpos(1)=mod(delta_pos(1)+half,one)    
                delta_cpos(2)=mod(delta_pos(2)+half,one)    
                delta_cpos(3)=delta_pos(3)
                found_sym = (delta_cpos(1).lt.toler).and.(delta_cpos(2)  &
     &                  .lt. toler) .and. (delta_cpos(3) .lt. toler)
                if (found_sym) exit sym_loop
              case('XZ')
                delta_cpos(1)=mod(delta_pos(1)+half,one)    
                delta_cpos(2)=delta_pos(2)
                delta_cpos(3)=mod(delta_pos(3)+half,one)    
                found_sym = (delta_cpos(1).lt.toler).and.(delta_cpos(2)  &
     &                  .lt. toler) .and. (delta_cpos(3) .lt. toler)
                if (found_sym) exit sym_loop
              case('YZ')
                delta_cpos(1)=delta_pos(1)
                delta_cpos(2)=mod(delta_pos(2)+half,one)    
                delta_cpos(3)=mod(delta_pos(3)+half,one)
                found_sym = (delta_cpos(1).lt.toler).and.(delta_cpos(2)  &
     &                  .lt. toler) .and. (delta_cpos(3) .lt. toler)
                if (found_sym) exit sym_loop
              end select
            end select    
          enddo sym_loop

          if(.not.found_sym) then 
            write(6,*) "ERROR in rotdef: no sym operation found for iat=",&
   &           iat,"ieq=",ieq,"idf",idf
            stop "ERROR in rotdef"
          endif 

          ncount = ncount + 1
          rotij(1:3,1:3,idf) = imat(1:3,1:3,isym)
!
!     redefine rotation matrix for non-orthogonal case
!     for  mon.cxz type rotation skipped (caution!!)
!
          if((.not. ortho) .and. (lattic(1:3).ne.'CXZ'))  &
            call locdef(rbas,gbas,rotij(1:3,1:3,idf))
!
        enddo ! ieq
      enddo ! iat
!
      if (ncount .ne. idf) goto 910
!
      return

  910 write(msg,*) "rotij not defined for all atoms of basis.",'ncount=',ncount
      call outerr(sname,msg)

      end subroutine rotdef
!EOC
