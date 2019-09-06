!BOP
! !ROUTINE:  dfactr
! !INTERFACE:
real(8) function dfactr(n)
! !INPUT/OUTPUT PARAMETERS:
!   n : numerator (in,integer)
! !DESCRIPTION:
!   Returns the double factor $n!!$.
!   n can be 0, positive number or negative odd number
!   No n checking.
!
! !REVISION HISTORY:
!   Created October 2019 (MYZ)
!EOP
!BOC
  implicit none
  ! arguments
  integer, intent(in) :: n
  ! local variables
  integer i
  if ((n.lt.0).and.(mod(n,2).eq.0)) then
    write(*,*)
    write(*,'("Error(dfactr): odd n < 0 : ",I8)') n
    write(*,*)
    stop
  end if
  dfactr = 1.0d0
  i = n
  if (n<-1) then
    do while(i<-1)
      i = i + 2
      dfactr = dfactr * real(i,8)
    enddo
  elseif (n.le.2) then
    dfactr = real(n,8)
  else
    do while(i.ge.1)
      dfactr = dfactr * real(i,8)
      i = i - 2
    enddo
  endif
return
end function
!EOC
