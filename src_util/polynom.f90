!BOP
!
! !ROUTINE: polynom
!
! !INTERFACE: 
real(8) function polynom(m,np,xa,ya,c,x)

! !INPUT PARAMETERS:
implicit none

integer, intent(in) :: m
integer, intent(in) :: np
real(8), intent(in) :: xa(np)
real(8), intent(in) :: ya(np)
real(8), intent(out) :: c(np)
real(8), intent(in) :: x
! !LOCAL VARIABLES:
integer i,j
real(8) x1,t1,sum
!EOP
!BOC
if (np.le.0) then
  write(6,*)
  write(6,'("Error(polynom): np <= 0 : ",I8)') np
  write(6,*)
  call outerr("polynom","np<=0")
end if
if (m.gt.np) then
  polynom=0.d0
  return
end if
x1=xa(1)
! find the polynomial coefficients in divided differences form
c(:)=ya(:)
do i=2,np
  do j=np,i,-1
    c(j)=(c(j)-c(j-1))/(xa(j)-xa(j+1-i))
  end do
end do
! special case m=0
if (m.eq.0) then
  sum=c(1)
  t1=1.d0
  do i=2,np
    t1=t1*(x-xa(i-1))
    sum=sum+c(i)*t1
  end do
  polynom=sum
  return
end if
! convert to standard form
do j=1,np-1
  do i=1,np-j
    c(np-i)=c(np-i)-(xa(np-i-j+1)-x1)*c(np-i+1)
  end do
end do
if (m.gt.0) then
! take the m'th derivative
  do j=1,m
    do i=m+1,np
      c(i)=c(i)*dble(i-j)
    end do
  end do
  t1=c(np)
  do i=np-1,m+1,-1
    t1=t1*(x-x1)+c(i)
  end do
  polynom=t1
else
! find the integral
  t1=c(np)/dble(np)
  do i=np-1,1,-1
    t1=t1*(x-x1)+c(i)/dble(i)
  end do
  polynom=t1*(x-x1)
end if
return
end function
!EOC
