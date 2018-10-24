      subroutine sdef(iio,nsym,lattic)
!
!    redefines symmetry operation to include cxz and cyz
!    bravais lattices (from cxy)
!
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      character*4 lattic
      dimension iio(3,3,48)
      IF(LATTIc(1:3).EQ.'CXZ'.or.LATTIc(1:3).EQ.'BO ') then
      do 2 j=1,nsym
      ihelp=iio(2,2,j)                                   
      iio(2,2,j)=iio(3,3,j)
      iio(3,3,j)=ihelp
      ihelp=iio(1,2,j)                                   
      iio(1,2,j)=iio(1,3,j)
      iio(1,3,j)=ihelp
      ihelp=iio(2,3,j)                                   
      iio(2,3,j)=iio(3,2,j)
      iio(3,2,j)=ihelp
      ihelp=iio(2,1,j)                                   
      iio(2,1,j)=iio(3,1,j)
      iio(3,1,j)=ihelp
   2  continue
      else IF(LATTIc(1:3).EQ.'CYZ'.or.LATTIc(1:3).EQ.'AO ') then
      do 1 j=1,nsym
      ihelp=iio(1,1,j)                                   
      iio(1,1,j)=iio(3,3,j)
      iio(3,3,j)=ihelp
      ihelp=iio(1,2,j)                                   
      iio(1,2,j)=iio(3,2,j)
      iio(3,2,j)=ihelp
      ihelp=iio(1,3,j)                                   
      iio(1,3,j)=iio(3,1,j)
      iio(3,1,j)=ihelp
      ihelp=iio(2,3,j)                                   
      iio(2,3,j)=iio(2,1,j)
      iio(2,1,j)=ihelp
   1  continue
      end if
      return
      end
