!! Copyright (C) 2003 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.

module liboct_parser
  implicit none

  ! Define the which routines can be seen from the outside
  private
  public :: loct_parse_init, loct_parse_putsym, loct_parse_input, loct_parse_end
  public :: loct_parse_isdef, loct_parse_int, loct_parse_float, loct_parse_cmplx, &
      loct_parse_string, loct_parse_logical
  public :: loct_parse_block_n, loct_parse_block_cols
  public :: loct_parse_block_int, loct_parse_block_float, loct_parse_block_cmplx, &
      loct_parse_block_string, loct_parse_block_logical
  public :: loct_parse_potential
  
  interface loct_parse_init
    integer function oct_parse_init(file_out)
      character(len=*), intent(in)  :: file_out
    end function oct_parse_init
  end interface

  interface loct_parse_putsym
    subroutine oct_parse_putsym_int(sym, i)
      character(len=*), intent(in)  :: sym
      integer, intent(in) :: i
    end subroutine oct_parse_putsym_int
    subroutine oct_parse_putsym_double(sym, d)
      character(len=*), intent(in)  :: sym
      real(8), intent(in) :: d
    end subroutine oct_parse_putsym_double
  end interface

  interface loct_parse_input
    integer function oct_parse_input(file_in)
      character(len=*), intent(in)  :: file_in
    end function oct_parse_input
  end interface

  interface loct_parse_end
    subroutine oct_parse_end()
    end subroutine oct_parse_end
  end interface

  interface loct_parse_isdef
    integer function oct_parse_isdef(name)
      character(len=*), intent(in) :: name
    end function oct_parse_isdef
  end interface

  interface loct_parse_int
    subroutine oct_parse_int(name, def, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: def
      integer, intent(out)         :: res
    end subroutine oct_parse_int
  end interface

  interface loct_parse_float
    subroutine oct_parse_double(name, def, res)
      character(len=*), intent(in)  :: name
      real(8),          intent(in)  :: def
      real(8),          intent(out) :: res
    end subroutine oct_parse_double
  end interface

  interface loct_parse_cmplx
    subroutine oct_parse_complex(name, def, res)
      character(len=*), intent(in) :: name
      complex(8), intent(in)       :: def
      complex(8), intent(out)      :: res
    end subroutine oct_parse_complex
  end interface

  interface loct_parse_string
    subroutine oct_parse_string(name, def, res)
      character(len=*), intent(in) :: name, def
      character(len=*), intent(out):: res
    end subroutine oct_parse_string
  end interface

  interface loct_parse_block_n
    integer function oct_parse_block_n(name)
      character(len=*), intent(in) :: name
    end function oct_parse_block_n
  end interface

  interface loct_parse_block_cols
    integer function oct_parse_block_cols(name, line)
      character(len=*), intent(in) :: name
      integer, intent(in) :: line
    end function oct_parse_block_cols
  end interface

  interface loct_parse_block_int
    subroutine oct_parse_block_int(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      integer, intent(out)         :: res
    end subroutine oct_parse_block_int

    subroutine oct_parse_block_int_new(name, l, c, res,def)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      integer, intent(in)         :: def
      integer, intent(out)         :: res
    end subroutine oct_parse_block_int_new
  end interface

  interface loct_parse_block_float
    subroutine oct_parse_block_double(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      real(8), intent(out)         :: res
    end subroutine oct_parse_block_double

    subroutine oct_parse_block_double_new(name, l, c, res,def)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      real(8), intent(in)         ::  def
      real(8), intent(out)         :: res
    end subroutine oct_parse_block_double_new
  end interface

  interface loct_parse_block_cmplx
    subroutine oct_parse_block_complex(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      complex(8), intent(out)      :: res
    end subroutine oct_parse_block_complex

    subroutine oct_parse_block_complex_new(name, l, c, res,def)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      complex(8), intent(out)      :: res
      complex(8), intent(in)      :: def
    end subroutine oct_parse_block_complex_new

  end interface

  interface loct_parse_block_string
    subroutine oct_parse_block_string(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      character(len=*), intent(out):: res
    end subroutine oct_parse_block_string

    subroutine oct_parse_block_string_new(name, l, c, res,def)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      character(len=*), intent(out):: res
      character(len=*), intent(in):: def
    end subroutine oct_parse_block_string_new
  end interface

  interface loct_parse_potential
    real(8) function oct_parse_potential(x, y, z, r, pot)
      real(8), intent(in) :: x, y, z, r
      character(len=*), intent(in) :: pot
    end function oct_parse_potential
  end interface

contains

  ! logical is a FORTRAN type, so we emulate the routine with integers
  subroutine loct_parse_logical(name, def, res)
    character(len=*), intent(IN) :: name
    logical, intent(in) :: def
    logical, intent(out) :: res
    
    integer :: idef, ires
    
    idef = 0
    if(def) idef = 1
    
    call oct_parse_int(name, idef, ires)
    res = (ires .ne. 0)
    
  end subroutine loct_parse_logical
  
  subroutine loct_parse_block_logical(name, l, c, res)
    character(len=*), intent(IN) :: name
    integer, intent(in)          :: l, c
    logical, intent(out)         :: res
    
    integer :: ires, r
    
    call oct_parse_block_int(name, l, c, ires)
    res = (ires .ne. 0)
    
  end subroutine loct_parse_block_logical

end module liboct_parser
