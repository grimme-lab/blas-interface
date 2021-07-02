module blas_d1
  use blas_kinds, only : ik, sp, dp
  implicit none
  private :: ik, sp, dp

  interface
    function dasum(n, dx, incx)
      import :: ik, dp
      integer(ik) :: n
      real(dp) :: dx(*)
      integer(ik) :: incx
      real(dp) :: dasum
    end function dasum
  end interface
  interface
    subroutine daxpy(n, da, dx, incx, dy, incy)
      import :: ik, dp
      integer(ik) :: n
      real(dp) :: da
      real(dp) :: dx(*)
      integer(ik) :: incx
      real(dp) :: dy(*)
      integer(ik) :: incy
    end subroutine daxpy
  end interface
  interface
    function dcabs1(z)
      import :: ik, dp
      complex(dp) :: z
      real(dp) :: dcabs1
    end function dcabs1
  end interface
  interface
    subroutine dcopy(n, dx, incx, dy, incy)
      import :: ik, dp
      integer(ik) :: n
      real(dp) :: dx(*)
      integer(ik) :: incx
      real(dp) :: dy(*)
      integer(ik) :: incy
    end subroutine dcopy
  end interface
  interface
    function ddot(n, dx, incx, dy, incy)
      import :: ik, dp
      integer(ik) :: n
      real(dp) :: dx(*)
      integer(ik) :: incx
      real(dp) :: dy(*)
      integer(ik) :: incy
      real(dp) :: ddot
    end function ddot
  end interface
  interface
    function dnrm2(n, x, incx)
      import :: ik, dp
      integer(ik) :: n
      real(dp) :: x(*)
      integer(ik) :: incx
      real(dp) :: dnrm2
    end function dnrm2
  end interface
  interface
    subroutine drot(n, dx, incx, dy, incy, c, s)
      import :: ik, dp
      integer(ik) :: n
      real(dp) :: dx(*)
      integer(ik) :: incx
      real(dp) :: dy(*)
      integer(ik) :: incy
      real(dp) :: c
      real(dp) :: s
    end subroutine drot
  end interface
  interface
    subroutine drotg(a, b, c, s)
      import :: ik, dp
      real(dp) :: a
      real(dp) :: b
      real(dp) :: c
      real(dp) :: s
    end subroutine drotg
  end interface
  interface
    subroutine drotm(n, dx, incx, dy, incy, dparam)
      import :: ik, dp
      integer(ik) :: n
      real(dp) :: dx(*)
      integer(ik) :: incx
      real(dp) :: dy(*)
      integer(ik) :: incy
      real(dp) :: dparam(5)
    end subroutine drotm
  end interface
  interface
    subroutine drotmg(dd1, dd2, dx1, dy1, dparam)
      import :: ik, dp
      real(dp) :: dd1
      real(dp) :: dd2
      real(dp) :: dx1
      real(dp) :: dy1
      real(dp) :: dparam(5)
    end subroutine drotmg
  end interface
  interface
    subroutine dscal(n, da, dx, incx)
      import :: ik, dp
      integer(ik) :: n
      real(dp) :: da
      real(dp) :: dx(*)
      integer(ik) :: incx
    end subroutine dscal
  end interface
  interface
    function dsdot(n, sx, incx, sy, incy)
      import :: ik, dp, sp
      integer(ik) :: n
      real(sp) :: sx(*)
      integer(ik) :: incx
      real(sp) :: sy(*)
      integer(ik) :: incy
      real(dp) :: dsdot
    end function dsdot
  end interface
  interface
    subroutine dswap(n, dx, incx, dy, incy)
      import :: ik, dp
      integer(ik) :: n
      real(dp) :: dx(*)
      integer(ik) :: incx
      real(dp) :: dy(*)
      integer(ik) :: incy
    end subroutine dswap
  end interface
  interface
    function idamax(n, dx, incx)
      import :: ik, dp
      integer(ik) :: n
      real(dp) :: dx(*)
      integer(ik) :: incx
      integer(ik) :: idamax
    end function idamax
  end interface
end module blas_d1
