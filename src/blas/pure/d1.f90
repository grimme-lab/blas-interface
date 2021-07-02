module blas_pure_d1
  use blas_kinds, only : ik, sp, dp
  implicit none
  private :: ik, sp, dp

  interface
    pure function dasum(n, dx, incx)
      import :: ik, dp
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: dx(*)
      integer(ik), intent(in) :: incx
      real(dp) :: dasum
    end function dasum
  end interface
  interface
    pure subroutine daxpy(n, da, dx, incx, dy, incy)
      import :: ik, dp
      integer(ik), intent(in) :: n
      real(dp) , intent(in):: da
      real(dp) , intent(in):: dx(*)
      integer(ik), intent(in) :: incx
      real(dp) , intent(inout):: dy(*)
      integer(ik), intent(in) :: incy
    end subroutine daxpy
  end interface
  interface
    pure subroutine dcopy(n, dx, incx, dy, incy)
      import :: ik, dp
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: dx(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(inout) :: dy(*)
      integer(ik), intent(in) :: incy
    end subroutine dcopy
  end interface
  interface
    pure function ddot(n, dx, incx, dy, incy)
      import :: ik, dp
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: dx(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(in) :: dy(*)
      integer(ik), intent(in) :: incy
      real(dp) :: ddot
    end function ddot
  end interface
  interface
    pure function dnrm2(n, x, incx)
      import :: ik, dp
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(dp) :: dnrm2
    end function dnrm2
  end interface
  interface
    pure subroutine drot(n, dx, incx, dy, incy, c, s)
      import :: ik, dp
      integer(ik), intent(in) :: n
      real(dp), intent(inout) :: dx(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(inout) :: dy(*)
      integer(ik), intent(in) :: incy
      real(dp), intent(in) :: c
      real(dp), intent(in) :: s
    end subroutine drot
  end interface
  interface
    pure subroutine drotg(a, b, c, s)
      import :: ik, dp
      real(dp), intent(in) :: a
      real(dp), intent(in) :: b
      real(dp), intent(out) :: c
      real(dp), intent(out) :: s
    end subroutine drotg
  end interface
  interface
    pure subroutine drotm(n, dx, incx, dy, incy, dparam)
      import :: ik, dp
      integer(ik), intent(in) :: n
      real(dp), intent(inout) :: dx(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(inout) :: dy(*)
      integer(ik), intent(in) :: incy
      real(dp), intent(in) :: dparam(5)
    end subroutine drotm
  end interface
  interface
    pure subroutine drotmg(dd1, dd2, dx1, dy1, dparam)
      import :: ik, dp
      real(dp), intent(inout) :: dd1
      real(dp), intent(inout) :: dd2
      real(dp), intent(inout) :: dx1
      real(dp), intent(in) :: dy1
      real(dp), intent(out) :: dparam(5)
    end subroutine drotmg
  end interface
  interface
    pure subroutine dscal(n, da, dx, incx)
      import :: ik, dp
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: da
      real(dp), intent(inout) :: dx(*)
      integer(ik), intent(in) :: incx
    end subroutine dscal
  end interface
  interface
    pure function dsdot(n, sx, incx, sy, incy)
      import :: ik, dp, sp
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: sx(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(in) :: sy(*)
      integer(ik), intent(in) :: incy
      real(dp) :: dsdot
    end function dsdot
  end interface
  interface
    pure subroutine dswap(n, dx, incx, dy, incy)
      import :: ik, dp
      integer(ik), intent(in) :: n
      real(dp), intent(inout) :: dx(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(inout) :: dy(*)
      integer(ik), intent(in) :: incy
    end subroutine dswap
  end interface
  interface
    pure function idamax(n, dx, incx)
      import :: ik, dp
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: dx(*)
      integer(ik), intent(in) :: incx
      integer(ik) :: idamax
    end function idamax
  end interface
end module blas_pure_d1
