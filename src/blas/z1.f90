module blas_z1
  use blas_kinds, only : ik, dp
  implicit none
  private :: ik, dp

  interface
    subroutine zaxpy(n, za, zx, incx, zy, incy)
      import :: ik, dp
      integer(ik) :: n
      complex(dp) :: za
      complex(dp) :: zx(*)
      integer(ik) :: incx
      complex(dp) :: zy(*)
      integer(ik) :: incy
    end subroutine zaxpy
  end interface
  interface
    subroutine zcopy(n, zx, incx, zy, incy)
      import :: ik, dp
      integer(ik) :: n
      complex(dp) :: zx(*)
      integer(ik) :: incx
      complex(dp) :: zy(*)
      integer(ik) :: incy
    end subroutine zcopy
  end interface
  interface
    function zdotc(n, zx, incx, zy, incy)
      import :: ik, dp
      integer(ik) :: n
      complex(dp) :: zx(*)
      integer(ik) :: incx
      complex(dp) :: zy(*)
      integer(ik) :: incy
      complex(dp) :: zdotc
    end function zdotc
  end interface
  interface
    function zdotu(n, zx, incx, zy, incy)
      import :: ik, dp
      integer(ik) :: n
      complex(dp) :: zx(*)
      integer(ik) :: incx
      complex(dp) :: zy(*)
      integer(ik) :: incy
      complex(dp) :: zdotu
    end function zdotu
  end interface
  interface
    subroutine zdrot(n, zx, incx, zy, incy, c, s)
      import :: ik, dp
      integer(ik) :: n
      complex(dp) :: zx(*)
      integer(ik) :: incx
      complex(dp) :: zy(*)
      integer(ik) :: incy
      real(dp) :: c
      real(dp) :: s
    end subroutine zdrot
  end interface
  interface
    subroutine zdscal(n, da, zx, incx)
      import :: ik, dp
      integer(ik) :: n
      real(dp) :: da
      complex(dp) :: zx(*)
      integer(ik) :: incx
    end subroutine zdscal
  end interface
  interface
    subroutine zrotg(a, b, c, s)
      import :: ik, dp
      complex(dp) :: a
      complex(dp) :: b
      real(dp) :: c
      complex(dp) :: s
    end subroutine zrotg
  end interface
  interface
    subroutine zscal(n, za, zx, incx)
      import :: ik, dp
      integer(ik) :: n
      complex(dp) :: za
      complex(dp) :: zx(*)
      integer(ik) :: incx
    end subroutine zscal
  end interface
  interface
    subroutine zswap(n, zx, incx, zy, incy)
      import :: ik, dp
      integer(ik) :: n
      complex(dp) :: zx(*)
      integer(ik) :: incx
      complex(dp) :: zy(*)
      integer(ik) :: incy
    end subroutine zswap
  end interface
  interface
    function dzasum(n, zx, incx)
      import :: ik, dp
      integer(ik) :: n
      complex(dp) :: zx(*)
      integer(ik) :: incx
      real(dp) :: dzasum
    end function dzasum
  end interface
  interface
    function dznrm2(n, x, incx)
      import :: ik, dp
      integer(ik) :: n
      complex(dp) :: x(*)
      integer(ik) :: incx
      real(dp) :: dznrm2
    end function dznrm2
  end interface
  interface
    function izamax(n, zx, incx)
      import :: ik, dp
      integer(ik) :: n
      complex(dp) :: zx(*)
      integer(ik) :: incx
      integer(ik) :: izamax
    end function izamax
  end interface
end module blas_z1
