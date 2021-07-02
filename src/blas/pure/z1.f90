module blas_pure_z1
  use blas_kinds, only : ik, dp
  implicit none
  private :: ik, dp

  interface
    pure subroutine zaxpy(n, za, zx, incx, zy, incy)
      import :: ik, dp
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: za
      complex(dp), intent(in) :: zx(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(inout) :: zy(*)
      integer(ik), intent(in) :: incy
    end subroutine zaxpy
  end interface
  interface
    pure subroutine zcopy(n, zx, incx, zy, incy)
      import :: ik, dp
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: zx(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(inout) :: zy(*)
      integer(ik), intent(in) :: incy
    end subroutine zcopy
  end interface
  interface
    pure function zdotc(n, zx, incx, zy, incy)
      import :: ik, dp
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: zx(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(in) :: zy(*)
      integer(ik), intent(in) :: incy
      complex(dp) :: zdotc
    end function zdotc
  end interface
  interface
    pure function zdotu(n, zx, incx, zy, incy)
      import :: ik, dp
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: zx(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(in) :: zy(*)
      integer(ik), intent(in) :: incy
      complex(dp) :: zdotu
    end function zdotu
  end interface
  interface
    pure subroutine zdrot(n, zx, incx, zy, incy, c, s)
      import :: ik, dp
      integer(ik), intent(in) :: n
      complex(dp), intent(inout) :: zx(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(inout) :: zy(*)
      integer(ik), intent(in) :: incy
      real(dp), intent(in) :: c
      real(dp), intent(in) :: s
    end subroutine zdrot
  end interface
  interface
    pure subroutine zdscal(n, da, zx, incx)
      import :: ik, dp
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: da
      complex(dp), intent(inout) :: zx(*)
      integer(ik), intent(in) :: incx
    end subroutine zdscal
  end interface
  interface
    pure subroutine zrotg(a, b, c, s)
      import :: ik, dp
      complex(dp), intent(in) :: a
      complex(dp), intent(in) :: b
      real(dp), intent(out) :: c
      complex(dp), intent(out) :: s
    end subroutine zrotg
  end interface
  interface
    pure subroutine zscal(n, za, zx, incx)
      import :: ik, dp
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: za
      complex(dp), intent(inout) :: zx(*)
      integer(ik), intent(in) :: incx
    end subroutine zscal
  end interface
  interface
    pure subroutine zswap(n, zx, incx, zy, incy)
      import :: ik, dp
      integer(ik), intent(in) :: n
      complex(dp), intent(inout) :: zx(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(inout) :: zy(*)
      integer(ik), intent(in) :: incy
    end subroutine zswap
  end interface
  interface
    pure function dzasum(n, zx, incx)
      import :: ik, dp
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: zx(*)
      integer(ik), intent(in) :: incx
      real(dp) :: dzasum
    end function dzasum
  end interface
  interface
    pure function dznrm2(n, x, incx)
      import :: ik, dp
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(dp) :: dznrm2
    end function dznrm2
  end interface
  interface
    pure function dcabs1(z)
      import :: ik, dp
      complex(dp), intent(in) :: z
      real(dp) :: dcabs1
    end function dcabs1
  end interface
  interface
    pure function izamax(n, zx, incx)
      import :: ik, dp
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: zx(*)
      integer(ik), intent(in) :: incx
      integer(ik) :: izamax
    end function izamax
  end interface
end module blas_pure_z1
