module blas_pure_c1
  use blas_kinds, only : ik, sp
  implicit none
  private :: ik, sp

  interface
    pure subroutine caxpy(n, ca, cx, incx, cy, incy)
      import :: ik, sp
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: ca
      complex(sp), intent(in) :: cx(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(inout) :: cy(*)
      integer(ik), intent(in) :: incy
    end subroutine caxpy
  end interface
  interface
    pure subroutine ccopy(n, cx, incx, cy, incy)
      import :: ik, sp
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: cx(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(inout) :: cy(*)
      integer(ik), intent(in) :: incy
    end subroutine ccopy
  end interface
  interface
    pure function cdotc(n, cx, incx, cy, incy)
      import :: ik, sp
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: cx(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(in) :: cy(*)
      integer(ik), intent(in) :: incy
      complex(sp) :: cdotc
    end function cdotc
  end interface
  interface
    pure function cdotu(n, cx, incx, cy, incy)
      import :: ik, sp
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: cx(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(in) :: cy(*)
      integer(ik), intent(in) :: incy
      complex(sp) :: cdotu
    end function cdotu
  end interface
  interface
    pure subroutine crotg(a, b, c, s)
      import :: ik, sp
      complex(sp), intent(in) :: a
      complex(sp), intent(in) :: b
      real(sp), intent(out) :: c
      complex(sp), intent(out) :: s
    end subroutine crotg
  end interface
  interface
    pure subroutine cscal(n, ca, cx, incx)
      import :: ik, sp
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: ca
      complex(sp), intent(inout) :: cx(*)
      integer(ik), intent(in) :: incx
    end subroutine cscal
  end interface
  interface
    pure subroutine csrot(n, cx, incx, cy, incy, c, s)
      import :: ik, sp
      integer(ik), intent(in) :: n
      complex(sp), intent(inout) :: cx(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(inout) :: cy(*)
      integer(ik), intent(in) :: incy
      real(sp), intent(in) :: c
      real(sp), intent(in) :: s
    end subroutine csrot
  end interface
  interface
    pure subroutine csscal(n, sa, cx, incx)
      import :: ik, sp
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: sa
      complex(sp), intent(inout) :: cx(*)
      integer(ik), intent(in) :: incx
    end subroutine csscal
  end interface
  interface
    pure subroutine cswap(n, cx, incx, cy, incy)
      import :: ik, sp
      integer(ik), intent(in) :: n
      complex(sp), intent(inout) :: cx(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(inout) :: cy(*)
      integer(ik), intent(in) :: incy
    end subroutine cswap
  end interface
  interface
    pure function scabs1(z)
      import :: ik, sp
      complex(sp), intent(in) :: z
      real(sp) :: scabs1
    end function scabs1
  end interface
  interface
    pure function scasum(n, cx, incx)
      import :: ik, sp
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: cx(*)
      integer(ik), intent(in) :: incx
      real(sp) :: scasum
    end function scasum
  end interface
  interface
    pure function scnrm2(n, x, incx)
      import :: ik, sp
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(sp) :: scnrm2
    end function scnrm2
  end interface
  interface
    pure function icamax(n, cx, incx)
      import :: ik, sp
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: cx(*)
      integer(ik), intent(in) :: incx
      integer(ik) :: icamax
    end function icamax
  end interface

end module blas_pure_c1
