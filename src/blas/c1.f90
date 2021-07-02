module blas_c1
  use blas_kinds, only : ik, sp
  implicit none
  private :: ik, sp

  interface
    subroutine caxpy(n, ca, cx, incx, cy, incy)
      import :: ik, sp
      integer(ik) :: n
      complex(sp) :: ca
      complex(sp) :: cx(*)
      integer(ik) :: incx
      complex(sp) :: cy(*)
      integer(ik) :: incy
    end subroutine caxpy
  end interface
  interface
    subroutine ccopy(n, cx, incx, cy, incy)
      import :: ik, sp
      integer(ik) :: n
      complex(sp) :: cx(*)
      integer(ik) :: incx
      complex(sp) :: cy(*)
      integer(ik) :: incy
    end subroutine ccopy
  end interface
  interface
    function cdotc(n, cx, incx, cy, incy)
      import :: ik, sp
      integer(ik) :: n
      complex(sp) :: cx(*)
      integer(ik) :: incx
      complex(sp) :: cy(*)
      integer(ik) :: incy
      complex(sp) :: cdotc
    end function cdotc
  end interface
  interface
    function cdotu(n, cx, incx, cy, incy)
      import :: ik, sp
      integer(ik) :: n
      complex(sp) :: cx(*)
      integer(ik) :: incx
      complex(sp) :: cy(*)
      integer(ik) :: incy
      complex(sp) :: cdotu
    end function cdotu
  end interface
  interface
    subroutine crotg(a, b, c, s)
      import :: ik, sp
      complex(sp) :: a
      complex(sp) :: b
      real(sp) :: c
      complex(sp) :: s
    end subroutine crotg
  end interface
  interface
    subroutine cscal(n, ca, cx, incx)
      import :: ik, sp
      integer(ik) :: n
      complex(sp) :: ca
      complex(sp) :: cx(*)
      integer(ik) :: incx
    end subroutine cscal
  end interface
  interface
    subroutine csrot(n, cx, incx, cy, incy, c, s)
      import :: ik, sp
      integer(ik) :: n
      complex(sp) :: cx(*)
      integer(ik) :: incx
      complex(sp) :: cy(*)
      integer(ik) :: incy
      real(sp) :: c
      real(sp) :: s
    end subroutine csrot
  end interface
  interface
    subroutine csscal(n, sa, cx, incx)
      import :: ik, sp
      integer(ik) :: n
      real(sp) :: sa
      complex(sp) :: cx(*)
      integer(ik) :: incx
    end subroutine csscal
  end interface
  interface
    subroutine cswap(n, cx, incx, cy, incy)
      import :: ik, sp
      integer(ik) :: n
      complex(sp) :: cx(*)
      integer(ik) :: incx
      complex(sp) :: cy(*)
      integer(ik) :: incy
    end subroutine cswap
  end interface
  interface
    function scabs1(z)
      import :: ik, sp
      complex(sp) :: z
      real(sp) :: scabs1
    end function scabs1
  end interface
  interface
    function scasum(n, cx, incx)
      import :: ik, sp
      integer(ik) :: n
      complex(sp) :: cx(*)
      integer(ik) :: incx
      real(sp) :: scasum
    end function scasum
  end interface
  interface
    function scnrm2(n, x, incx)
      import :: ik, sp
      integer(ik) :: n
      complex(sp) :: x(*)
      integer(ik) :: incx
      real(sp) :: scnrm2
    end function scnrm2
  end interface
  interface
    function icamax(n, cx, incx)
      import :: ik, sp
      integer(ik) :: n
      complex(sp) :: cx(*)
      integer(ik) :: incx
      integer(ik) :: icamax
    end function icamax
  end interface

end module blas_c1
