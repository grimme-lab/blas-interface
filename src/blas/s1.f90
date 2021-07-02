module blas_s1
  use blas_kinds, only : ik, sp
  implicit none
  private :: ik, sp

  interface
    function sasum(n, sx, incx)
      import :: ik, sp
      integer(ik) :: n
      real(sp) :: sx(*)
      integer(ik) :: incx
      real(sp) :: sasum
    end function sasum
  end interface
  interface
    subroutine saxpy(n, sa, sx, incx, sy, incy)
      import :: ik, sp
      integer(ik) :: n
      real(sp) :: sa
      real(sp) :: sx(*)
      integer(ik) :: incx
      real(sp) :: sy(*)
      integer(ik) :: incy
    end subroutine saxpy
  end interface
  interface
    subroutine scopy(n, sx, incx, sy, incy)
      import :: ik, sp
      integer(ik) :: n
      real(sp) :: sx(*)
      integer(ik) :: incx
      real(sp) :: sy(*)
      integer(ik) :: incy
    end subroutine scopy
  end interface
  interface
    function sdot(n, sx, incx, sy, incy)
      import :: ik, sp
      integer(ik) :: n
      real(sp) :: sx(*)
      integer(ik) :: incx
      real(sp) :: sy(*)
      integer(ik) :: incy
      real(sp) :: sdot
    end function sdot
  end interface
  interface
    function sdsdot(n, sb, sx, incx, sy, incy)
      import :: ik, sp
      integer(ik) :: n
      real(sp) :: sb
      real(sp) :: sx(*)
      integer(ik) :: incx
      real(sp) :: sy(*)
      integer(ik) :: incy
      real(sp) :: sdsdot
    end function sdsdot
  end interface
  interface
    function snrm2(n, x, incx)
      import :: ik, sp
      integer(ik) :: n
      real(sp) :: x(*)
      integer(ik) :: incx
      real(sp) :: snrm2
    end function snrm2
  end interface
  interface
    subroutine srot(n, sx, incx, sy, incy, c, s)
      import :: ik, sp
      integer(ik) :: n
      real(sp) :: sx(*)
      integer(ik) :: incx
      real(sp) :: sy(*)
      integer(ik) :: incy
      real(sp) :: c
      real(sp) :: s
    end subroutine srot
  end interface
  interface
    subroutine srotg(a, b, c, s)
      import :: ik, sp
      real(sp) :: a
      real(sp) :: b
      real(sp) :: c
      real(sp) :: s
    end subroutine srotg
  end interface
  interface
    subroutine srotm(n, sx, incx, sy, incy, sparam)
      import :: ik, sp
      integer(ik) :: n
      real(sp) :: sx(*)
      integer(ik) :: incx
      real(sp) :: sy(*)
      integer(ik) :: incy
      real(sp) :: sparam(5)
    end subroutine srotm
  end interface
  interface
    subroutine srotmg(sd1, sd2, sx1, sy1, sparam)
      import :: ik, sp
      real(sp) :: sd1
      real(sp) :: sd2
      real(sp) :: sx1
      real(sp) :: sy1
      real(sp) :: sparam(5)
    end subroutine srotmg
  end interface
  interface
    subroutine sscal(n, sa, sx, incx)
      import :: ik, sp
      integer(ik) :: n
      real(sp) :: sa
      real(sp) :: sx(*)
      integer(ik) :: incx
    end subroutine sscal
  end interface
  interface
    subroutine sswap(n, sx, incx, sy, incy)
      import :: ik, sp
      integer(ik) :: n
      real(sp) :: sx(*)
      integer(ik) :: incx
      real(sp) :: sy(*)
      integer(ik) :: incy
    end subroutine sswap
  end interface
  interface
    function isamax(n, sx, incx)
      import :: ik, sp
      integer(ik) :: n
      real(sp) :: sx(*)
      integer(ik) :: incx
      integer(ik) :: isamax
    end function isamax
  end interface
end module blas_s1
