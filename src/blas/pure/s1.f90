module blas_pure_s1
  use blas_kinds, only : ik, sp
  implicit none
  private :: ik, sp

  interface
    pure function sasum(n, sx, incx)
      import :: ik, sp
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: sx(*)
      integer(ik), intent(in) :: incx
      real(sp) :: sasum
    end function sasum
  end interface
  interface
    pure subroutine saxpy(n, sa, sx, incx, sy, incy)
      import :: ik, sp
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: sa
      real(sp), intent(in) :: sx(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(inout) :: sy(*)
      integer(ik), intent(in) :: incy
    end subroutine saxpy
  end interface
  interface
    pure subroutine scopy(n, sx, incx, sy, incy)
      import :: ik, sp
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: sx(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(inout) :: sy(*)
      integer(ik), intent(in) :: incy
    end subroutine scopy
  end interface
  interface
    pure function sdot(n, sx, incx, sy, incy)
      import :: ik, sp
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: sx(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(in) :: sy(*)
      integer(ik), intent(in) :: incy
      real(sp) :: sdot
    end function sdot
  end interface
  interface
    pure function sdsdot(n, sb, sx, incx, sy, incy)
      import :: ik, sp
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: sb
      real(sp), intent(in) :: sx(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(in) :: sy(*)
      integer(ik), intent(in) :: incy
      real(sp) :: sdsdot
    end function sdsdot
  end interface
  interface
    pure function snrm2(n, x, incx)
      import :: ik, sp
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(sp) :: snrm2
    end function snrm2
  end interface
  interface
    pure subroutine srot(n, sx, incx, sy, incy, c, s)
      import :: ik, sp
      integer(ik), intent(in) :: n
      real(sp), intent(inout) :: sx(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(inout) :: sy(*)
      integer(ik), intent(in) :: incy
      real(sp), intent(in) :: c
      real(sp), intent(in) :: s
    end subroutine srot
  end interface
  interface
    pure subroutine srotg(a, b, c, s)
      import :: ik, sp
      real(sp), intent(in) :: a
      real(sp), intent(in) :: b
      real(sp), intent(out) :: c
      real(sp), intent(out) :: s
    end subroutine srotg
  end interface
  interface
    pure subroutine srotm(n, sx, incx, sy, incy, sparam)
      import :: ik, sp
      integer(ik), intent(in) :: n
      real(sp), intent(inout) :: sx(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(inout) :: sy(*)
      integer(ik), intent(in) :: incy
      real(sp), intent(in) :: sparam(5)
    end subroutine srotm
  end interface
  interface
    pure subroutine srotmg(sd1, sd2, sx1, sy1, sparam)
      import :: ik, sp
      real(sp), intent(inout) :: sd1
      real(sp), intent(inout) :: sd2
      real(sp), intent(inout) :: sx1
      real(sp), intent(in) :: sy1
      real(sp), intent(out) :: sparam(5)
    end subroutine srotmg
  end interface
  interface
    pure subroutine sscal(n, sa, sx, incx)
      import :: ik, sp
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: sa
      real(sp), intent(inout) :: sx(*)
      integer(ik), intent(in) :: incx
    end subroutine sscal
  end interface
  interface
    pure subroutine sswap(n, sx, incx, sy, incy)
      import :: ik, sp
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: sx(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(inout) :: sy(*)
      integer(ik), intent(in) :: incy
    end subroutine sswap
  end interface
  interface
    pure function isamax(n, sx, incx)
      import :: ik, sp
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: sx(*)
      integer(ik), intent(in) :: incx
      integer(ik) :: isamax
    end function isamax
  end interface
end module blas_pure_s1
