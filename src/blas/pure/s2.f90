module blas_pure_s2
  use blas_kinds, only : ik, sp
  implicit none
  private :: ik, sp

  interface
    pure subroutine sgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: kl
      integer(ik), intent(in) :: ku
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(in) :: beta
      real(sp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine sgbmv
  end interface
  interface
    pure subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(in) :: beta
      real(sp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine sgemv
  end interface
  interface
    pure subroutine sger(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(sp), intent(inout) :: a(lda, *)
    end subroutine sger
  end interface
  interface
    pure subroutine ssbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(in) :: beta
      real(sp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine ssbmv
  end interface
  interface
    pure subroutine sspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
      import :: ik, sp
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: ap(*)
      real(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(in) :: beta
      real(sp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine sspmv
  end interface
  interface
    pure subroutine sspr2(uplo, n, alpha, x, incx, y, incy, ap)
      import :: ik, sp
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(sp), intent(inout) :: ap(*)
    end subroutine sspr2
  end interface
  interface
    pure subroutine sspr(uplo, n, alpha, x, incx, ap)
      import :: ik, sp
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(inout) :: ap(*)
    end subroutine sspr
  end interface
  interface
    pure subroutine ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(in) :: beta
      real(sp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine ssymv
  end interface
  interface
    pure subroutine ssyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(sp), intent(inout) :: a(lda, *)
    end subroutine ssyr2
  end interface
  interface
    pure subroutine ssyr(uplo, n, alpha, x, incx, a, lda)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(inout) :: a(lda, *)
    end subroutine ssyr
  end interface
  interface
    pure subroutine stbmv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine stbmv
  end interface
  interface
    pure subroutine stbsv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine stbsv
  end interface
  interface
    pure subroutine stpmv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, sp
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: ap(*)
      real(sp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine stpmv
  end interface
  interface
    pure subroutine stpsv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, sp
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: ap(*)
      real(sp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine stpsv
  end interface
  interface
    pure subroutine strmv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine strmv
  end interface
  interface
    pure subroutine strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, sp
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: side
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: transa
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(inout) :: b(ldb, *)
    end subroutine strsm
  end interface
  interface
    pure subroutine strsv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine strsv
  end interface

end module blas_pure_s2
