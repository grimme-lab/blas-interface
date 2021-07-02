module blas_s2
  use blas_kinds, only : ik, sp
  implicit none
  private :: ik, sp

  interface
    subroutine sgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: trans
      integer(ik) :: m
      integer(ik) :: n
      integer(ik) :: kl
      integer(ik) :: ku
      real(sp) :: alpha
      real(sp) :: a(lda, *)
      real(sp) :: x(*)
      integer(ik) :: incx
      real(sp) :: beta
      real(sp) :: y(*)
      integer(ik) :: incy
    end subroutine sgbmv
  end interface
  interface
    subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: trans
      integer(ik) :: m
      integer(ik) :: n
      real(sp) :: alpha
      real(sp) :: a(lda, *)
      real(sp) :: x(*)
      integer(ik) :: incx
      real(sp) :: beta
      real(sp) :: y(*)
      integer(ik) :: incy
    end subroutine sgemv
  end interface
  interface
    subroutine sger(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer(ik) :: lda
      integer(ik) :: m
      integer(ik) :: n
      real(sp) :: alpha
      real(sp) :: x(*)
      integer(ik) :: incx
      real(sp) :: y(*)
      integer(ik) :: incy
      real(sp) :: a(lda, *)
    end subroutine sger
  end interface
  interface
    subroutine ssbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      integer(ik) :: k
      real(sp) :: alpha
      real(sp) :: a(lda, *)
      real(sp) :: x(*)
      integer(ik) :: incx
      real(sp) :: beta
      real(sp) :: y(*)
      integer(ik) :: incy
    end subroutine ssbmv
  end interface
  interface
    subroutine sspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
      import :: ik, sp
      character(len=1) :: uplo
      integer(ik) :: n
      real(sp) :: alpha
      real(sp) :: ap(*)
      real(sp) :: x(*)
      integer(ik) :: incx
      real(sp) :: beta
      real(sp) :: y(*)
      integer(ik) :: incy
    end subroutine sspmv
  end interface
  interface
    subroutine sspr2(uplo, n, alpha, x, incx, y, incy, ap)
      import :: ik, sp
      character(len=1) :: uplo
      integer(ik) :: n
      real(sp) :: alpha
      real(sp) :: x(*)
      integer(ik) :: incx
      real(sp) :: y(*)
      integer(ik) :: incy
      real(sp) :: ap(*)
    end subroutine sspr2
  end interface
  interface
    subroutine sspr(uplo, n, alpha, x, incx, ap)
      import :: ik, sp
      character(len=1) :: uplo
      integer(ik) :: n
      real(sp) :: alpha
      real(sp) :: x(*)
      integer(ik) :: incx
      real(sp) :: ap(*)
    end subroutine sspr
  end interface
  interface
    subroutine ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      real(sp) :: alpha
      real(sp) :: a(lda, *)
      real(sp) :: x(*)
      integer(ik) :: incx
      real(sp) :: beta
      real(sp) :: y(*)
      integer(ik) :: incy
    end subroutine ssymv
  end interface
  interface
    subroutine ssyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      real(sp) :: alpha
      real(sp) :: x(*)
      integer(ik) :: incx
      real(sp) :: y(*)
      integer(ik) :: incy
      real(sp) :: a(lda, *)
    end subroutine ssyr2
  end interface
  interface
    subroutine ssyr(uplo, n, alpha, x, incx, a, lda)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      real(sp) :: alpha
      real(sp) :: x(*)
      integer(ik) :: incx
      real(sp) :: a(lda, *)
    end subroutine ssyr
  end interface
  interface
    subroutine stbmv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      integer(ik) :: k
      real(sp) :: a(lda, *)
      real(sp) :: x(*)
      integer(ik) :: incx
    end subroutine stbmv
  end interface
  interface
    subroutine stbsv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      integer(ik) :: k
      real(sp) :: a(lda, *)
      real(sp) :: x(*)
      integer(ik) :: incx
    end subroutine stbsv
  end interface
  interface
    subroutine stpmv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, sp
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      real(sp) :: ap(*)
      real(sp) :: x(*)
      integer(ik) :: incx
    end subroutine stpmv
  end interface
  interface
    subroutine stpsv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, sp
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      real(sp) :: ap(*)
      real(sp) :: x(*)
      integer(ik) :: incx
    end subroutine stpsv
  end interface
  interface
    subroutine strmv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      real(sp) :: a(lda, *)
      real(sp) :: x(*)
      integer(ik) :: incx
    end subroutine strmv
  end interface
  interface
    subroutine strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, sp
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: side
      character(len=1) :: uplo
      character(len=1) :: transa
      character(len=1) :: diag
      integer(ik) :: m
      integer(ik) :: n
      real(sp) :: alpha
      real(sp) :: a(lda, *)
      real(sp) :: b(ldb, *)
    end subroutine strsm
  end interface
  interface
    subroutine strsv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      real(sp) :: a(lda, *)
      real(sp) :: x(*)
      integer(ik) :: incx
    end subroutine strsv
  end interface

end module blas_s2
