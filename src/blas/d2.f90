module blas_d2
  use blas_kinds, only : ik, dp
  implicit none
  private :: ik, dp

  interface
    subroutine dgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: trans
      integer(ik) :: m
      integer(ik) :: n
      integer(ik) :: kl
      integer(ik) :: ku
      real(dp) :: alpha
      real(dp) :: a(lda, *)
      real(dp) :: x(*)
      integer(ik) :: incx
      real(dp) :: beta
      real(dp) :: y(*)
      integer(ik) :: incy
    end subroutine dgbmv
  end interface
  interface
    subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: trans
      integer(ik) :: m
      integer(ik) :: n
      real(dp) :: alpha
      real(dp) :: a(lda, *)
      real(dp) :: x(*)
      integer(ik) :: incx
      real(dp) :: beta
      real(dp) :: y(*)
      integer(ik) :: incy
    end subroutine dgemv
  end interface
  interface
    subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer(ik) :: lda
      integer(ik) :: m
      integer(ik) :: n
      real(dp) :: alpha
      real(dp) :: x(*)
      integer(ik) :: incx
      real(dp) :: y(*)
      integer(ik) :: incy
      real(dp) :: a(lda, *)
    end subroutine dger
  end interface
  interface
    subroutine dsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      integer(ik) :: k
      real(dp) :: alpha
      real(dp) :: a(lda, *)
      real(dp) :: x(*)
      integer(ik) :: incx
      real(dp) :: beta
      real(dp) :: y(*)
      integer(ik) :: incy
    end subroutine dsbmv
  end interface
  interface
    subroutine dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
      import :: ik, dp
      character(len=1) :: uplo
      integer(ik) :: n
      real(dp) :: alpha
      real(dp) :: ap(*)
      real(dp) :: x(*)
      integer(ik) :: incx
      real(dp) :: beta
      real(dp) :: y(*)
      integer(ik) :: incy
    end subroutine dspmv
  end interface
  interface
    subroutine dspr2(uplo, n, alpha, x, incx, y, incy, ap)
      import :: ik, dp
      character(len=1) :: uplo
      integer(ik) :: n
      real(dp) :: alpha
      real(dp) :: x(*)
      integer(ik) :: incx
      real(dp) :: y(*)
      integer(ik) :: incy
      real(dp) :: ap(*)
    end subroutine dspr2
  end interface
  interface
    subroutine dspr(uplo, n, alpha, x, incx, ap)
      import :: ik, dp
      character(len=1) :: uplo
      integer(ik) :: n
      real(dp) :: alpha
      real(dp) :: x(*)
      integer(ik) :: incx
      real(dp) :: ap(*)
    end subroutine dspr
  end interface
  interface
    subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      real(dp) :: alpha
      real(dp) :: a(lda, *)
      real(dp) :: x(*)
      integer(ik) :: incx
      real(dp) :: beta
      real(dp) :: y(*)
      integer(ik) :: incy
    end subroutine dsymv
  end interface
  interface
    subroutine dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      real(dp) :: alpha
      real(dp) :: x(*)
      integer(ik) :: incx
      real(dp) :: y(*)
      integer(ik) :: incy
      real(dp) :: a(lda, *)
    end subroutine dsyr2
  end interface
  interface
    subroutine dsyr(uplo, n, alpha, x, incx, a, lda)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      real(dp) :: alpha
      real(dp) :: x(*)
      integer(ik) :: incx
      real(dp) :: a(lda, *)
    end subroutine dsyr
  end interface
  interface
    subroutine dtbmv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      integer(ik) :: k
      real(dp) :: a(lda, *)
      real(dp) :: x(*)
      integer(ik) :: incx
    end subroutine dtbmv
  end interface
  interface
    subroutine dtbsv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      integer(ik) :: k
      real(dp) :: a(lda, *)
      real(dp) :: x(*)
      integer(ik) :: incx
    end subroutine dtbsv
  end interface
  interface
    subroutine dtpmv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, dp
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      real(dp) :: ap(*)
      real(dp) :: x(*)
      integer(ik) :: incx
    end subroutine dtpmv
  end interface
  interface
    subroutine dtpsv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, dp
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      real(dp) :: ap(*)
      real(dp) :: x(*)
      integer(ik) :: incx
    end subroutine dtpsv
  end interface
  interface
    subroutine dtrmv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      real(dp) :: a(lda, *)
      real(dp) :: x(*)
      integer(ik) :: incx
    end subroutine dtrmv
  end interface
  interface
    subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, dp
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: side
      character(len=1) :: uplo
      character(len=1) :: transa
      character(len=1) :: diag
      integer(ik) :: m
      integer(ik) :: n
      real(dp) :: alpha
      real(dp) :: a(lda, *)
      real(dp) :: b(ldb, *)
    end subroutine dtrsm
  end interface
  interface
    subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      real(dp) :: a(lda, *)
      real(dp) :: x(*)
      integer(ik) :: incx
    end subroutine dtrsv
  end interface
end module blas_d2
