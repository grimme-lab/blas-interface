module blas_pure_d2
  use blas_kinds, only : ik, dp
  implicit none
  private :: ik, dp

  interface
    pure subroutine dgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: kl
      integer(ik), intent(in) :: ku
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(in) :: beta
      real(dp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine dgbmv
  end interface
  interface
    pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(in) :: beta
      real(dp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine dgemv
  end interface
  interface
    pure subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(dp), intent(inout) :: a(lda, *)
    end subroutine dger
  end interface
  interface
    pure subroutine dsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(in) :: beta
      real(dp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine dsbmv
  end interface
  interface
    pure subroutine dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
      import :: ik, dp
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: ap(*)
      real(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(in) :: beta
      real(dp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine dspmv
  end interface
  interface
    pure subroutine dspr2(uplo, n, alpha, x, incx, y, incy, ap)
      import :: ik, dp
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(dp), intent(inout) :: ap(*)
    end subroutine dspr2
  end interface
  interface
    pure subroutine dspr(uplo, n, alpha, x, incx, ap)
      import :: ik, dp
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(inout) :: ap(*)
    end subroutine dspr
  end interface
  interface
    pure subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(in) :: beta
      real(dp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine dsymv
  end interface
  interface
    pure subroutine dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(dp), intent(inout) :: a(lda, *)
    end subroutine dsyr2
  end interface
  interface
    pure subroutine dsyr(uplo, n, alpha, x, incx, a, lda)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(dp), intent(inout) :: a(lda, *)
    end subroutine dsyr
  end interface
  interface
    pure subroutine dtbmv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine dtbmv
  end interface
  interface
    pure subroutine dtbsv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine dtbsv
  end interface
  interface
    pure subroutine dtpmv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, dp
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: ap(*)
      real(dp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine dtpmv
  end interface
  interface
    pure subroutine dtpsv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, dp
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: ap(*)
      real(dp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine dtpsv
  end interface
  interface
    pure subroutine dtrmv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine dtrmv
  end interface
  interface
    pure subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, dp
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: side
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: transa
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(inout) :: b(ldb, *)
    end subroutine dtrsm
  end interface
  interface
    pure subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine dtrsv
  end interface
end module blas_pure_d2
