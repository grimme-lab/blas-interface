module blas_z2
  use blas_kinds, only : ik, dp
  implicit none
  private :: ik, dp

  interface
    subroutine zgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: trans
      integer(ik) :: m
      integer(ik) :: n
      integer(ik) :: kl
      integer(ik) :: ku
      complex(dp) :: alpha
      complex(dp) :: a(lda, *)
      complex(dp) :: x(*)
      integer(ik) :: incx
      complex(dp) :: beta
      complex(dp) :: y(*)
      integer(ik) :: incy
    end subroutine zgbmv
  end interface
  interface
    subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: trans
      integer(ik) :: m
      integer(ik) :: n
      complex(dp) :: alpha
      complex(dp) :: a(lda, *)
      complex(dp) :: x(*)
      integer(ik) :: incx
      complex(dp) :: beta
      complex(dp) :: y(*)
      integer(ik) :: incy
    end subroutine zgemv
  end interface
  interface
    subroutine zgerc(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer(ik) :: lda
      integer(ik) :: m
      integer(ik) :: n
      complex(dp) :: alpha
      complex(dp) :: x(*)
      integer(ik) :: incx
      complex(dp) :: y(*)
      integer(ik) :: incy
      complex(dp) :: a(lda, *)
    end subroutine zgerc
  end interface
  interface
    subroutine zgeru(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer(ik) :: lda
      integer(ik) :: m
      integer(ik) :: n
      complex(dp) :: alpha
      complex(dp) :: x(*)
      integer(ik) :: incx
      complex(dp) :: y(*)
      integer(ik) :: incy
      complex(dp) :: a(lda, *)
    end subroutine zgeru
  end interface
  interface
    subroutine zhbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      integer(ik) :: k
      complex(dp) :: alpha
      complex(dp) :: a(lda, *)
      complex(dp) :: x(*)
      integer(ik) :: incx
      complex(dp) :: beta
      complex(dp) :: y(*)
      integer(ik) :: incy
    end subroutine zhbmv
  end interface
  interface
    subroutine zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      complex(dp) :: alpha
      complex(dp) :: a(lda, *)
      complex(dp) :: x(*)
      integer(ik) :: incx
      complex(dp) :: beta
      complex(dp) :: y(*)
      integer(ik) :: incy
    end subroutine zhemv
  end interface
  interface
    subroutine zher2(uplo, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      complex(dp) :: alpha
      complex(dp) :: x(*)
      integer(ik) :: incx
      complex(dp) :: y(*)
      integer(ik) :: incy
      complex(dp) :: a(lda, *)
    end subroutine zher2
  end interface
  interface
    subroutine zher(uplo, n, alpha, x, incx, a, lda)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      real(dp) :: alpha
      complex(dp) :: x(*)
      integer(ik) :: incx
      complex(dp) :: a(lda, *)
    end subroutine zher
  end interface
  interface
    subroutine zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
      import :: ik, dp
      character(len=1) :: uplo
      integer(ik) :: n
      complex(dp) :: alpha
      complex(dp) :: ap(*)
      complex(dp) :: x(*)
      integer(ik) :: incx
      complex(dp) :: beta
      complex(dp) :: y(*)
      integer(ik) :: incy
    end subroutine zhpmv
  end interface
  interface
    subroutine zhpr(uplo, n, alpha, x, incx, ap)
      import :: ik, dp
      character(len=1) :: uplo
      integer(ik) :: n
      real(dp) :: alpha
      complex(dp) :: x(*)
      integer(ik) :: incx
      complex(dp) :: ap(*)
    end subroutine zhpr
  end interface
  interface
    subroutine zhpr2(uplo, n, alpha, x, incx, y, incy, ap)
      import :: ik, dp
      character(len=1) :: uplo
      integer(ik) :: n
      complex(dp) :: alpha
      complex(dp) :: x(*)
      integer(ik) :: incx
      complex(dp) :: y(*)
      integer(ik) :: incy
      complex(dp) :: ap(*)
    end subroutine zhpr2
  end interface
  interface
    subroutine ztbmv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      integer(ik) :: k
      complex(dp) :: a(lda, *)
      complex(dp) :: x(*)
      integer(ik) :: incx
    end subroutine ztbmv
  end interface
  interface
    subroutine ztbsv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      integer(ik) :: k
      complex(dp) :: a(lda, *)
      complex(dp) :: x(*)
      integer(ik) :: incx
    end subroutine ztbsv
  end interface
  interface
    subroutine ztpmv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, dp
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      complex(dp) :: ap(*)
      complex(dp) :: x(*)
      integer(ik) :: incx
    end subroutine ztpmv
  end interface
  interface
    subroutine ztpsv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, dp
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      complex(dp) :: ap(*)
      complex(dp) :: x(*)
      integer(ik) :: incx
    end subroutine ztpsv
  end interface
  interface
    subroutine ztrmv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      complex(dp) :: a(lda, *)
      complex(dp) :: x(*)
      integer(ik) :: incx
    end subroutine ztrmv
  end interface
  interface
    subroutine ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, dp
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: side
      character(len=1) :: uplo
      character(len=1) :: transa
      character(len=1) :: diag
      integer(ik) :: m
      integer(ik) :: n
      complex(dp) :: alpha
      complex(dp) :: a(lda, *)
      complex(dp) :: b(ldb, *)
    end subroutine ztrsm
  end interface
  interface
    subroutine ztrsv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, dp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      complex(dp) :: a(lda, *)
      complex(dp) :: x(*)
      integer(ik) :: incx
    end subroutine ztrsv
  end interface
end module blas_z2
