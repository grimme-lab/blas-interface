module blas_pure_z2
  use blas_kinds, only : ik, dp
  implicit none
  private :: ik, dp

  interface
    pure subroutine zgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: kl
      integer(ik), intent(in) :: ku
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(in) :: beta
      complex(dp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine zgbmv
  end interface
  interface
    pure subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(in) :: beta
      complex(dp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine zgemv
  end interface
  interface
    pure subroutine zgerc(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(dp), intent(inout) :: a(lda, *)
    end subroutine zgerc
  end interface
  interface
    pure subroutine zgeru(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(dp), intent(inout) :: a(lda, *)
    end subroutine zgeru
  end interface
  interface
    pure subroutine zhbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(in) :: beta
      complex(dp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine zhbmv
  end interface
  interface
    pure subroutine zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(in) :: beta
      complex(dp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine zhemv
  end interface
  interface
    pure subroutine zher2(uplo, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(dp), intent(inout) :: a(lda, *)
    end subroutine zher2
  end interface
  interface
    pure subroutine zher(uplo, n, alpha, x, incx, a, lda)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: alpha
      complex(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(inout) :: a(lda, *)
    end subroutine zher
  end interface
  interface
    pure subroutine zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
      import :: ik, dp
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: ap(*)
      complex(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(in) :: beta
      complex(dp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine zhpmv
  end interface
  interface
    pure subroutine zhpr(uplo, n, alpha, x, incx, ap)
      import :: ik, dp
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: alpha
      complex(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(inout) :: ap(*)
    end subroutine zhpr
  end interface
  interface
    pure subroutine zhpr2(uplo, n, alpha, x, incx, y, incy, ap)
      import :: ik, dp
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(dp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(dp), intent(inout) :: ap(*)
    end subroutine zhpr2
  end interface
  interface
    pure subroutine ztbmv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ztbmv
  end interface
  interface
    pure subroutine ztbsv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ztbsv
  end interface
  interface
    pure subroutine ztpmv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, dp
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: ap(*)
      complex(dp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ztpmv
  end interface
  interface
    pure subroutine ztpsv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, dp
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: ap(*)
      complex(dp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ztpsv
  end interface
  interface
    pure subroutine ztrmv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ztrmv
  end interface
  interface
    pure subroutine ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, dp
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: side
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: transa
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(inout) :: b(ldb, *)
    end subroutine ztrsm
  end interface
  interface
    pure subroutine ztrsv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, dp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ztrsv
  end interface
end module blas_pure_z2
