module blas_pure_c2
  use blas_kinds, only : ik, sp
  implicit none
  private :: ik, sp

  interface
    pure subroutine cgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: kl
      integer(ik), intent(in) :: ku
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(in) :: beta
      complex(sp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine cgbmv
  end interface
  interface
    pure subroutine cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(in) :: beta
      complex(sp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine cgemv
  end interface
  interface
    pure subroutine cgerc(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(sp), intent(inout) :: a(lda, *)
    end subroutine cgerc
  end interface
  interface
    pure subroutine cgeru(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(sp), intent(inout) :: a(lda, *)
    end subroutine cgeru
  end interface
  interface
    pure subroutine chbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(in) :: beta
      complex(sp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine chbmv
  end interface
  interface
    pure subroutine chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(in) :: beta
      complex(sp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine chemv
  end interface
  interface
    pure subroutine cher2(uplo, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(sp), intent(inout) :: a(lda, *)
    end subroutine cher2
  end interface
  interface
    pure subroutine cher(uplo, n, alpha, x, incx, a, lda)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: alpha
      complex(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(inout) :: a(lda, *)
    end subroutine cher
  end interface
  interface
    pure subroutine chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
      import :: ik, sp
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: ap(*)
      complex(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(in) :: beta
      complex(sp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine chpmv
  end interface
  interface
    pure subroutine chpr2(uplo, n, alpha, x, incx, y, incy, ap)
      import :: ik, sp
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(sp), intent(inout) :: ap(*)
    end subroutine chpr2
  end interface
  interface
    pure subroutine chpr(uplo, n, alpha, x, incx, ap)
      import :: ik, sp
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: alpha
      complex(sp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(sp), intent(inout) :: ap(*)
    end subroutine chpr
  end interface
  interface
    pure subroutine ctbmv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ctbmv
  end interface
  interface
    pure subroutine ctbsv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ctbsv
  end interface
  interface
    pure subroutine ctpmv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, sp
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: ap(*)
      complex(sp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ctpmv
  end interface
  interface
    pure subroutine ctpsv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, sp
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: ap(*)
      complex(sp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ctpsv
  end interface
  interface
    pure subroutine ctrmv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ctrmv
  end interface
  interface
    pure subroutine ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, sp
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: side
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: transa
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(inout) :: b(ldb, *)
    end subroutine ctrsm
  end interface
  interface
    pure subroutine ctrsv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, sp
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ctrsv
  end interface

end module blas_pure_c2
