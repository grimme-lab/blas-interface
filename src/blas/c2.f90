module blas_c2
  use blas_kinds, only : ik, sp
  implicit none
  private :: ik, sp

  interface
    subroutine cgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: trans
      integer(ik) :: m
      integer(ik) :: n
      integer(ik) :: kl
      integer(ik) :: ku
      complex(sp) :: alpha
      complex(sp) :: a(lda, *)
      complex(sp) :: x(*)
      integer(ik) :: incx
      complex(sp) :: beta
      complex(sp) :: y(*)
      integer(ik) :: incy
    end subroutine cgbmv
  end interface
  interface
    subroutine cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: trans
      integer(ik) :: m
      integer(ik) :: n
      complex(sp) :: alpha
      complex(sp) :: a(lda, *)
      complex(sp) :: x(*)
      integer(ik) :: incx
      complex(sp) :: beta
      complex(sp) :: y(*)
      integer(ik) :: incy
    end subroutine cgemv
  end interface
  interface
    subroutine cgerc(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer(ik) :: lda
      integer(ik) :: m
      integer(ik) :: n
      complex(sp) :: alpha
      complex(sp) :: x(*)
      integer(ik) :: incx
      complex(sp) :: y(*)
      integer(ik) :: incy
      complex(sp) :: a(lda, *)
    end subroutine cgerc
  end interface
  interface
    subroutine cgeru(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer(ik) :: lda
      integer(ik) :: m
      integer(ik) :: n
      complex(sp) :: alpha
      complex(sp) :: x(*)
      integer(ik) :: incx
      complex(sp) :: y(*)
      integer(ik) :: incy
      complex(sp) :: a(lda, *)
    end subroutine cgeru
  end interface
  interface
    subroutine chbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      integer(ik) :: k
      complex(sp) :: alpha
      complex(sp) :: a(lda, *)
      complex(sp) :: x(*)
      integer(ik) :: incx
      complex(sp) :: beta
      complex(sp) :: y(*)
      integer(ik) :: incy
    end subroutine chbmv
  end interface
  interface
    subroutine chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      complex(sp) :: alpha
      complex(sp) :: a(lda, *)
      complex(sp) :: x(*)
      integer(ik) :: incx
      complex(sp) :: beta
      complex(sp) :: y(*)
      integer(ik) :: incy
    end subroutine chemv
  end interface
  interface
    subroutine cher2(uplo, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      complex(sp) :: alpha
      complex(sp) :: x(*)
      integer(ik) :: incx
      complex(sp) :: y(*)
      integer(ik) :: incy
      complex(sp) :: a(lda, *)
    end subroutine cher2
  end interface
  interface
    subroutine chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
      import :: ik, sp
      character(len=1) :: uplo
      integer(ik) :: n
      complex(sp) :: alpha
      complex(sp) :: ap(*)
      complex(sp) :: x(*)
      integer(ik) :: incx
      complex(sp) :: beta
      complex(sp) :: y(*)
      integer(ik) :: incy
    end subroutine chpmv
  end interface
  interface
    subroutine chpr2(uplo, n, alpha, x, incx, y, incy, ap)
      import :: ik, sp
      character(len=1) :: uplo
      integer(ik) :: n
      complex(sp) :: alpha
      complex(sp) :: x(*)
      integer(ik) :: incx
      complex(sp) :: y(*)
      integer(ik) :: incy
      complex(sp) :: ap(*)
    end subroutine chpr2
  end interface
  interface
    subroutine chpr(uplo, n, alpha, x, incx, ap)
      import :: ik, sp
      character(len=1) :: uplo
      integer(ik) :: n
      real(sp) :: alpha
      complex(sp) :: x(*)
      integer(ik) :: incx
      complex(sp) :: ap(*)
    end subroutine chpr
  end interface
  interface
    subroutine cher(uplo, n, alpha, x, incx, a, lda)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      integer(ik) :: n
      real(sp) :: alpha
      complex(sp) :: x(*)
      integer(ik) :: incx
      complex(sp) :: a(lda, *)
    end subroutine cher
  end interface
  interface
    subroutine ctbmv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      integer(ik) :: k
      complex(sp) :: a(lda, *)
      complex(sp) :: x(*)
      integer(ik) :: incx
    end subroutine ctbmv
  end interface
  interface
    subroutine ctbsv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      integer(ik) :: k
      complex(sp) :: a(lda, *)
      complex(sp) :: x(*)
      integer(ik) :: incx
    end subroutine ctbsv
  end interface
  interface
    subroutine ctpmv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, sp
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      complex(sp) :: ap(*)
      complex(sp) :: x(*)
      integer(ik) :: incx
    end subroutine ctpmv
  end interface
  interface
    subroutine ctpsv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, sp
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      complex(sp) :: ap(*)
      complex(sp) :: x(*)
      integer(ik) :: incx
    end subroutine ctpsv
  end interface
  interface
    subroutine ctrmv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      complex(sp) :: a(lda, *)
      complex(sp) :: x(*)
      integer(ik) :: incx
    end subroutine ctrmv
  end interface
  interface
    subroutine ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, sp
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: side
      character(len=1) :: uplo
      character(len=1) :: transa
      character(len=1) :: diag
      integer(ik) :: m
      integer(ik) :: n
      complex(sp) :: alpha
      complex(sp) :: a(lda, *)
      complex(sp) :: b(ldb, *)
    end subroutine ctrsm
  end interface
  interface
    subroutine ctrsv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, sp
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      character(len=1) :: diag
      integer(ik) :: n
      complex(sp) :: a(lda, *)
      complex(sp) :: x(*)
      integer(ik) :: incx
    end subroutine ctrsv
  end interface

end module blas_c2
