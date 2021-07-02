module blas_c3
  use blas_kinds, only : ik, sp
  implicit none
  private :: ik, sp

  interface
    subroutine cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: transa
      character(len=1) :: transb
      integer(ik) :: m
      integer(ik) :: n
      integer(ik) :: k
      complex(sp) :: alpha
      complex(sp) :: a(lda, *)
      complex(sp) :: b(ldb, *)
      complex(sp) :: beta
      complex(sp) :: c(ldc, *)
    end subroutine cgemm
  end interface
  interface
    subroutine chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: side
      character(len=1) :: uplo
      integer(ik) :: m
      integer(ik) :: n
      complex(sp) :: alpha
      complex(sp) :: a(lda, *)
      complex(sp) :: b(ldb, *)
      complex(sp) :: beta
      complex(sp) :: c(ldc, *)
    end subroutine chemm
  end interface
  interface
    subroutine cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      integer(ik) :: n
      integer(ik) :: k
      complex(sp) :: alpha
      complex(sp) :: a(lda, *)
      complex(sp) :: b(ldb, *)
      real(sp) :: beta
      complex(sp) :: c(ldc, *)
    end subroutine cher2k
  end interface
  interface
    subroutine cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, sp
      integer(ik) :: ldc
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      integer(ik) :: n
      integer(ik) :: k
      real(sp) :: alpha
      complex(sp) :: a(lda, *)
      real(sp) :: beta
      complex(sp) :: c(ldc, *)
    end subroutine cherk
  end interface
  interface
    subroutine csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: side
      character(len=1) :: uplo
      integer(ik) :: m
      integer(ik) :: n
      complex(sp) :: alpha
      complex(sp) :: a(lda, *)
      complex(sp) :: b(ldb, *)
      complex(sp) :: beta
      complex(sp) :: c(ldc, *)
    end subroutine csymm
  end interface
  interface
    subroutine csyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      integer(ik) :: n
      integer(ik) :: k
      complex(sp) :: alpha
      complex(sp) :: a(lda, *)
      complex(sp) :: b(ldb, *)
      complex(sp) :: beta
      complex(sp) :: c(ldc, *)
    end subroutine csyr2k
  end interface
  interface
    subroutine csyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, sp
      integer(ik) :: ldc
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      integer(ik) :: n
      integer(ik) :: k
      complex(sp) :: alpha
      complex(sp) :: a(lda, *)
      complex(sp) :: beta
      complex(sp) :: c(ldc, *)
    end subroutine csyrk
  end interface
  interface
    subroutine ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
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
    end subroutine ctrmm
  end interface

end module blas_c3
