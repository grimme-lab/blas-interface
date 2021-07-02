module blas_z3
  use blas_kinds, only : ik, dp
  implicit none
  private :: ik, dp

  interface
    subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: transa
      character(len=1) :: transb
      integer(ik) :: m
      integer(ik) :: n
      integer(ik) :: k
      complex(dp) :: alpha
      complex(dp) :: a(lda, *)
      complex(dp) :: b(ldb, *)
      complex(dp) :: beta
      complex(dp) :: c(ldc, *)
    end subroutine zgemm
  end interface
  interface
    subroutine zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: side
      character(len=1) :: uplo
      integer(ik) :: m
      integer(ik) :: n
      complex(dp) :: alpha
      complex(dp) :: a(lda, *)
      complex(dp) :: b(ldb, *)
      complex(dp) :: beta
      complex(dp) :: c(ldc, *)
    end subroutine zhemm
  end interface
  interface
    subroutine zher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      integer(ik) :: n
      integer(ik) :: k
      complex(dp) :: alpha
      complex(dp) :: a(lda, *)
      complex(dp) :: b(ldb, *)
      real(dp) :: beta
      complex(dp) :: c(ldc, *)
    end subroutine zher2k
  end interface
  interface
    subroutine zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, dp
      integer(ik) :: ldc
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      integer(ik) :: n
      integer(ik) :: k
      real(dp) :: alpha
      complex(dp) :: a(lda, *)
      real(dp) :: beta
      complex(dp) :: c(ldc, *)
    end subroutine zherk
  end interface
  interface
    subroutine zsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: side
      character(len=1) :: uplo
      integer(ik) :: m
      integer(ik) :: n
      complex(dp) :: alpha
      complex(dp) :: a(lda, *)
      complex(dp) :: b(ldb, *)
      complex(dp) :: beta
      complex(dp) :: c(ldc, *)
    end subroutine zsymm
  end interface
  interface
    subroutine zsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      integer(ik) :: n
      integer(ik) :: k
      complex(dp) :: alpha
      complex(dp) :: a(lda, *)
      complex(dp) :: b(ldb, *)
      complex(dp) :: beta
      complex(dp) :: c(ldc, *)
    end subroutine zsyr2k
  end interface
  interface
    subroutine zsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, dp
      integer(ik) :: ldc
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      integer(ik) :: n
      integer(ik) :: k
      complex(dp) :: alpha
      complex(dp) :: a(lda, *)
      complex(dp) :: beta
      complex(dp) :: c(ldc, *)
    end subroutine zsyrk
  end interface
  interface
    subroutine ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
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
    end subroutine ztrmm
  end interface
end module blas_z3
