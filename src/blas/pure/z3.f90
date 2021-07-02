module blas_pure_z3
  use blas_kinds, only : ik, dp
  implicit none
  private :: ik, dp

  interface
    pure subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: transa
      character(len=1), intent(in) :: transb
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(in) :: b(ldb, *)
      complex(dp), intent(in) :: beta
      complex(dp), intent(inout) :: c(ldc, *)
    end subroutine zgemm
  end interface
  interface
    pure subroutine zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: side
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(in) :: b(ldb, *)
      complex(dp), intent(in) :: beta
      complex(dp), intent(inout) :: c(ldc, *)
    end subroutine zhemm
  end interface
  interface
    pure subroutine zher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(in) :: b(ldb, *)
      real(dp), intent(in) :: beta
      complex(dp), intent(inout) :: c(ldc, *)
    end subroutine zher2k
  end interface
  interface
    pure subroutine zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(dp), intent(in) :: alpha
      complex(dp), intent(in) :: a(lda, *)
      real(dp), intent(in) :: beta
      complex(dp), intent(inout) :: c(ldc, *)
    end subroutine zherk
  end interface
  interface
    pure subroutine zsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: side
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(in) :: b(ldb, *)
      complex(dp), intent(in) :: beta
      complex(dp), intent(inout) :: c(ldc, *)
    end subroutine zsymm
  end interface
  interface
    pure subroutine zsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(in) :: b(ldb, *)
      complex(dp), intent(in) :: beta
      complex(dp), intent(inout) :: c(ldc, *)
    end subroutine zsyr2k
  end interface
  interface
    pure subroutine zsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(dp), intent(in) :: alpha
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(in) :: beta
      complex(dp), intent(inout) :: c(ldc, *)
    end subroutine zsyrk
  end interface
  interface
    pure subroutine ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
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
    end subroutine ztrmm
  end interface
end module blas_pure_z3
