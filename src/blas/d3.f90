module blas_d3
  use blas_kinds, only : ik, dp
  implicit none

  interface
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: transa
      character(len=1) :: transb
      integer(ik) :: m
      integer(ik) :: n
      integer(ik) :: k
      real(dp) :: alpha
      real(dp) :: a(lda, *)
      real(dp) :: b(ldb, *)
      real(dp) :: beta
      real(dp) :: c(ldc, *)
    end subroutine dgemm
  end interface
  interface
    subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: side
      character(len=1) :: uplo
      integer(ik) :: m
      integer(ik) :: n
      real(dp) :: alpha
      real(dp) :: a(lda, *)
      real(dp) :: b(ldb, *)
      real(dp) :: beta
      real(dp) :: c(ldc, *)
    end subroutine dsymm
  end interface
  interface
    subroutine dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      integer(ik) :: n
      integer(ik) :: k
      real(dp) :: alpha
      real(dp) :: a(lda, *)
      real(dp) :: b(ldb, *)
      real(dp) :: beta
      real(dp) :: c(ldc, *)
    end subroutine dsyr2k
  end interface
  interface
    subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, dp
      integer(ik) :: ldc
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      integer(ik) :: n
      integer(ik) :: k
      real(dp) :: alpha
      real(dp) :: a(lda, *)
      real(dp) :: beta
      real(dp) :: c(ldc, *)
    end subroutine dsyrk
  end interface
  interface
    subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
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
    end subroutine dtrmm
  end interface
end module blas_d3
