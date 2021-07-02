module blas_pure_d3
  use blas_kinds, only : ik, dp
  implicit none

  interface
    pure subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: transa
      character(len=1), intent(in) :: transb
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(in) :: b(ldb, *)
      real(dp), intent(in) :: beta
      real(dp), intent(inout) :: c(ldc, *)
    end subroutine dgemm
  end interface
  interface
    pure subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: side
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(in) :: b(ldb, *)
      real(dp), intent(in) :: beta
      real(dp), intent(inout) :: c(ldc, *)
    end subroutine dsymm
  end interface
  interface
    pure subroutine dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(in) :: b(ldb, *)
      real(dp), intent(in) :: beta
      real(dp), intent(inout) :: c(ldc, *)
    end subroutine dsyr2k
  end interface
  interface
    pure subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(in) :: beta
      real(dp), intent(inout) :: c(ldc, *)
    end subroutine dsyrk
  end interface
  interface
    pure subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
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
    end subroutine dtrmm
  end interface
end module blas_pure_d3
