module blas_pure_s3
  use blas_kinds, only : ik, sp
  implicit none
  private :: ik, sp

  interface
    pure subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: transa
      character(len=1), intent(in) :: transb
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(in) :: b(ldb, *)
      real(sp), intent(in) :: beta
      real(sp), intent(inout) :: c(ldc, *)
    end subroutine sgemm
  end interface
  interface
    pure subroutine ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: side
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(in) :: b(ldb, *)
      real(sp), intent(in) :: beta
      real(sp), intent(inout) :: c(ldc, *)
    end subroutine ssymm
  end interface
  interface
    pure subroutine ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(in) :: b(ldb, *)
      real(sp), intent(in) :: beta
      real(sp), intent(inout) :: c(ldc, *)
    end subroutine ssyr2k
  end interface
  interface
    pure subroutine ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(in) :: beta
      real(sp), intent(inout) :: c(ldc, *)
    end subroutine ssyrk
  end interface
  interface
    pure subroutine strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, sp
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: side
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: transa
      character(len=1), intent(in) :: diag
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: a(lda, *)
      real(sp), intent(inout) :: b(ldb, *)
    end subroutine strmm
  end interface

end module blas_pure_s3
