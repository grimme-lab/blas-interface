module blas_s3
  use blas_kinds, only : ik, sp
  implicit none
  private :: ik, sp

  interface
    subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: transa
      character(len=1) :: transb
      integer(ik) :: m
      integer(ik) :: n
      integer(ik) :: k
      real(sp) :: alpha
      real(sp) :: a(lda, *)
      real(sp) :: b(ldb, *)
      real(sp) :: beta
      real(sp) :: c(ldc, *)
    end subroutine sgemm
  end interface
  interface
    subroutine ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: side
      character(len=1) :: uplo
      integer(ik) :: m
      integer(ik) :: n
      real(sp) :: alpha
      real(sp) :: a(lda, *)
      real(sp) :: b(ldb, *)
      real(sp) :: beta
      real(sp) :: c(ldc, *)
    end subroutine ssymm
  end interface
  interface
    subroutine ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik) :: ldc
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      integer(ik) :: n
      integer(ik) :: k
      real(sp) :: alpha
      real(sp) :: a(lda, *)
      real(sp) :: b(ldb, *)
      real(sp) :: beta
      real(sp) :: c(ldc, *)
    end subroutine ssyr2k
  end interface
  interface
    subroutine ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, sp
      integer(ik) :: ldc
      integer(ik) :: lda
      character(len=1) :: uplo
      character(len=1) :: trans
      integer(ik) :: n
      integer(ik) :: k
      real(sp) :: alpha
      real(sp) :: a(lda, *)
      real(sp) :: beta
      real(sp) :: c(ldc, *)
    end subroutine ssyrk
  end interface
  interface
    subroutine strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, sp
      integer(ik) :: ldb
      integer(ik) :: lda
      character(len=1) :: side
      character(len=1) :: uplo
      character(len=1) :: transa
      character(len=1) :: diag
      integer(ik) :: m
      integer(ik) :: n
      real(sp) :: alpha
      real(sp) :: a(lda, *)
      real(sp) :: b(ldb, *)
    end subroutine strmm
  end interface

end module blas_s3
