module blas_pure_c3
  use blas_kinds, only : ik, sp
  implicit none
  private :: ik, sp

  interface
    pure subroutine cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: transa
      character(len=1), intent(in) :: transb
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(in) :: b(ldb, *)
      complex(sp), intent(in) :: beta
      complex(sp), intent(inout) :: c(ldc, *)
    end subroutine cgemm
  end interface
  interface
    pure subroutine chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: side
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(in) :: b(ldb, *)
      complex(sp), intent(in) :: beta
      complex(sp), intent(inout) :: c(ldc, *)
    end subroutine chemm
  end interface
  interface
    pure subroutine cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(in) :: b(ldb, *)
      real(sp), intent(in) :: beta
      complex(sp), intent(inout) :: c(ldc, *)
    end subroutine cher2k
  end interface
  interface
    pure subroutine cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(sp), intent(in) :: alpha
      complex(sp), intent(in) :: a(lda, *)
      real(sp), intent(in) :: beta
      complex(sp), intent(inout) :: c(ldc, *)
    end subroutine cherk
  end interface
  interface
    pure subroutine csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: side
      character(len=1), intent(in) :: uplo
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(in) :: b(ldb, *)
      complex(sp), intent(in) :: beta
      complex(sp), intent(inout) :: c(ldc, *)
    end subroutine csymm
  end interface
  interface
    pure subroutine csyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(in) :: b(ldb, *)
      complex(sp), intent(in) :: beta
      complex(sp), intent(inout) :: c(ldc, *)
    end subroutine csyr2k
  end interface
  interface
    pure subroutine csyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: lda
      character(len=1), intent(in) :: uplo
      character(len=1), intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(sp), intent(in) :: alpha
      complex(sp), intent(in) :: a(lda, *)
      complex(sp), intent(in) :: beta
      complex(sp), intent(inout) :: c(ldc, *)
    end subroutine csyrk
  end interface
  interface
    pure subroutine ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
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
    end subroutine ctrmm
  end interface

end module blas_pure_c3
