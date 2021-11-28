!> Interfaces to level 3 BLAS routines.
module blas_pure_level3
  use blas_kinds, only : ik, sp, dp
  implicit none
  private

  public :: blas_gemm, blas_hemm, blas_her2k, blas_herk, blas_symm, blas_syr2k, blas_syrk, &
    & blas_trsm, blas_trmm

  !> Performs one of the matrix-matrix operations
  !>
  !>    C := alpha*op( A )*op( B ) + beta*C,
  !>
  !> where  op( X ) is one of
  !>
  !>    op( X ) = X   or   op( X ) = X**T,
  !>
  !> alpha and beta are scalars, and A, B and C are matrices, with op( A )
  !> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  interface blas_gemm
    pure subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: transa
      character, intent(in) :: transb
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: b(ldb, *)
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: c(ldc, *)
    end subroutine sgemm
    pure subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: transa
      character, intent(in) :: transb
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: b(ldb, *)
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: c(ldc, *)
    end subroutine dgemm

    pure subroutine cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: transa
      character, intent(in) :: transb
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: b(ldb, *)
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine cgemm
    pure subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: transa
      character, intent(in) :: transb
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: b(ldb, *)
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine zgemm

    module procedure :: blas_gemm_rsp
    module procedure :: blas_gemm_csp
    module procedure :: blas_gemm_rdp
    module procedure :: blas_gemm_cdp
  end interface blas_gemm

  interface blas_hemm
    pure subroutine chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: b(ldb, *)
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine chemm
    pure subroutine zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: b(ldb, *)
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine zhemm

    module procedure :: blas_hemm_csp
    module procedure :: blas_hemm_cdp
  end interface blas_hemm

  interface blas_her2k
    pure subroutine cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: b(ldb, *)
      real(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine cher2k
    pure subroutine zher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: b(ldb, *)
      real(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine zher2k

    module procedure :: blas_her2k_csp
    module procedure :: blas_her2k_cdp
  end interface blas_her2k

  interface blas_herk
    pure subroutine cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine cherk
    pure subroutine zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine zherk

    module procedure :: blas_herk_csp
    module procedure :: blas_herk_cdp
  end interface blas_herk

  !> Performs one of the matrix-matrix operations
  !>
  !>    C := alpha*A*B + beta*C,
  !>
  !> or
  !>
  !>    C := alpha*B*A + beta*C,
  !>
  !> where alpha and beta are scalars,  A is a symmetric matrix and  B and
  !> C are  m by n matrices.
  interface blas_symm
    pure subroutine ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: b(ldb, *)
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: c(ldc, *)
    end subroutine ssymm
    pure subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: b(ldb, *)
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: c(ldc, *)
    end subroutine dsymm

    pure subroutine csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: b(ldb, *)
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine csymm
    pure subroutine zsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: b(ldb, *)
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine zsymm

    module procedure :: blas_symm_rsp
    module procedure :: blas_symm_csp
    module procedure :: blas_symm_rdp
    module procedure :: blas_symm_cdp
  end interface blas_symm

  !> Performs one of the symmetric rank 2k operations
  !>
  !>    C := alpha*A*B**T + alpha*B*A**T + beta*C,
  !>
  !> or
  !>
  !>    C := alpha*A**T*B + alpha*B**T*A + beta*C,
  !>
  !> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
  !> and  A and B  are  n by k  matrices  in the  first  case  and  k by n
  !> matrices in the second case.
  interface blas_syr2k
    pure subroutine ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: b(ldb, *)
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: c(ldc, *)
    end subroutine ssyr2k
    pure subroutine dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: b(ldb, *)
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: c(ldc, *)
    end subroutine dsyr2k

    pure subroutine csyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: b(ldb, *)
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine csyr2k
    pure subroutine zsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: b(ldb, *)
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine zsyr2k

    module procedure :: blas_syr2k_rsp
    module procedure :: blas_syr2k_csp
    module procedure :: blas_syr2k_rdp
    module procedure :: blas_syr2k_cdp
  end interface blas_syr2k

  !> Performs one of the symmetric rank k operations
  !>
  !>    C := alpha*A*A**T + beta*C,
  !>
  !> or
  !>
  !>    C := alpha*A**T*A + beta*C,
  !>
  !> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
  !> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
  !> in the second case.
  interface blas_syrk
    pure subroutine ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: c(ldc, *)
    end subroutine ssyrk
    pure subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: c(ldc, *)
    end subroutine dsyrk

    pure subroutine csyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine csyrk
    pure subroutine zsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldc
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: c(ldc, *)
    end subroutine zsyrk

    module procedure :: blas_syrk_rsp
    module procedure :: blas_syrk_csp
    module procedure :: blas_syrk_rdp
    module procedure :: blas_syrk_cdp
  end interface blas_syrk

  !> Performs one of the matrix-matrix operations
  !>
  !>    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
  !>
  !> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
  !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
  !>
  !>    op( A ) = A   or   op( A ) = A**T.
  interface blas_trmm
    pure subroutine strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      character, intent(in) :: transa
      character, intent(in) :: diag
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: b(ldb, *)
    end subroutine strmm
    pure subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      character, intent(in) :: transa
      character, intent(in) :: diag
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: b(ldb, *)
    end subroutine dtrmm

    pure subroutine ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      character, intent(in) :: transa
      character, intent(in) :: diag
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: b(ldb, *)
    end subroutine ctrmm
    pure subroutine ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      character, intent(in) :: transa
      character, intent(in) :: diag
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: b(ldb, *)
    end subroutine ztrmm

    module procedure :: blas_trmm_rsp
    module procedure :: blas_trmm_csp
    module procedure :: blas_trmm_rdp
    module procedure :: blas_trmm_cdp
  end interface blas_trmm

  !> Solves one of the matrix equations
  !>
  !>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
  !>
  !> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
  !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
  !>
  !>    op( A ) = A   or   op( A ) = A**T.
  !>
  !> The matrix X is overwritten on B.
  interface blas_trsm
    pure subroutine strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      character, intent(in) :: transa
      character, intent(in) :: diag
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: b(ldb, *)
    end subroutine strsm
    pure subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      character, intent(in) :: transa
      character, intent(in) :: diag
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: b(ldb, *)
    end subroutine dtrsm

    pure subroutine ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      character, intent(in) :: transa
      character, intent(in) :: diag
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: b(ldb, *)
    end subroutine ctrsm
    pure subroutine ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: lda
      character, intent(in) :: side
      character, intent(in) :: uplo
      character, intent(in) :: transa
      character, intent(in) :: diag
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: b(ldb, *)
    end subroutine ztrsm

    module procedure :: blas_trsm_rsp
    module procedure :: blas_trsm_csp
    module procedure :: blas_trsm_rdp
    module procedure :: blas_trsm_cdp
  end interface blas_trsm

contains

  pure subroutine blas_gemm_rsp(amat, bmat, cmat, transa, transb, alpha, beta)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: bmat(:, :)
    real(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: transa
    character, intent(in), optional :: transb
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: tra, trb
    real(wp) :: a, b
    integer(ik) :: m, n, k, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(transa)) tra = transa
    trb = 'n'
    if (present(transb)) trb = transb
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    m = size(cmat, 1)
    n = size(cmat, 2)
    call blas_gemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_gemm_rsp

  pure subroutine blas_gemm_csp(amat, bmat, cmat, transa, transb, alpha, beta)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: bmat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: transa
    character, intent(in), optional :: transb
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: tra, trb
    complex(wp) :: a, b
    integer(ik) :: m, n, k, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(transa)) tra = transa
    trb = 'n'
    if (present(transb)) trb = transb
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    m = size(cmat, 1)
    n = size(cmat, 2)
    call blas_gemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_gemm_csp

  pure subroutine blas_gemm_rdp(amat, bmat, cmat, transa, transb, alpha, beta)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: bmat(:, :)
    real(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: transa
    character, intent(in), optional :: transb
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: tra, trb
    real(wp) :: a, b
    integer(ik) :: m, n, k, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(transa)) tra = transa
    trb = 'n'
    if (present(transb)) trb = transb
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    m = size(cmat, 1)
    n = size(cmat, 2)
    call blas_gemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_gemm_rdp

  pure subroutine blas_gemm_cdp(amat, bmat, cmat, transa, transb, alpha, beta)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: bmat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: transa
    character, intent(in), optional :: transb
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: tra, trb
    complex(wp) :: a, b
    integer(ik) :: m, n, k, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(transa)) tra = transa
    trb = 'n'
    if (present(transb)) trb = transb
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    m = size(cmat, 1)
    n = size(cmat, 2)
    call blas_gemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_gemm_cdp


  pure subroutine blas_symm_rsp(amat, bmat, cmat, side, uplo, alpha, beta)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: bmat(:, :)
    real(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: lra, ula
    real(wp) :: a, b
    integer(ik) :: m, n, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    lra = 'l'
    if (present(side)) lra = side
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    m = size(cmat, 1)
    n = size(cmat, 2)
    call blas_symm(lra, ula, m, n, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_symm_rsp

  pure subroutine blas_symm_csp(amat, bmat, cmat, side, uplo, alpha, beta)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: bmat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: lra, ula
    complex(wp) :: a, b
    integer(ik) :: m, n, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    lra = 'l'
    if (present(side)) lra = side
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    m = size(cmat, 1)
    n = size(cmat, 2)
    call blas_symm(lra, ula, m, n, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_symm_csp

  pure subroutine blas_hemm_csp(amat, bmat, cmat, side, uplo, alpha, beta)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: bmat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: lra, ula
    complex(wp) :: a, b
    integer(ik) :: m, n, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    lra = 'l'
    if (present(side)) lra = side
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    m = size(cmat, 1)
    n = size(cmat, 2)
    call blas_hemm(lra, ula, m, n, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_hemm_csp

  pure subroutine blas_symm_rdp(amat, bmat, cmat, side, uplo, alpha, beta)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: bmat(:, :)
    real(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: lra, ula
    real(wp) :: a, b
    integer(ik) :: m, n, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    lra = 'l'
    if (present(side)) lra = side
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    m = size(cmat, 1)
    n = size(cmat, 2)
    call blas_symm(lra, ula, m, n, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_symm_rdp

  pure subroutine blas_symm_cdp(amat, bmat, cmat, side, uplo, alpha, beta)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: bmat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: lra, ula
    complex(wp) :: a, b
    integer(ik) :: m, n, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    lra = 'l'
    if (present(side)) lra = side
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    m = size(cmat, 1)
    n = size(cmat, 2)
    call blas_symm(lra, ula, m, n, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_symm_cdp

  pure subroutine blas_hemm_cdp(amat, bmat, cmat, side, uplo, alpha, beta)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: bmat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: lra, ula
    complex(wp) :: a, b
    integer(ik) :: m, n, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    lra = 'l'
    if (present(side)) lra = side
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    m = size(cmat, 1)
    n = size(cmat, 2)
    call blas_hemm(lra, ula, m, n, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_hemm_cdp


  pure subroutine blas_syrk_rsp(amat, cmat, uplo, trans, alpha, beta)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula, tra
    real(wp) :: a, b
    integer(ik) :: n, k, lda, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    ula = 'u'
    if (present(uplo)) ula = uplo
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldc = max(1, size(cmat, 1))
    n = size(cmat, 2)
    call blas_syrk(ula, tra, n, k, a, amat, lda, b, cmat, ldc)
  end subroutine blas_syrk_rsp

  pure subroutine blas_syrk_csp(amat, cmat, uplo, trans, alpha, beta)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: ula, tra
    complex(wp) :: a, b
    integer(ik) :: n, k, lda, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    ula = 'u'
    if (present(uplo)) ula = uplo
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldc = max(1, size(cmat, 1))
    n = size(cmat, 2)
    call blas_syrk(ula, tra, n, k, a, amat, lda, b, cmat, ldc)
  end subroutine blas_syrk_csp

  pure subroutine blas_herk_csp(amat, cmat, uplo, trans, alpha, beta)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula, tra
    real(wp) :: a, b
    integer(ik) :: n, k, lda, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    ula = 'u'
    if (present(uplo)) ula = uplo
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldc = max(1, size(cmat, 1))
    n = size(cmat, 2)
    call blas_herk(ula, tra, n, k, a, amat, lda, b, cmat, ldc)
  end subroutine blas_herk_csp

  pure subroutine blas_syrk_rdp(amat, cmat, uplo, trans, alpha, beta)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula, tra
    real(wp) :: a, b
    integer(ik) :: n, k, lda, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    ula = 'u'
    if (present(uplo)) ula = uplo
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldc = max(1, size(cmat, 1))
    n = size(cmat, 2)
    call blas_syrk(ula, tra, n, k, a, amat, lda, b, cmat, ldc)
  end subroutine blas_syrk_rdp

  pure subroutine blas_syrk_cdp(amat, cmat, uplo, trans, alpha, beta)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: ula, tra
    complex(wp) :: a, b
    integer(ik) :: n, k, lda, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    ula = 'u'
    if (present(uplo)) ula = uplo
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldc = max(1, size(cmat, 1))
    n = size(cmat, 2)
    call blas_syrk(ula, tra, n, k, a, amat, lda, b, cmat, ldc)
  end subroutine blas_syrk_cdp

  pure subroutine blas_herk_cdp(amat, cmat, uplo, trans, alpha, beta)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula, tra
    real(wp) :: a, b
    integer(ik) :: n, k, lda, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    ula = 'u'
    if (present(uplo)) ula = uplo
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldc = max(1, size(cmat, 1))
    n = size(cmat, 2)
    call blas_herk(ula, tra, n, k, a, amat, lda, b, cmat, ldc)
  end subroutine blas_herk_cdp


  pure subroutine blas_syr2k_rsp(amat, bmat, cmat, uplo, trans, alpha, beta)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: bmat(:, :)
    real(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula, tra
    real(wp) :: a, b
    integer(ik) :: n, k, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    ula = 'u'
    if (present(uplo)) ula = uplo
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    n = size(cmat, 2)
    call blas_syr2k(ula, tra, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_syr2k_rsp

  pure subroutine blas_syr2k_csp(amat, bmat, cmat, uplo, trans, alpha, beta)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: bmat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: ula, tra
    complex(wp) :: a, b
    integer(ik) :: n, k, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    ula = 'u'
    if (present(uplo)) ula = uplo
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    n = size(cmat, 2)
    call blas_syr2k(ula, tra, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_syr2k_csp

  pure subroutine blas_her2k_csp(amat, bmat, cmat, uplo, trans, alpha, beta)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: bmat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    complex(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula, tra
    complex(wp) :: a
    real(wp) :: b
    integer(ik) :: n, k, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    ula = 'u'
    if (present(uplo)) ula = uplo
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    n = size(cmat, 2)
    call blas_her2k(ula, tra, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_her2k_csp

  pure subroutine blas_syr2k_rdp(amat, bmat, cmat, uplo, trans, alpha, beta)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: bmat(:, :)
    real(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula, tra
    real(wp) :: a, b
    integer(ik) :: n, k, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    ula = 'u'
    if (present(uplo)) ula = uplo
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    n = size(cmat, 2)
    call blas_syr2k(ula, tra, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_syr2k_rdp

  pure subroutine blas_syr2k_cdp(amat, bmat, cmat, uplo, trans, alpha, beta)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: bmat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: ula, tra
    complex(wp) :: a, b
    integer(ik) :: n, k, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    ula = 'u'
    if (present(uplo)) ula = uplo
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    n = size(cmat, 2)
    call blas_syr2k(ula, tra, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_syr2k_cdp

  pure subroutine blas_her2k_cdp(amat, bmat, cmat, uplo, trans, alpha, beta)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: bmat(:, :)
    complex(wp), contiguous, intent(inout) :: cmat(:, :)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    complex(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula, tra
    complex(wp) :: a
    real(wp) :: b
    integer(ik) :: n, k, lda, ldb, ldc

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    ula = 'u'
    if (present(uplo)) ula = uplo
    if ((tra == 'n'.or.tra == 'N')) then
      k = size(amat, 2)
    else
      k = size(amat, 1)
    end if
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    ldc = max(1, size(cmat, 1))
    n = size(cmat, 2)
    call blas_her2k(ula, tra, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
  end subroutine blas_her2k_cdp


  pure subroutine blas_trsm_rsp(amat, bmat, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(inout) :: bmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    character, intent(in), optional :: transa
    character, intent(in), optional :: diag
    real(wp), intent(in), optional :: alpha

    character :: lra, ula, tra, dia
    real(wp) :: a
    integer(ik) :: m, n, lda, ldb

    a = 1.0_wp
    dia = 'n'
    lra = 'l'
    tra = 'n'
    ula = 'u'
    if (present(alpha)) a = alpha
    if (present(diag)) dia = diag
    if (present(side)) lra = side
    if (present(transa)) tra = transa
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    m = size(bmat, 1)
    n = size(bmat, 2)
    call blas_trsm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
  end subroutine blas_trsm_rsp

  pure subroutine blas_trsm_csp(amat, bmat, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(inout) :: bmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    character, intent(in), optional :: transa
    character, intent(in), optional :: diag
    complex(wp), intent(in), optional :: alpha

    character :: lra, ula, tra, dia
    complex(wp) :: a
    integer(ik) :: m, n, lda, ldb

    a = 1.0_wp
    dia = 'n'
    lra = 'l'
    tra = 'n'
    ula = 'u'
    if (present(alpha)) a = alpha
    if (present(diag)) dia = diag
    if (present(side)) lra = side
    if (present(transa)) tra = transa
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    m = size(bmat, 1)
    n = size(bmat, 2)
    call blas_trsm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
  end subroutine blas_trsm_csp

  pure subroutine blas_trsm_rdp(amat, bmat, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(inout) :: bmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    character, intent(in), optional :: transa
    character, intent(in), optional :: diag
    real(wp), intent(in), optional :: alpha

    character :: lra, ula, tra, dia
    real(wp) :: a
    integer(ik) :: m, n, lda, ldb

    a = 1.0_wp
    dia = 'n'
    lra = 'l'
    tra = 'n'
    ula = 'u'
    if (present(alpha)) a = alpha
    if (present(diag)) dia = diag
    if (present(side)) lra = side
    if (present(transa)) tra = transa
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    m = size(bmat, 1)
    n = size(bmat, 2)
    call blas_trsm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
  end subroutine blas_trsm_rdp

  pure subroutine blas_trsm_cdp(amat, bmat, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(inout) :: bmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    character, intent(in), optional :: transa
    character, intent(in), optional :: diag
    complex(wp), intent(in), optional :: alpha

    character :: lra, ula, tra, dia
    complex(wp) :: a
    integer(ik) :: m, n, lda, ldb

    a = 1.0_wp
    dia = 'n'
    lra = 'l'
    tra = 'n'
    ula = 'u'
    if (present(alpha)) a = alpha
    if (present(diag)) dia = diag
    if (present(side)) lra = side
    if (present(transa)) tra = transa
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    m = size(bmat, 1)
    n = size(bmat, 2)
    call blas_trsm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
  end subroutine blas_trsm_cdp


  pure subroutine blas_trmm_rsp(amat, bmat, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(inout) :: bmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    character, intent(in), optional :: transa
    character, intent(in), optional :: diag
    real(wp), intent(in), optional :: alpha

    character :: lra, ula, tra, dia
    real(wp) :: a
    integer(ik) :: m, n, lda, ldb

    a = 1.0_wp
    dia = 'n'
    lra = 'l'
    tra = 'n'
    ula = 'u'
    if (present(alpha)) a = alpha
    if (present(diag)) dia = diag
    if (present(side)) lra = side
    if (present(transa)) tra = transa
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    m = size(bmat, 1)
    n = size(bmat, 2)
    call blas_trmm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
  end subroutine blas_trmm_rsp

  pure subroutine blas_trmm_csp(amat, bmat, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(inout) :: bmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    character, intent(in), optional :: transa
    character, intent(in), optional :: diag
    complex(wp), intent(in), optional :: alpha

    character :: lra, ula, tra, dia
    complex(wp) :: a
    integer(ik) :: m, n, lda, ldb

    a = 1.0_wp
    dia = 'n'
    lra = 'l'
    tra = 'n'
    ula = 'u'
    if (present(alpha)) a = alpha
    if (present(diag)) dia = diag
    if (present(side)) lra = side
    if (present(transa)) tra = transa
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    m = size(bmat, 1)
    n = size(bmat, 2)
    call blas_trmm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
  end subroutine blas_trmm_csp

  pure subroutine blas_trmm_rdp(amat, bmat, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(inout) :: bmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    character, intent(in), optional :: transa
    character, intent(in), optional :: diag
    real(wp), intent(in), optional :: alpha

    character :: lra, ula, tra, dia
    real(wp) :: a
    integer(ik) :: m, n, lda, ldb

    a = 1.0_wp
    dia = 'n'
    lra = 'l'
    tra = 'n'
    ula = 'u'
    if (present(alpha)) a = alpha
    if (present(diag)) dia = diag
    if (present(side)) lra = side
    if (present(transa)) tra = transa
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    m = size(bmat, 1)
    n = size(bmat, 2)
    call blas_trmm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
  end subroutine blas_trmm_rdp

  pure subroutine blas_trmm_cdp(amat, bmat, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(inout) :: bmat(:, :)
    character, intent(in), optional :: side
    character, intent(in), optional :: uplo
    character, intent(in), optional :: transa
    character, intent(in), optional :: diag
    complex(wp), intent(in), optional :: alpha

    character :: lra, ula, tra, dia
    complex(wp) :: a
    integer(ik) :: m, n, lda, ldb

    a = 1.0_wp
    dia = 'n'
    lra = 'l'
    tra = 'n'
    ula = 'u'
    if (present(alpha)) a = alpha
    if (present(diag)) dia = diag
    if (present(side)) lra = side
    if (present(transa)) tra = transa
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    m = size(bmat, 1)
    n = size(bmat, 2)
    call blas_trmm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
  end subroutine blas_trmm_cdp


end module blas_pure_level3
