module blas_utils
  use blas_kinds, only : ik, lk
  implicit none
  private

  public :: lsame, xerbla, xerbla_array

  interface
    function lsame(ca, cb)
      import :: lk
      character(len=1) :: ca
      character(len=1) :: cb
      logical(lk) :: lsame
    end function lsame
  end interface

  interface
    subroutine xerbla(srname, info)
      import :: ik
      character(*) :: srname
      integer(ik) :: info
    end subroutine xerbla
  end interface

  interface
    subroutine xerbla_array(srname_array, srname_len, info)
      import :: ik
      integer(ik) :: srname_len
      character(len=1) :: srname_array(srname_len)
      integer(ik) :: info
    end subroutine xerbla_array
  end interface

end module blas_utils
