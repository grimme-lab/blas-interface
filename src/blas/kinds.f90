module blas_kinds
  implicit none

  integer, parameter :: sp = selected_real_kind(7)

  integer, parameter :: dp = selected_real_kind(16)

  integer, parameter :: ik = kind(1)

  integer, parameter :: lk = kind(.true.)

end module blas_kinds
