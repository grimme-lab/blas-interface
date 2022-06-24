! This file is part of blas-interface.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

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
