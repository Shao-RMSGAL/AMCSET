! matrix_multiply.f90
subroutine matrix_multiply(a, b, c) bind(C)
    implicit none
    real(8), intent(in) :: a(6,6), b(6,6)
    real(8), intent(out) :: c(6,6)
    integer :: i, j, k

    c = matmul(a,b)
end subroutine matrix_multiply

! initialize_identity.f90
subroutine initialize_identity(a) bind(C)
    implicit none
    real(8), intent(out) :: a(6,6)
    integer :: i, j

    a = 0.0
    do i = 1, 6
        a(i,i) = 2.0
    end do


end subroutine initialize_identity