subroutine diagonalization(matrix, dimen, spectrum)
implicit none

integer dimen,LWORK,INFO
real*8 WORK(3*dimen-1)
real*8 matrix(dimen,dimen),spectrum(dimen),matrix1(dimen,dimen)


matrix1 = matrix
LWORK = 3*dimen-1

call DSYEV( 'V', 'U', dimen, matrix1, dimen , spectrum, WORK, LWORK, INFO)

end subroutine
