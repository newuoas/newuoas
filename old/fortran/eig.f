      SUBROUTINE EIG(A, N, K, LAMBDA, V)
         IMPLICIT NONE
         INTEGER(kind=4), INTENT(in) :: N, K
         REAL(kind=8), INTENT(in) :: A(N,N) 
         REAL(kind=8), INTENT(OUT) :: LAMBDA(K), V(N,K) 

         REAL(kind=8) :: WORK(1+6*N+2*N**2),TEMPA(N,N), W(N)
         INTEGER(kind=4) :: IWORK(3+5*N) 
         INTEGER(kind=4) :: LWORK, LIWORK, INFO, I
         CHARACTER :: JOBZ, UPLO
         JOBZ = 'V'
         UPLO = 'U'
         LWORK = 1+6*N+2*N**2
         LIWORK = 3+5*N
         TEMPA = A
C         write(*,*) N, K
         CALL DSYEVD( JOBZ, UPLO, N, TEMPA, N, W, WORK, LWORK, IWORK,
     1   LIWORK, INFO )
         IF (INFO /= 0) write(*,*) "ERROR IN EIG. INFO=", INFO
         DO I = 1, K
            LAMBDA(I) = W(N-I+1)
            V(:,I) = TEMPA(:, N-I+1)
         END DO
      END
