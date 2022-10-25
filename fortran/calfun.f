      SUBROUTINE CALFUN(N,X,F)
        USE PROBLEM
        USE PERMUTE
        USE STATISTICS
        IMPLICIT NONE
        INTEGER(KIND=4), INTENT(IN) :: N
        REAL(KIND=8), INTENT(IN) :: X(N)
        REAL(KIND=8), INTENT(OUT) :: F

        INTEGER(KIND=4) :: I, J, K, M, J1, J2, J3, J4, J5, IW, ML, MU
        REAL(KIND=8) :: TEMP, TEMP1, TEMP2, TEMP3, SUM
        REAL(KIND=8) :: YCHEBQ(N+1,N), DP1, DP2, PD, R, S
        REAL(KIND=8) :: ALPHA, BETA, GAMMA, DELTA, F1, F2, F3, F4, F5 
        REAL(KIND=8) :: Y(0:N+1), h, T(N)
        REAL(KIND=8) :: C, BETAI, BETAJ
        REAL(KIND=8) :: ZETA, KAPPA, p, PI, SCALE
        REAL(KIND=8) :: YTEMP(N), PAR(N) 

        Y(0) = 0.0D0
        Y(N+1) = 0.0D0
        DO I = 1, N
            Y(I) = X(PERM(I))
        END DO
        F = 0.0D0

        
        IF (TRIM(PROB) .EQ. 'ARGLINA') THEN
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 32, also CUTEr.
            M = 2*N
            TEMP = 0.0D0
            DO I = 1, N
                TEMP = TEMP+Y(I)
            END DO
            F = 0.0D0
            DO I = 1, N
                F = F+(Y(I)-2.0D0/DFLOAT(M)*TEMP-1.0D0)**2
            END DO
            DO I = N+1, M 
                F = F+(-2.0D0/DFLOAT(M)*TEMP-1.0D0)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'ARGLINA4') THEN
CCCC !!!!!! DIFFERENT FROM ARGLINA. Quartic instead of  Quadratic. !!!!
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 32, also CUTEr.
            M = 2*N
            TEMP = 0.0D0
            DO I = 1, N
                TEMP = TEMP+Y(I)
            END DO
            F = 0.0D0
            DO I = 1, N
                F = F+(Y(I)-2.0D0/DFLOAT(M)*TEMP-1.0D0)**4
            END DO
            DO I = N+1, M 
                F = F+(-2.0D0/DFLOAT(M)*TEMP-1.0D0)**4
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'ARGLINB') THEN
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 33, also CUTEr.
            M = 2*N
            TEMP = 0.0D0
            DO I = 1, N
                TEMP = TEMP+DFLOAT(I)*Y(I)
            END DO
            F = 0.0D0
            DO I = 1,M 
                F = F + (DFLOAT(I)*TEMP-1.0D0)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'ARGLINC') THEN
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 34, also CUTEr.
            M = 2*N
            TEMP = 0.0D0
            DO I = 2, N-1
                TEMP = TEMP+DFLOAT(I)*Y(I)
            END DO
            F = 1.0D0
            DO I = 2,M-1 
                F = F + (DFLOAT(I-1)*TEMP-1.0D0)**2
            END DO
            F = F + 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'ARGTRIG') THEN
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 26, also CUTEr.
            F = 0.0D0
            DO I = 1, N
                TEMP = DFLOAT(N)+DFLOAT(I)*(1.0D0-COS(Y(I)))-SIN(Y(I))
                DO J = 1, N
                    TEMP = TEMP - COS(Y(J))
                END DO
                F = F + TEMP*TEMP
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'ARWHEAD') THEN
            F = 0.0D0
            do i = 1,n-1
                F = F + (Y(I)*Y(I)+Y(N)*Y(N))**2 - 4.0D0*Y(I) +3.0D0
            end do
        ELSE IF (TRIM(PROB) .EQ. 'BDQRTIC') THEN
C See Powell, 200, On trust region methods for unconstrained
C minimization without derivatives, Problem 2. Also see CUTEr.
C NOTICE: the definitions in the above references are different,
C the difference being -4Y(I)+3 V.S. (-4Y(I)+3)^2. We call the former 
C one BDQRTICP (see below). The original definition
C is in Conn, Gould Lescrenoer, Toint, 1994, Performance of a 
C Multifrontal scheme for partially seperable optimization, Problem
C 61. But the paper is not available for now.
            F = 0.0D0
            do i = 1,n-4
                F = F + (Y(I)**2+2.0D0*Y(I+1)**2+3.0D0*Y(I+2)**2
     1                + 4.0D0*Y(I+3)**2+5.0D0*Y(N)**2)**2 
     1                + (- 4.0D0*Y(I) + 3.0D0)**2
            end do
        ELSE IF (TRIM(PROB) .EQ. 'BDQRTICP') THEN
C See Powell, 200, On trust region methods for unconstrained
C minimization without derivatives, Problem 2. Also see CUTEr.
C NOTICE: the definitions in the above references are different,
C the difference being -4Y(I)+3 V.S. (-4Y(I)+3)^2. Here is Powell's
C version.
            F = 0.0D0
            do i = 1,n-4
                F = F + (Y(I)**2+2.0D0*Y(I+1)**2+3.0D0*Y(I+2)**2
     1                + 4.0D0*Y(I+3)**2+5.0D0*Y(N)**2)**2 
     1                - 4.0D0*Y(I) + 3.0D0
            end do
        ELSE IF (TRIM(PROB) .EQ. 'BDVALUE') THEN
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 28, also CUTEr.
            F = 0.0D0
            h = 1.0D0/DFLOAT(N+1) 
            DO I = 1, N
            F = F+(2.0D0*Y(I)-Y(I-1)-Y(I+1)
     1      +0.5D0*h*h*(Y(I)+DFLOAT(I)/DFLOAT(N+1)+1.0D0)**3)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'BROWNAL') THEN
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 27, also CUTEr.
            F = 0.0D0
            TEMP1 = 0.0D0
            TEMP2 = 1.0D0
            DO I = 1, N
                TEMP1 = TEMP1 + Y(I)
                TEMP2 = TEMP2*Y(I)
            END DO
            F = (TEMP2 - 1.0D0)**2
            DO I = 1, N-1
                F = F + (Y(I) + TEMP1 - DFLOAT(N+1))**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'BROYDN3D') THEN
C See CUTEr.
            F = 0.0D0
            DO I = 1, N
            F = F+((3.0D0-2.0D0*Y(I))*Y(I)-Y(I-1)-2*Y(I+1)+1.0D0)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'BROYDN7D') THEN
C See Ph. L. Toint "Some Numerical Results Using a Sparse Matrix
C Updating Formula in Unconstrained Optimization", Problem 3.4, also
C CUTEr. 
            F = 0.0D0
            DO I = 1, N
            F =F+DABS(Y(I-1)-Y(I)*(3.0D0-0.5D0*Y(I))+2.0D0*Y(I+1)-1.0D0)
     1       **(7.0D0/3.0D0)
            END DO
            DO I = 1, N/2
                F = F + DABS(Y(I)+Y(I+N/2))**(7.0D0/3.0D0)
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'BRYBND') THEN
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 31, also CUTEr.
            ML = 5
            MU = 1
            F = 0.0D0
            DO I = 1, N
    
                TEMP = 0.0D0
                DO J = MAX(1, I-ML), MIN(N,I+MU) 
                    IF (J/=I) THEN 
                        TEMP = TEMP+Y(J)*(1.0D0+Y(J))
                    END IF
                END DO
               F = F + (Y(I)*(2.0D0+5.0D00*Y(I)*Y(I))+1.0D0-TEMP)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'CHAINWOO') THEN
C See CUTEr.
            F = 1.0D0
            DO I = 1, N/2-1
                J = 2*I
                F = F + 1.0D2*(Y(J)-Y(J-1)**2)**2+(1.0D0-Y(J-1))**2
     1            + 9.0D1*(Y(J+2)-Y(J+1)**2)**2 + (1.0D0-Y(J+1))**2
     1            + 1.0D1*(Y(J)+Y(J+2)-2.0D0)**2 +
     1            1.0D-1*(Y(J)-Y(J+2))**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'CHEBQUAD') THEN
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 35, also CUTEr.
            DO 10 J=1,N
            YCHEBQ(1,J)=1.0D0
   10       YCHEBQ(2,J)=2.0D0*Y(J)-1.0D0
            DO 20 I=2,N
            DO 20 J=1,N
   20       YCHEBQ(I+1,J)=2.0D0*YCHEBQ(2,J)*YCHEBQ(I,J)-YCHEBQ(I-1,J)
            F=0.0D0
            IW=1
            DO 40 I=1,N+1
            SUM=0.0D0
            DO 30 J=1,N
   30       SUM=SUM+YCHEBQ(I,J)
            SUM=SUM/DFLOAT(N)
            IF (IW .GT. 0) SUM=SUM+1.0D0/DFLOAT(I*I-2*I)
            IW=-IW
   40       F=F+SUM*SUM
        ELSE IF (TRIM(PROB) .EQ. 'CHPOWELLB') THEN
C Chained version of Powell Badly Scaled Function.
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 3.
            F = 0.0D0
            DO I = 1, N-1 
                F = F + (1.0D4*Y(I)*Y(I+1)-1.0D0)**2 +
     1          (EXP(-Y(I))+EXP(-Y(I+1)) - 1.0001D0)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'CHPOWELLS') THEN
C Chained version of Powell Singular Function. 
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 13, and Ying-Jie Li, Dong-Hui Li, Truncated regularized Newton
C Method for Convex Minimization, Problem 6, and Luksan, Vicek,
C Sparse and partially separable test problems for unconstrained and
C equality constrained optimization.
            F = 0.0D0
            DO J = 1, (N-2)/2
                I = 2*J
                F = F + (Y(I-1)+10.0D0*Y(I))**2 +
     1          5.0D0*(Y(I+1)-Y(I+2))**2 +
     1          (Y(I)-2.0D0*Y(I+1))**4 +
     1          10.0D0*(Y(I-1)-Y(I+2))**4 
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'CHNROSNB') THEN
C See CUTEr. MODIFIED by Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." It is another version of chained
C ROSENBROCK. Compare with CHROSEN and ROSENBROCK.
            ALPHA = 16.0D0*(1.5D0+SIN(DFLOAT(I)))**2
            F = 0.0D0
            DO I = 2, N
                F = F + ALPHA*(Y(I-1)-Y(I)**2)**2+(Y(I)-1)**2 
            END DO

        ELSE IF (TRIM(PROB) .EQ. 'CHROSEN') THEN
            F = 0.0D0
            do i = 1,n-1
                F = F + (4.0D0)*(Y(I)-Y(I+1)*Y(I+1))*
     1       (Y(I)-Y(I+1)*Y(I+1)) +
     1       (1.0D0-Y(I+1))*(1.0D0-Y(I+1))
            end do
        ELSE IF (TRIM(PROB) .EQ. 'COSINE') THEN
C See CUTEr.
            F = 0.0D0
            DO I = 1, N-1
                F = F + COS(Y(I)**2-0.5D0*Y(I+1))
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'CRAGGLVY') THEN
C See CUTEr.
            F = 0.0D0
            DO J = 1, (N-2)/2
                I = 2*J
                F = F + (EXP(Y(I-1))-Y(I))**4
     1                + 1.0D2*(Y(I)-Y(I+1))**6
     1                + (TAN(Y(I+1)-Y(I+2))+Y(I+1)-Y(I+2))**4
     1                + Y(I-1)**8 + (Y(I+2)-1.0D0)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'CUBE') THEN
C See CUTEr.
            F = (Y(1)-1.0D0)**2
            DO I = 2, N
                F = F + 100.0D0 * (Y(I)-Y(I-1)**3)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'CURLY10' .OR. 
     1   TRIM(PROB) .EQ. 'CURLY20' .OR. TRIM(PROB) .EQ. 'CURLY30') THEN
C See CUTEr.
            IF (TRIM(PROB) .EQ. 'CURLY10') THEN
                K = 10
            ELSE IF (TRIM(PROB) .EQ. 'CURLY20') THEN
                K = 20
            ELSE
                K = 30 
            END IF
            F = 0.0D0
            DO I = 1, N
                TEMP = 0.0D0
                DO J = I, MIN(I+K,N)
                    TEMP = TEMP + Y(J)
                END DO
             F = F + TEMP*(TEMP*(TEMP**2-2.0D1)-1.0D-1)
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'DIXMAANE') THEN
            M = N/3
            ALPHA = 1.0D0
!!!            BETA = 0.0D0
            GAMMA = 0.125D0
            DELTA = 0.125D0
            F1 = 0.0D0
            DO I = 1, N
                F1 = F1 + DFLOAT(I)/DFLOAT(N)*Y(I)**2
            END DO
            F1 = ALPHA*F1
!!!            F2 = 0.0D0
!!!            DO I = 1, N-1
!!!                F2 = F2 + Y(I)**2*(Y(I+1)+Y(I+1)**2)**2
!!!            END DO
!!!            F2 = BETA*F2
            F3 = 0.0D0
            DO I = 1, 2*M
                F3 = F3 + Y(I)**2*Y(I+M)**4
            END DO
            F3 = GAMMA*F3
            F4 = 0.0D0
            DO I = 1, M
                F4 = F4 + DFLOAT(I)/DFLOAT(N)*Y(I)*Y(I+2*M)
            END DO
            F4 = DELTA*F4
!!!            F = 1.0D0+F1+F2+F3+F4
            F = 1.0D0+F1+F3+F4
        ELSE IF (TRIM(PROB) .EQ. 'DIXMAANF') THEN
            M = N/3
            ALPHA = 1.0D0
            BETA = 0.0625D0
            GAMMA = 0.0625D0
            DELTA = 0.0625D0
            F1 = 0.0D0
            DO I = 1, N
                F1 = F1 + DFLOAT(I)/DFLOAT(N)*Y(I)**2
            END DO
            F1 = ALPHA*F1
            F2 = 0.0D0
            DO I = 1, N-1
                F2 = F2 + Y(I)**2*(Y(I+1)+Y(I+1)**2)**2
            END DO
            F2 = BETA*F2
            F3 = 0.0D0
            DO I = 1, 2*M
                F3 = F3 + Y(I)**2*Y(I+M)**4
            END DO
            F3 = GAMMA*F3
            F4 = 0.0D0
            DO I = 1, M
                F4 = F4 + DFLOAT(I)/DFLOAT(N)*Y(I)*Y(I+2*M)
            END DO
            F4 = DELTA*F4
            F = 1.0D0+F1+F2+F3+F4

        ELSE IF (TRIM(PROB) .EQ. 'DIXMAANG') THEN
            M = N/3
            ALPHA = 1.0D0
            BETA = 0.125D0
            GAMMA = 0.125D0
            DELTA = 0.125D0
            F1 = 0.0D0
            DO I = 1, N
                F1 = F1 + DFLOAT(I)/DFLOAT(N)*Y(I)**2
            END DO
            F1 = ALPHA*F1
            F2 = 0.0D0
            DO I = 1, N-1
                F2 = F2 + Y(I)**2*(Y(I+1)+Y(I+1)**2)**2
            END DO
            F2 = BETA*F2
            F3 = 0.0D0
            DO I = 1, 2*M
                F3 = F3 + Y(I)**2*Y(I+M)**4
            END DO
            F3 = GAMMA*F3
            F4 = 0.0D0
            DO I = 1, M
                F4 = F4 + DFLOAT(I)/DFLOAT(N)*Y(I)*Y(I+2*M)
            END DO
            F4 = DELTA*F4
            F = 1.0D0+F1+F2+F3+F4
        ELSE IF (TRIM(PROB) .EQ. 'DIXMAANH') THEN
            M = N/3
            ALPHA = 1.0D0
            BETA = 0.26D0
            GAMMA = 0.26D0
            DELTA = 0.26D0
            F1 = 0.0D0
            DO I = 1, N
                F1 = F1 + DFLOAT(I)/DFLOAT(N)*Y(I)**2
            END DO
            F1 = ALPHA*F1
            F2 = 0.0D0
            DO I = 1, N-1
                F2 = F2 + Y(I)**2*(Y(I+1)+Y(I+1)**2)**2
            END DO
            F2 = BETA*F2
            F3 = 0.0D0
            DO I = 1, 2*M
                F3 = F3 + Y(I)**2*Y(I+M)**4
            END DO
            F3 = GAMMA*F3
            F4 = 0.0D0
            DO I = 1, M
                F4 = F4 + DFLOAT(I)/DFLOAT(N)*Y(I)*Y(I+2*M)
            END DO
            F4 = DELTA*F4
            F = 1.0D0+F1+F2+F3+F4
        ELSE IF (TRIM(PROB) .EQ. 'DIXMAANI') THEN
            M = N/3
            ALPHA = 1.0D0
!!!            BETA = 0.0D0
            GAMMA = 0.125D0
            DELTA = 0.125D0
            F1 = 0.0D0
            DO I = 1, N
                F1 = F1 + (DFLOAT(I)/DFLOAT(N)*Y(I))**2
            END DO
            F1 = ALPHA*F1
!!!            F2 = 0.0D0
!!!            DO I = 1, N-1
!!!                F2 = F2 + Y(I)**2*(Y(I+1)+Y(I+1)**2)**2
!!!            END DO
!!!            F2 = BETA*F2
            F3 = 0.0D0
            DO I = 1, 2*M
                F3 = F3 + Y(I)**2*Y(I+M)**4
            END DO
            F3 = GAMMA*F3
            F4 = 0.0D0
            DO I = 1, M
                F4 = F4 + (DFLOAT(I)/DFLOAT(N))**2*Y(I)*Y(I+2*M)
            END DO
            F4 = DELTA*F4
!!!            F = 1.0D0+F1+F2+F3+F4
            F = 1.0D0+F1+F3+F4
        ELSE IF (TRIM(PROB) .EQ. 'DIXMAANJ') THEN
            M = N/3
            ALPHA = 1.0D0
            BETA = 0.0625D0
            GAMMA = 0.0625D0
            DELTA = 0.0625D0
            F1 = 0.0D0
            DO I = 1, N
                F1 = F1 + (DFLOAT(I)/DFLOAT(N)*Y(I))**2
            END DO
            F1 = ALPHA*F1
            F2 = 0.0D0
            DO I = 1, N-1
                F2 = F2 + Y(I)**2*(Y(I+1)+Y(I+1)**2)**2
            END DO
            F2 = BETA*F2
            F3 = 0.0D0
            DO I = 1, 2*M
                F3 = F3 + Y(I)**2*Y(I+M)**4
            END DO
            F3 = GAMMA*F3
            F4 = 0.0D0
            DO I = 1, M
                F4 = F4 + (DFLOAT(I)/DFLOAT(N))**2*Y(I)*Y(I+2*M)
            END DO
            F4 = DELTA*F4
            F = 1.0D0+F1+F2+F3+F4
        ELSE IF (TRIM(PROB) .EQ. 'DIXMAANK') THEN
            M = N/3
            ALPHA = 1.0D0
            BETA = 0.125D0
            GAMMA = 0.125D0
            DELTA = 0.125D0
            F1 = 0.0D0
            DO I = 1, N
                F1 = F1 + (DFLOAT(I)/DFLOAT(N)*Y(I))**2
            END DO
            F1 = ALPHA*F1
            F2 = 0.0D0
            DO I = 1, N-1
                F2 = F2 + Y(I)**2*(Y(I+1)+Y(I+1)**2)**2
            END DO
            F2 = BETA*F2
            F3 = 0.0D0
            DO I = 1, 2*M
                F3 = F3 + Y(I)**2*Y(I+M)**4
            END DO
            F3 = GAMMA*F3
            F4 = 0.0D0
            DO I = 1, M
                F4 = F4 + (DFLOAT(I)/DFLOAT(N))**2*Y(I)*Y(I+2*M)
            END DO
            F4 = DELTA*F4
            F = 1.0D0+F1+F2+F3+F4
        ELSE IF (TRIM(PROB) .EQ. 'DIXMAANL') THEN
            M = N/3
            ALPHA = 1.0D0
            BETA = 0.26D0
            GAMMA = 0.26D0
            DELTA = 0.26D0
            F1 = 0.0D0
            DO I = 1, N
                F1 = F1 + (DFLOAT(I)/DFLOAT(N)*Y(I))**2
            END DO
            F1 = ALPHA*F1
            F2 = 0.0D0
            DO I = 1, N-1
                F2 = F2 + Y(I)**2*(Y(I+1)+Y(I+1)**2)**2
            END DO
            F2 = BETA*F2
            F3 = 0.0D0
            DO I = 1, 2*M
                F3 = F3 + Y(I)**2*Y(I+M)**4
            END DO
            F3 = GAMMA*F3
            F4 = 0.0D0
            DO I = 1, M
                F4 = F4 + (DFLOAT(I)/DFLOAT(N))**2*Y(I)*Y(I+2*M)
            END DO
            F4 = DELTA*F4
            F = 1.0D0+F1+F2+F3+F4
        ELSE IF (TRIM(PROB) .EQ. 'DIXMAANM') THEN
            M = N/3
            ALPHA = 1.0D0
!!!            BETA = 0.0D0
            GAMMA = 0.125D0
            DELTA = 0.125D0
            F1 = 0.0D0
            DO I = 1, N
                F1 = F1 + (DFLOAT(I)/DFLOAT(N)*Y(I))**2
            END DO
            F1 = ALPHA*F1
!!!            F2 = 0.0D0
!!!            DO I = 1, N-1
!!!                F2 = F2 + (DFLOAT(I)/DFLOAT(N))*Y(I)**2*
!!!     1           (Y(I+1)+Y(I+1)**2)**2
!!!            END DO
!!!            F2 = BETA*F2
            F3 = 0.0D0
            DO I = 1, 2*M
                F3 = F3 + (DFLOAT(I)/DFLOAT(N))*Y(I)**2*Y(I+M)**4
            END DO
            F3 = GAMMA*F3
            F4 = 0.0D0
            DO I = 1, M
                F4 = F4 + (DFLOAT(I)/DFLOAT(N))**2*Y(I)*Y(I+2*M)
            END DO
            F4 = DELTA*F4
!!!            F = 1.0D0+F1+F2+F3+F4
            F = 1.0D0+F1+F3+F4
        ELSE IF (TRIM(PROB) .EQ. 'DIXMAANN') THEN
            M = N/3
            ALPHA = 1.0D0
            BETA = 0.0625D0
            GAMMA = 0.0625D0
            DELTA = 0.0625D0
            F1 = 0.0D0
            DO I = 1, N
                F1 = F1 + (DFLOAT(I)/DFLOAT(N)*Y(I))**2
            END DO
            F1 = ALPHA*F1
            F2 = 0.0D0
            DO I = 1, N-1
                F2 = F2 + (DFLOAT(I)/DFLOAT(N))*Y(I)**2*
     1           (Y(I+1)+Y(I+1)**2)**2
            END DO
            F2 = BETA*F2
            F3 = 0.0D0
            DO I = 1, 2*M
                F3 = F3 + (DFLOAT(I)/DFLOAT(N))*Y(I)**2*Y(I+M)**4
            END DO
            F3 = GAMMA*F3
            F4 = 0.0D0
            DO I = 1, M
                F4 = F4 + (DFLOAT(I)/DFLOAT(N))**2*Y(I)*Y(I+2*M)
            END DO
            F4 = DELTA*F4
            F = 1.0D0+F1+F2+F3+F4
        ELSE IF (TRIM(PROB) .EQ. 'DIXMAANO') THEN
            M = N/3
            ALPHA = 1.0D0
            BETA = 0.125D0
            GAMMA = 0.125D0
            DELTA = 0.125D0
            F1 = 0.0D0
            DO I = 1, N
                F1 = F1 + (DFLOAT(I)/DFLOAT(N)*Y(I))**2
            END DO
            F1 = ALPHA*F1
            F2 = 0.0D0
            DO I = 1, N-1
                F2 = F2 + (DFLOAT(I)/DFLOAT(N))*Y(I)**2*
     1           (Y(I+1)+Y(I+1)**2)**2
            END DO
            F2 = BETA*F2
            F3 = 0.0D0
            DO I = 1, 2*M
                F3 = F3 + (DFLOAT(I)/DFLOAT(N))*Y(I)**2*Y(I+M)**4
            END DO
            F3 = GAMMA*F3
            F4 = 0.0D0
            DO I = 1, M
                F4 = F4 + (DFLOAT(I)/DFLOAT(N))**2*Y(I)*Y(I+2*M)
            END DO
            F4 = DELTA*F4
            F = 1.0D0+F1+F2+F3+F4
        ELSE IF (TRIM(PROB) .EQ. 'DIXMAANP') THEN
            M = N/3
            ALPHA = 1.0D0
            BETA = 0.26D0
            GAMMA = 0.26D0
            DELTA = 0.26D0
            F1 = 0.0D0
            DO I = 1, N
                F1 = F1 + (DFLOAT(I)/DFLOAT(N)*Y(I))**2
            END DO
            F1 = ALPHA*F1
            F2 = 0.0D0
            DO I = 1, N-1
                F2 = F2 + (DFLOAT(I)/DFLOAT(N))*Y(I)**2*
     1           (Y(I+1)+Y(I+1)**2)**2
            END DO
            F2 = BETA*F2
            F3 = 0.0D0
            DO I = 1, 2*M
                F3 = F3 + (DFLOAT(I)/DFLOAT(N))*Y(I)**2*Y(I+M)**4
            END DO
            F3 = GAMMA*F3
            F4 = 0.0D0
            DO I = 1, M
                F4 = F4 + (DFLOAT(I)/DFLOAT(N))**2*Y(I)*Y(I+2*M)
            END DO
            F4 = DELTA*F4
            F = 1.0D0+F1+F2+F3+F4
        ELSE IF (TRIM(PROB) .EQ. 'DQRTIC') THEN
C See CUTEr.
            F = 0.0D0
            DO I = 1, N
                F = F + (Y(I)-DFLOAT(I))**4
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'EDENSCH') THEN
C See CUTEr.
            F = 16.0D0
            DO I = 1, N-1 
                F = F +
     1          (Y(I)-2.0D0)**4+
     1          (Y(I)*Y(I+1)-2.0D0*Y(I+1))**2 +
     1          (Y(I+1)+1.0D0)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'EG2') THEN
C See CUTEr.
            F = 0.0D0
            DO I = 1, N-1
                F = F + SIN(Y(1)+Y(I)**2-1.0D0)
            END DO
            F = F + 0.5D0*SIN(Y(N)**2)
        ELSE IF (TRIM(PROB) .EQ. 'ENGVAL1') THEN
C See CUTEr. Compare with ARWHEAD.
            F = 0.0D0
            DO I = 1, N-1
                F = F + (Y(I)**2+Y(I+1)**2)**2 - 4.0D0*Y(I) + 3.0D0
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'ERRINROS') THEN
C See CUTEr. MODIFIED by Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." It is another version of chained
C ROSENBROCK. Compare with CHNROSNB, CHROSEN, and ROSENBROCK.
            ALPHA = 16.0D0*(1.5D0+SIN(DFLOAT(I)))**2
            F = 0.0D0
            DO I = 2, N
                F = F + (Y(I-1)-ALPHA*Y(I)**2)**2+(Y(I)-1)**2 
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'EXTROSNB') THEN
C See CUTEr and Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." It is another version of chained
C ROSENBROCK. Compare with CHNROSNB, CHROSEN, ERRINROS, and ROSENBROCK.
            F = (Y(1)-1.0D0)**2
            DO I = 2, N
                F = F + 1.0D2*(Y(I)-Y(I-1)**2)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'EXTTET') THEN
C See N. Andrei, 2008, An Unconstrained Optimization Test Functions
C Collection, Extended Three Term Exponentials.
            F = 0.0D0
            DO J = 1, N/2
                I = 2*J
                F = F + EXP(Y(I-1)+3.0D0*Y(I)-0.1D0)
     1                + EXP(Y(I-1)-3.0D0*Y(I)-0.1D0)
     1                + EXP(-Y(I-1)-0.1D0)
                END DO

        ELSE IF (TRIM(PROB) .EQ. 'FIROSE') THEN
C Five-Diagonal ROSENBROCK. The Jacobian of the corresponding nonlinear
C equations is five-diagonal. See Example 5.3 of "The Secant/Finite 
C Difference Algorithms for Solving Sparse Nonlinear Systems of Equations"
C by Guangye Li (SIAM Journal on Numerical Analysis, 1988).
            F = 0.0D0
            IF (N >= 4) THEN
            F = (4.0D0*(Y(1)-Y(2)**2) + Y(2)-Y(3)**2)**2
     1        + (8.0D0*Y(2)*(Y(2)**2-Y(1))- 2.0D0*(1.0D0-Y(2))
     1        + 4.0D0*(Y(2)-Y(3)**2) + Y(3)-Y(4)**2)**2
            DO I = 3, N-2
                F = F + (8.0D0*Y(I)*(Y(I)**2-Y(I-1))-2.0D0*(1.0D0-Y(I))
     1            + 4.0D0*(Y(I)-Y(I+1)**2)+Y(I-1)**2-Y(I-2)
     1            + Y(I+1)-Y(I+2)**2)**2
            END DO
            F = F + (8.0D0*Y(N-1)*(Y(N-1)**2-Y(N-2))
     1            - 2.0D0*(1.0D0-Y(N-1))+4.0D0*(Y(N-1)-Y(N)**2)
     1            + Y(N-2)**2-Y(N-3))**2
            F = F + (8.0D0*Y(N)*(Y(N)**2-Y(N-1))-2.0D0*(1.0D0-Y(N))
     1            + Y(N-1)**2-Y(N-2))**2
            END IF
        ELSE IF (TRIM(PROB) .EQ. 'FLETCBV2') THEN
C See CUTEr. 
            h = 1.0D0/DFLOAT(N+1)
            TEMP1 = Y(1)**2
            DO I = 1, N-1 
                TEMP1 = TEMP1 + (Y(I)-Y(I+1))**2
            END DO
            TEMP1 = TEMP1 + Y(N)**2 
            TEMP2 = 0.0D0
            DO I = 1, N
                TEMP2 = TEMP2 + 2.0D0*Y(I)+COS(Y(I))
            END DO
            F = 0.5D0*TEMP1 - h**2*TEMP2 -Y(N)

        ELSE IF (TRIM(PROB) .EQ. 'FLETCBV3') THEN
C See CUTEr. MODIFIED by Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." 
            p = 1.0D-8
            h = 1.0D0/DFLOAT(N+1)
            KAPPA = 1.0D0
            TEMP1 = Y(1)**2
            DO I = 1, N-1
                TEMP1 = TEMP1 + (Y(I)-Y(I+1))**2
            END DO
            TEMP1 = TEMP1 + Y(N)**2 
            TEMP2 = 0.0D0
            DO I = 1,N
                TEMP2 = TEMP2 +
     1           1.0D2*SIN(Y(I)*1.0D-2)*(1.0D0+2.0D0/h**2) +
     1           1.0D0/h**2*KAPPA*COS(Y(I)) 
            END DO
            F = 0.5D0*p*TEMP1 - p*TEMP2
            
        ELSE IF (TRIM(PROB) .EQ. 'FLETCHCR') THEN
C See CUTEr.
            F = 0.0D0
            DO I = 1, N-1
                F = F + (Y(I+1)-Y(I)+1.0D0-Y(I)**2)**2
            END DO
            F = 1.0D2*F
        ELSE IF (TRIM(PROB) .EQ. 'FREUROTH') THEN
            TEMP1 = 0.0D0
            TEMP2 = 0.0D0
            DO I = 1, N-1 
                TEMP1 = TEMP1 + ((5.0D0-Y(I+1))*Y(I+1)**2 + Y(I) -
     1           2.0D0*Y(I+1)-13.0D0)**2
                TEMP2 = TEMP2 +
     1          ((1.0D0+Y(I+1))*Y(I+1)**2+Y(I)-14.0D0*Y(I+1)-29.0D0)**2
            END DO
            F = TEMP1 + TEMP2
        ELSE IF (TRIM(PROB) .EQ. 'GENBROWN') THEN
C Generalized Brown Function 1. See Ying-Jie Li, Dong-Hui Li, Truncated
C regularized Newton Method for Convex Minimization, Problem 6, 
C and Luksan, Vicek, Sparse and partially separable test problems for 
C unconstrained and equality constrained optimization.
            F = 0.0D0
            DO I = 2, N 
                F = F+(Y(I-1)-3.0D0)**2+(Y(I-1)-Y(I))**2
     1            + EXP(20.0D0*(Y(I-1)-Y(I)))
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'GENHUMPS') THEN
C See CUTEr.
            ZETA = 2.0D0
            F = 0.0D0
            DO I = 1, N-1
                F = F + SIN(ZETA*Y(I))**2*SIN(ZETA*Y(I+1))**2
     1            + 0.05D0*(Y(I)**2+Y(I+1)**2)
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'GENROSE') THEN
C See CUTEr. 
C Compare with CHROSEN and ROSENBROCK.
C The starting point is (I/(N+1)).
            F = 1.0D0
            DO I = 2, N
                F = F + 1.0D2*(Y(I)-Y(I-1)**2)**2+(Y(I)-1.0D0)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'INDEF') THEN
C See CUTEr. MODIFIED by Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." 
            TEMP1 = 0.0D0
            DO I = 1, N
                TEMP1 = TEMP1 + 1.0D2*SIN(1.0D-2*Y(I))
            END DO
            TEMP2 = 0.0D0
            DO I = 2, N-1
                 TEMP2 = TEMP2 + COS(2.0D0*Y(I)-Y(N)-Y(1)) 
            END DO
            F = TEMP1 + 0.5D0*TEMP2
        ELSE IF (TRIM(PROB) .EQ. 'INTEGREQ') THEN
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 29 and also CUTEr.
            h = 1.0D0/DFLOAT(N+1)
            DO I =1, N
                T(I) = DFLOAT(I)/DFLOAT(N+1)
            END DO
            F = 0.0D0
            DO I = 1, N
                TEMP1 = 0.0D0
                DO J = 1, I
                    TEMP1 = TEMP1 + T(J)*(Y(J) + T(J) + 1.0D0)**3
                END DO
                TEMP2 = 0.0D0
                DO J = I+1, N
                TEMP2 = TEMP2 + (1.0D0-T(J))*(Y(J)+T(J)+1.0D0)**3
                END DO
                F = F + (Y(I)+
     1          0.5D0*h*((1-T(I))*TEMP1+T(I)*TEMP2))**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'LIARWHD') THEN
            F = 0.0D0
            DO I =1 ,N
                F = F + 4.0D0*(Y(I)**2-Y(1))**2+(Y(I)-1.0D0)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'LILIFUN3') THEN
C See Ying-Jie Li, Dong-Hui Li, Truncated
C regularized Newton Method for Convex Minimization, Problem 3.
            F = 0.0D0
            DO I = 2, N
                F = F + EXP(Y(I)-Y(I-1))**2 + (Y(I)-Y(I-1))**2
     1            + 2.0D0*Y(I)**4+4.0D0*Y(I-1)**4
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'LILIFUN4') THEN
C See Ying-Jie Li, Dong-Hui Li, Truncated
C regularized Newton Method for Convex Minimization, Problem 4.
            F = 0.0D0
            DO I = 2, N
                F = F +
     1           0.5D0*(Y(I)-Y(I-1))**2
     1           +SIN(Y(I)-Y(I-1))
     1           +2.0D0*(2.0D0*Y(I)+3.0D0*Y(I-1)-15.0D0)**4
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'MOREBV' .OR. 
     1   TRIM(PROB) .EQ. 'MOREBVL') THEN
C See CUTEr. MOREBVL uses the start point suggested by Luksan, Matonoha,
C Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." 
            h = 1.0D0/DFLOAT(N+1)
            F = 0.0D0
            DO I = 1, N
                F = F + (2.0D0*Y(I)-Y(I-1)-Y(I+1)+
     1           0.5D0*h*h*(Y(I)+DFLOAT(I)*h+1.0D0)**3)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'NCB20') THEN
C See CUTEr. Corrected by Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." 
            F = 0.0D0
            IF (N >= 20) THEN
                TEMP1 =0.0D0
                DO I = 1, N-30
                    TEMP = 0.0D0
                    DO J = 1, 20
                        TEMP = TEMP + Y(I+J-1)/(1+Y(I+J-1)**2)
                    END DO
                    TEMP1 = TEMP1 + 10.0D0/DFLOAT(I)*TEMP**2
                    TEMP = 0.0D0
                    DO J = 1, 20
                        TEMP = TEMP + Y(I+J-1)
                    END DO
                    TEMP1 = TEMP1 - 0.2D0*TEMP
                END DO
                TEMP2 = 0.0D0
                DO I = 1, N-10
                    TEMP2 = TEMP2 + Y(I)**4+2.0D0 
                END DO
                TEMP3 = 0.0D0
                DO I = 1, 10
                    TEMP3 = TEMP3 + 
     1               Y(I)*Y(I+10)*Y(I+N-10) +2.0D0*Y(I+N-10)**2
                END DO
                F = 2.0D0 + TEMP1 + TEMP2 + 1.0D-4*TEMP3
            END IF
        ELSE IF (TRIM(PROB) .EQ. 'NCB20B') THEN
C See CUTEr. Corrected by Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." 
            F = 0.0D0
            IF (N >= 20) THEN
                TEMP1 = 0.0D0
                DO I = 1, N-19
                    TEMP = 0.0D0
                    DO J = 1, 20
                        TEMP = TEMP + Y(I+J-1)/(1+Y(I+J-1)**2)
                    END DO
                    TEMP1 = TEMP1 + 10.0D0/DFLOAT(I)*TEMP**2
                    TEMP = 0.0D0
                    DO J = 1, 20
                        TEMP = TEMP + Y(I+J-1)
                    END DO
                    TEMP1 = TEMP1 - 0.2D0*TEMP
                END DO 
                TEMP2 = 0.0D0
                DO I = 1, N
                    TEMP2 = TEMP2 + (1.0D2*Y(I)**4+2.0D0)
                END DO
                F = TEMP1 +TEMP2
            END IF
        ELSE IF(TRIM(PROB) .EQ. 'NONCVXU2') THEN
C See CUTEr.
            F = 0.0D0
            DO I = 1, N
                J = MOD(3*I-2,N)
                K = MOD(7*I-3,N)
                TEMP = Y(I) + Y(J+1) + Y(K+1)
                F = F + TEMP**2 + 4.0D0*COS(TEMP) 
            END DO
            
        ELSE IF(TRIM(PROB) .EQ. 'NONCVXUN') THEN
C See CUTEr.
            F = 0.0D0
            DO I = 1, N
                J = MOD(2*I-1,N)
                K = MOD(3*I-1,N)
                TEMP = Y(I) + Y(J+1) + Y(K+1)
                F = F + TEMP**2 + 4.0D0*COS(TEMP) 
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'NONDIA') THEN
C See CUTEr.
            F = (Y(1)-1.0D0)**2
            DO I = 2, N
                F = F + (1.0D2*Y(1)-Y(I-1)**2)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'NONDQUAR') THEN
            F = (Y(1)-Y(2))**2 + (Y(N-1)-Y(N))**2
            DO I = 1, N-2
                F = F + (Y(I) + Y(I+1) + Y(N))**4
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'PENALTY1') THEN
            DP1 = 0.0D0
            DP2 = 0.0D0
            DO I = 1, N
                DP1 = DP1 + (Y(I)-1.0D0)**2
                DP2 = DP2 + Y(I)**2
            END DO
            F=1.0D-5*DP1 +(0.25D0-DP2)**2
        ELSE IF (TRIM(PROB) .EQ. 'PENALTY2') THEN
            F = 0.0D0
            temp = 0.0D0
            do i = 2,n
                F = F + 
     1           (exp(Y(I-1)*0.1D0)+exp(Y(I)*0.1D0)-
     1            exp(DFLOAT(i-1)*0.1D0)
     1           -exp(DFLOAT(i)*0.1D0))**2 + (exp(Y(I)*0.1D0)
     1           -exp(-0.1D0))**2
            end do
            do i =1,n
                temp = temp + DFLOAT(n-i+1)*Y(I)*Y(I)
            end do
            F = F + (1.0D0-temp)*(1.0D0-temp) + (Y(1)-0.2D0)*
     1         (Y(1)-0.2D0)
        ELSE IF (TRIM(PROB) .EQ. 'PENALTY3' .OR.
     1   TRIM(PROB) .EQ. 'PENALTY3P') THEN
            F = 0.0D0
            R = 0.0D0
            S = 0.0D0
            PD = 0.0D0
            DO I = 1, N
                PD = PD + (Y(I)**2-DFLOAT(N))
            END DO
            IF (TRIM(PROB) .EQ. 'PENALTY3') THEN
                ALPHA = 1.0D0 
            ELSE
                ALPHA = 1.0D-3  ! Suggested by Powell, "The NEWUOA Software
C            for unconstrained optimization without derivatives", 2004,
C            section 8.
            END IF
            do i = 1,n-2
                R = R + (Y(I) + 2.0D0* Y(I+1) + 
     1            10.0D0*Y(I+2) - 1.0D0)**2
                S = S + (2.0D0*Y(I) + Y(I+1) - 3.0D0)**2
            end do
            do i = 1, n/2
                F = F + (Y(I) - 1.0D0)**2
            end do
            F = F + PD**2 + 
     1       ALPHA*(1.0D0+R*exp(Y(N)) + S*exp(Y(N-1)) + R*S)
        ELSE IF (TRIM(PROB) .EQ. 'POWELLSG') THEN
            F = 0.0D0
            DO I = 1, N/4
               J = 4*(I-1) + 1
               F = F + (Y(J)+1.0D1*Y(J+1))**2 +5.0D0*(Y(J+2)-Y(J+3))**2
     1               + (Y(J+1)-2.0D0*Y(J+2))**4+1.0D1*(Y(J)-Y(J+3))**4
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'POWER') THEN
C See CUTEr.
            F = 0.0D0
            DO I = 1, N
                F = F + (DFLOAT(I)*Y(I))**2
C                F = F + (Y(I))**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'ROSENBROCK') THEN
C In claasical Rosenbrock function, alpha = 100.0D0
C When alpha = 4.0D0, this function is essentially the same as Chrosen,
C except the order of the variables.
            alpha = 100.0D0
            F = 0.0D0
            DO I = 1, N-1
                F = F + (1.0D0-Y(I))*(1.0D0-Y(I)) + 
     1       alpha*(Y(I+1)-Y(I)*Y(I))*(Y(I+1)
     1       -Y(I)*Y(I))
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'SBRYBND' 
     1 .OR. TRIM(PROB) .EQ. 'SBRYBNDL') THEN
C See CUTEr. And see Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." The scaling parameter of Luksan is 6
C instead of 12. There is a typo in Luksan: the square is missed, which
C makes the function unbounded from below.
            ML = 5
            MU = 1
            IF (TRIM(PROB) .EQ. 'SBRYBND') THEN
                SCALE = 12.0D0
            ELSE
                SCALE = 6.0D0
            END IF
            F = 0.0D0
            IF (N>=2) THEN
                DO I = 1, N
                    YTEMP(I) =EXP(SCALE*(DFLOAT(I-1)/DFLOAT(N-1)))*Y(I)
                END DO
                DO I = 1, N
                    TEMP = 0.0D0
                    DO J = MAX(1, I-ML), MIN(N,I+MU) 
                        IF (J/=I) THEN 
                            TEMP = TEMP+YTEMP(J)*(1.0D0+YTEMP(J))
                        END IF
                    END DO
                   F = F + (YTEMP(I)*(2.0D0+5.0D00*YTEMP(I)*YTEMP(I))
     1                  +1.0D0-TEMP)**2
                END DO
            END IF
        ELSE IF (TRIM(PROB) .EQ. 'SCHMVETT') THEN
C See CUTEr.
            PI = ACOS(-1.0D0)
            F = 0.0D0
            DO I = 1, N -2 
                F = F - 1.0D0/(1.0D0+(Y(I)-Y(I+1))**2)
     1               - SIN(0.5D0*(PI*Y(I+1)+Y(I+2)))
     1               - EXP(-((Y(I)+Y(I+2))/Y(I+1)-2.0D0)**2)
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'SCOSINE' .OR.
     1   TRIM(PROB) .EQ. 'SCOSINEL') THEN
C See CUTEr. And see Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." The scaling parameter of Luksan is 6
C instead of 12. 
            IF (TRIM(PROB) .EQ. 'SCOSINE') THEN
                SCALE = 12.0D0
            ELSE
                SCALE = 6.0D0
            END IF
            F = 0.0D0
            IF (N >= 2) THEN
                DO I = 1, N
                    YTEMP(I) =EXP(SCALE*(DFLOAT(I-1)/DFLOAT(N-1)))*Y(I)
                END DO
                DO I = 1, N-1
                    F = F +COS(YTEMP(I)**2-0.5D0*YTEMP(I+1)) 
                END DO
            END IF
        ELSE IF (TRIM(PROB) .EQ. 'SEROSE') THEN
C Seven-Diagonal ROSENBROCK. The Jacobian of the corresponding nonlinear
C equations is 7-diagonal. See Example 5.4 of "The Secant/Finite 
C Difference Algorithms for Solving Sparse Nonlinear Systems of Equations"
C by Guangye Li (SIAM Journal on Numerical Analysis, 1988).
            F = 0.0D0
            IF (N >= 6) THEN
            F = (4.0D0*(Y(1)-Y(2)**2) + Y(2)-Y(3)**2 + Y(3)-Y(4)**2)**2
            DO I = 2, N-1
                F = F + (8.0D0*Y(I)*(Y(I)**2-Y(I-1))-2.0D0*(1.0D0-Y(I))
     1            + 4.0D0*(Y(I)-Y(I+1)**2)
     1            + Y(I-1)**2-Y(I-2)+Y(I+1)-Y(I+2)**2
     1            + Y(I-2)**2-Y(I-3)+Y(I+2)-Y(I+3)**2)**2
            END DO
            F = F + (8.0D0*Y(N)*(Y(N)**2-Y(N-1))-2.0D0*(1.0D0-Y(N))
     1            + Y(N-1)**2-Y(N-2)+Y(N-2)**2-Y(N-3))**2
            END IF
        ELSE IF(TRIM(PROB) .EQ. 'SINQUAD') THEN
C See CUTEr.
            F = (Y(1)-1.0D0)**4
            DO I = 2, N-1
                F = F + (SIN(Y(I)-Y(N))-Y(1)**2+Y(I)**2)**2
            END DO
            F = F +(Y(N)**2-Y(1)**2)**2
    
        ELSE IF (TRIM(PROB) .EQ. 'SPHRPTS') THEN
            F = 0.0D0
            IF (MOD(N,2) /= 0) THEN
                WRITE(*,*) "ERROR in CALFUN: Dimension for SPHRPTS
     1should be even."
            ELSE
            do i = 2, n/2
                do j = 1, i-1
                    F = F + 
     1              1.0D0/((cos(Y(2*I-1))*cos(Y(2*I))
     1              -cos(Y(2*J-1))*cos(Y(2*J)))**2
     1              +(sin(Y(2*I-1))*cos(Y(2*I))
     1              -sin(Y(2*J-1))*cos(Y(2*J)))**2
     1              +(sin(Y(2*I))-sin(Y(2*J)))**2)
                end do
            end do
            END IF
        ELSE IF (TRIM(PROB) .EQ. 'SPARSINE') THEN
            F = 0.0D0
            DO I = 1, N
                J1 = MOD(2*I-1,N)+1
                J2 = MOD(3*I-1,N)+1
                J3 = MOD(5*I-1,N)+1
                J4 = MOD(7*I-1,N)+1
                J5 = MOD(11*I-1,N)+1
                F = F + DFLOAT(I)*(SIN(Y(I))+SIN(Y(J1))+SIN(Y(J2))
     1            + SIN(Y(J3))+SIN(Y(J4))+SIN(Y(J5)))**2
            END DO
            F = 0.5D0*F
        ELSE IF (TRIM(PROB) .EQ. 'SPARSQUR') THEN
            F = 0.0D0
            DO I = 1, N
                J1 = MOD(2*I-1,N)+1
                J2 = MOD(3*I-1,N)+1
                J3 = MOD(5*I-1,N)+1
                J4 = MOD(7*I-1,N)+1
                J5 = MOD(11*I-1,N)+1
                F = F + DFLOAT(I)*(Y(I)**2+Y(J1)**2+Y(J2)**2
     1            + Y(J3)**2+Y(J4)**2+Y(J5)**2)**2
            END DO
            F = 0.125D0*F
        ELSE IF (TRIM(PROB) .EQ. 'SPMSRTLS') THEN
            DO I = 1, N
                PAR(I) = SIN(DFLOAT(I)**2)
            END DO
            M = (N+2)/3
            F = 0.0D0
            DO I = 1, M
                F1 = 0.0D0
                F2 = 0.0D0
                F3 = 0.0D0
                F4 = 0.0D0
                F5 = 0.0D0
                J = 3*(I-1) + 1
                IF (I>2) THEN
                    F1 = (Y(J-4)*Y(J-1)-PAR(J-4)*PAR(J-1))**2
                END IF
                F3 = (Y(J)**2-PAR(J)**2)**2
                IF (I>1) THEN
                    F2 = (Y(J-3)*Y(J-1)+Y(J-1)*Y(J)
     1               -PAR(J-3)*PAR(J-1)-PAR(J-1)*PAR(J))**2
                    F3 = F3 + (Y(J-2)*Y(J-1)-PAR(J-2)*PAR(J-1))**2
                END IF
                IF (I<M) THEN
                    F3 = F3 + (Y(J+2)*Y(J+1)-PAR(J+2)*PAR(J+1))**2
                    F4 = (Y(J+3)*Y(J+1)+Y(J+1)*Y(J)
     1               -PAR(J+3)*PAR(J+1)-PAR(J+1)*PAR(J))**2
                END IF
                IF (I <M-1) THEN
                    F5 = (Y(J+4)*Y(J+1)-PAR(J+4)*PAR(J+1))**2
                END IF
                F = F + F1 + F2 + F3 + F4 + F5 
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'SROSENBR') THEN
            F = 0.0D0
            DO I = 1, N/2
                F = F + 1.0D2*(Y(2*I)-Y(2*I-1)**2)**2 
     1           + (Y(2*I-1)-1.0D0)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'TOINTGSS') THEN
            F = 0.0D0
            DO I = 1, N-2
                F = F + (1.0D1/DFLOAT(N+2)+Y(I+2)**2)*(2.0D0-
     1           EXP(-(Y(I)-Y(I+1))**2/(1.0D-1+Y(I+2)**2)))
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'TOINTTRIG') THEN
C See Ph. L. Toint, "Some Numerical Results Using a Sparse Matrix
C Updating Formula in Unconstrained Optimization", Example 3.6. 
            F = 0.0D0
            DO I = 2, N
                BETAI = 1.0D0 + 1.0D-1*DFLOAT(I)
                DO J = 1, I-1
                    ALPHA = DFLOAT(5*(1+MOD(I,5)+MOD(J,5)))
                    BETAJ = 1.0D0 + 1.0D-1*DFLOAT(J)
                    C = DFLOAT(I+J)*1.0D-1
                    F = F + ALPHA*SIN(BETAI*Y(I)+BETAJ*Y(J)+C) 
                END DO
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'TQUARTIC') THEN
            F = (Y(1)-1.0D0)**2
            DO I = 1, N-1
                F = F + (Y(1)**2-Y(I+1)**2)**2
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'TRIROSE') THEN
C Tri-Diagonal ROSENBROCK. The Jacobian of the corresponding nonlinear
C equations is 3-diagonal. See Example 5.2 of "The Secant/Finite 
C Difference Algorithms for Solving Sparse Nonlinear Systems of Equations"
C by Guangye Li (SIAM Journal on Numerical Analysis, 1988).
            F = 0.0D0
            IF (N >= 2) THEN
            F = 16.0D0*(Y(1)-Y(2)**2)**2
            DO I = 2, N-1
                F = F + (8.0D0*Y(I)*(Y(I)**2-Y(I-1))
     1           -2.0D0*(1.0D0-Y(I))
     1           +4.0D0*(Y(I)-Y(I+1)**2))**2
            END DO
            F = F + (8.0D0*Y(N)*(Y(N)**2-Y(N-1))
     1           -2.0D0*(1.0D0-Y(N)))**2
            END IF
        ELSE IF (TRIM(PROB) .EQ. 'VARDIM') THEN
            F = 0.0D0
            do i = 1, N
                F = F + DFLOAT(i)*(Y(I)-1.0D0)
            end do
            F = F*F
            F = F + F*F
            do i = 1, N
                F = F + (Y(I) - 1.0D0)*(Y(I) - 1.0D0)
            end do
        ELSE IF (TRIM(PROB) .EQ. 'WOODS') THEN
            F = 0.0D0
            DO I = 1, N/4
                J = 4*I
                F = F + 1.0D2*(Y(J-2)-Y(J-3)**2)**2+(1.0D0-Y(J-3))**2
     1           + 9.0D1*(Y(J)-Y(J-1)**2)**2 + (1.0D0-Y(J-1))**2 +
     1           1.0D1*(Y(J-2)+Y(J)-2.0D0)**2+1.0D-1*(Y(J-2)-Y(J))**2
            END DO
        ELSE
            WRITE(*,*) "ERROR in CALFUN: Unknown problem name."
            STOP
        END IF 
        
        NFEVL = NFEVL + 1
        FOPTIMAL = DMIN1(F,FOPTIMAL)
        FREC(NFEVL) = F
        FOPTREC(NFEVL) = FOPTIMAL

        RETURN
      END
