      SUBROUTINE HSSD(N, NPT, HD, HQ, PQ, XPT, D)

      IMPLICIT NONE
      INTEGER(kind = 4), intent(in) :: N, NPT
      REAL(kind = 8), intent(out) :: HD(N)
      REAL(kind = 8), intent(in) :: HQ(N*(N+1)/2), PQ(NPT), XPT(NPT,N)
      REAL(kind = 8), intent(in) :: D(N)

      INTEGER(kind = 4) I, J, K, IH
      REAL(kind = 8) :: TEMP,ZERO
      ZERO = 0.0D0

      HD=ZERO
      DO  K=1,NPT
          TEMP=ZERO
        DO  J=1,N
          TEMP=TEMP+XPT(K,J)*D(J)
        END DO
        TEMP=TEMP*PQ(K)
        DO  I=1,N
          HD(I)=HD(I)+TEMP*XPT(K,I)
        END DO
      END DO
      IH=0
      DO  J=1,N
        DO  I=1,J
            IH=IH+1
            IF (I .LT. J) HD(J)=HD(J)+HQ(IH)*D(I)
            HD(I)=HD(I)+HQ(IH)*D(J)
        END DO
      END DO
      END

      SUBROUTINE ONEDIMSEARCH(X, D, F0, F1, N, NF, NFTEST, LONG, IPRINT)
        IMPLICIT NONE
        INTEGER(KIND = 4), INTENT(IN) :: N, NFTEST, IPRINT
        INTEGER(KIND = 4), INTENT(INOUT) :: NF, LONG
        REAL(KIND = 8), INTENT(INOUT) :: D(N), X(N), F1, F0
C Search from X along D for a better point. The function value at 
C X-D is F0, and the one at X is F1 (F1<F0). N is the dimension of the 
C problem. The objective function is calculated by calling CALFUN, and
C NF is increased by 1 every time the objective function calculated.
C When the subroutine exits, X should be the point found, and F1 should 
C be the corresponding function value.
C Termination:
C Suppose the iteratives generated are X_0 (=X-D), X_1 (=X), X_2, X_3, ...,
C and the corresponding function values are F_0 (=F0), F_1 (=F1), F_2, F_3, ...
C Then the subroutine will terminate if one of the following holds:
C 1. ||X_K - X_(K-1)|| < 1/3*||D||
C 2. NF > NFTEST;
C 3. F_K - F_(K-1) > 1/3*(F1-F0);
C 4. There is evidence that X_K is nearly optimal (e.g., X_K is the global
C    minimizer of the trustregion subproblem, and ARED/PRED <= 1.5).

        REAL(KIND = 8) :: RED  ! RED = F1 - F0
        REAL(KIND = 8) :: XNEW(N), F, ARED, STEP
        REAL(KIND = 8) :: PRED, RATIO, DELTA, ALPHA, BETA
        REAL(KIND = 8) :: T0, T1, T
        INTEGER(KIND = 4) :: I

C During the iterations, XNEW is the new trial point; 
C F=F_K, F1=F_(K-1), F0=F(K-2); 
C ARED = F_K-F_(K-1);
C X_K = X + T_K*D (thus T_0 = -1, T_1 = 0);
C STEP = T_K-T_(K-1);
C PRED is the reduction predicted by model, RATIO = ARED/PRED,
C and DELTA is the trust-region radius (for variable T).


        RED = F1 - F0
        T0 = -1
        T1 = 0
        T = 0
        DELTA = 1.0D0
        PRED = DELTA * RED
        STEP = DELTA
        LONG = 0
        DO WHILE (1 == 1) 
            WRITE(*,*) "TRY LONG STEP."
            XNEW = X + STEP * D
            NF = NF + 1
            IF (NF >= NFTEST) THEN
                EXIT
            END IF

            CALL CALFUN (N,XNEW,F)
            IF (IPRINT .EQ. 3) THEN
                PRINT 100, NF,F,(XNEW(I),I=1,N)
  100           FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1              '    The corresponding X is:'/(2X,5D15.6))
            END IF

            IF (F < F1) THEN
                WRITE(*,*) "SUCCESSFUL LONG STEP."
                LONG = LONG + 1
                X = XNEW
                T = T + STEP
            END IF
            ARED = F - F1
C            IF (ARED >= DFLOAT(1)/DFLOAT(3) * RED) THEN
            IF (ARED >= 0.3D0 * RED) THEN
                EXIT
            END IF
            
C            ARED = F - F1
C            IF (ARED <= DFLOAT(1)/DFLOAT(3) * RED) THEN
C                WRITE(*,*) "SUCCESSFUL LONG STEP."
C                LONG = LONG + 1
C                X = XNEW
C                T = T + STEP
C            ELSE
C                F = F1
C                EXIT
C            END IF
            
            RATIO = ARED/PRED

            IF (DABS(STEP) < DELTA .AND. RATIO <= 1.5D0) THEN
                EXIT
            END IF

C            IF (RATIO >= 0.6D0) THEN
                DELTA = 1.1D0*DELTA
C            END IF
            
C DOES SCALING BY RED HELP TO IMPROVE ACCURACY???????
C            ALPHA = 2.0D0*(((F0-F)/RED)*(T1-T)-((F1-F)/RED)*(T0-T))/
C     1              ((T0-T)*(T1-T)*(T0-T1))*RED
C            BETA = (((F0-F)/RED)*(T1-T)**2-((F1-F)/RED)*(T0-T)**2)/
C     1              ((T0-T)*(T1-T)*(T1-T0))*RED
C            STEP =0.5D0*(((F0-F)/RED)*(T1-T)**2-((F1-F)/RED)*(T0-T)**2)/
C     1              ((F0-F)/RED*(T1-T)-(F1-F)/RED*(T0-T))
            ALPHA = 2.0D0*((F0-F)*(T1-T)-(F1-F)*(T0-T))/
     1              ((T0-T)*(T1-T)*(T0-T1))
            BETA = ((F0-F)*(T1-T)**2-(F1-F)*(T0-T)**2)/
     1              ((T0-T)*(T1-T)*(T1-T0))

            IF (ALPHA <= 0) THEN
                IF (BETA >= 0) THEN
                    WRITE(*,*) "BETA ERROR!"
                    WRITE(*,*) PRED/RED
                    WRITE(*,*) "ALPHA = ",ALPHA
                    WRITE(*,*) "BETA =",BETA
                    WRITE(*,*) "STEP=", STEP
                    WRITE(*,*) "T0 = ", T0
                    WRITE(*,*) "F0 = ", F0
                    WRITE(*,*) "T1 = ", T1
                    WRITE(*,*) "F1 = ", F1
                    WRITE(*,*) "T = ", T
                    WRITE(*,*) "F = ", F
                    EXIT
                END IF
                STEP = DELTA
            ELSE
                STEP =0.5D0*((F0-F)*(T1-T)**2-(F1-F)*(T0-T)**2)/
     1              ((F0-F)*(T1-T)-(F1-F)*(T0-T))
                STEP = DMAX1(-DELTA, DMIN1(DELTA, STEP))
            END IF

            IF (DABS(STEP) < DFLOAT(1)/DFLOAT(3)) THEN
C                WRITE(*,*) "SHORT STEP IN LONGSTEP."
                EXIT
            END IF

            PRED = (0.5D0*ALPHA*STEP*STEP + BETA*STEP)
            IF (PRED >= 1.0D-5*RED) THEN
                WRITE(*,*) "ERROR: FAIL TO REDUCE MODEL IN LONGSTEP."
                WRITE(*,*) PRED/RED
                WRITE(*,*) "ALPHA = ",ALPHA
                WRITE(*,*) "BETA =",BETA
                WRITE(*,*) "STEP=", STEP
                WRITE(*,*) "T0 = ", T0
                WRITE(*,*) "F0 = ", F0
                WRITE(*,*) "T1 = ", T1
                WRITE(*,*) "F1 = ", F1
                WRITE(*,*) "T = ", T
                WRITE(*,*) "F = ", F
                EXIT
            END IF


            F0 = F1
            F1 = F
            T0 = T1
            T1 = T

        END DO
            D = D*(T+1.0D0)
            F1 = DMIN1(F, F1)

      END



      SUBROUTINE CHECKINT(XPT,XOPT,XBASE,HQ,PQ,GQ,FOPT,
     1           FVAL,N,NPT,KOPT,ERROR)
        IMPLICIT NONE
        INTEGER(KIND = 4), INTENT(IN) :: N, NPT, KOPT
        REAL(KIND = 8), INTENT(IN) :: XPT(NPT,N), XOPT(N), XBASE(N),
     1   HQ(N*(N+1)/2),PQ(NPT),GQ(N),FOPT,FVAL(NPT)
        REAL(KIND = 8), INTENT(OUT) :: ERROR(NPT)

        REAL(KIND = 8) :: q(NPT), FTEST, HALF,ZERO
        INTEGER(KIND = 4) :: I, J, K, L, IH
        HALF = 0.5D0
        ZERO = 0.0D0
        ERROR = 0.0D0
        

        DO I = 1, NPT
            CALL CALFUN(N, XPT(I,:)+XBASE, FTEST)
            IF(DABS(FTEST-FVAL(I))/DABS(FTEST)>1D-8.OR.
     1        (FTEST<FOPT .AND. I /= KOPT)) THEN
                ERROR(I) = 999999999.0D0
                WRITE(*,*) "INTERPOLTION ERROR!! FVAL(K) /=F(XPT(K))!!"
                WRITE(*,*) (FTEST-FVAL(I))/DABS(FTEST)
                WRITE(*,*) FTEST, FVAL(I), FOPT
                return
            END IF
        END DO
        CALL CALFUN(N, XOPT+XBASE, FTEST)
        IF (DABS(FTEST-FOPT)/DABS(FTEST)>1D-8) THEN
            ERROR = 999999999.0D0
            WRITE(*,*) "INTERPOLTION ERROR!! F(XOPT) /= FOPT!!"
            WRITE(*,*) (FTEST-FOPT)/DABS(FTEST)
            WRITE(*,*) FTEST, FOPT
            return
        END IF
        IF(dot_product(XOPT-XPT(KOPT,:),XOPT-XPT(KOPT,:)) 
     1      /dot_product(XOPT,XOPT) > 1D-16) THEN
            ERROR = -999999999.0D0
            WRITE(*,*) "INTERPOLTION ERROR!! XOPT /= XPT(KOPT)!!"
            WRITE(*,*) dot_product(XOPT-XPT(KOPT,:),XOPT-XPT(KOPT,:))
     1        /dot_product(XOPT,XOPT)
            WRITE(*,*) XOPT, XPT(KOPT, :)
            return
        END IF

        
        q = 0.0D0
        do i = 1, NPT
          IH = 1
          do k = 1, N
              do l = 1, k 
                  if (k == l) then
                     q(i) = q(i) + HALF * XPT(i, k) * HQ(IH) * XPT(i, l)
                  else
                     q(i) = q(i) + XPT(i, k) * HQ(IH) * XPT(i, l)
                  end if
                  IH = IH +1
              end do
          end do
          do j = 1, NPT
            q(i) = q(i) +HALF* PQ(j)*(dot_product(XPT(i,:),XPT(j,:))**2)
          end do
          q(i) = q(i) + dot_product(XPT(i,:), GQ)
          ERROR(i) = FVAL(I)-Q(I)
        end do
        ERROR = ERROR - (FVAL(KOPT) - Q(KOPT))
C        WRITE(*,*) "ERROR = ", ERROR
        DO I = 1, NPT
          ERROR(I) = ERROR(I)/DMAX1(ABS(FVAL(I)-FVAL(KOPT)),1D-16)
        END DO
        IF (MAXVAL(DABS(ERROR))>1.0D-5) THEN
          WRITE(*,*) "INTERPOLTION ERROR!!!!"
          WRITE(*,*) "ERROR = ", MAXVAL(DABS(ERROR))
        END IF
      END

      SUBROUTINE HESSIAN(H, HQ, PQ, XPT, N, NPT)
        IMPLICIT NONE
        INTEGER(KIND = 4), INTENT(IN) :: N, NPT
        REAL(KIND = 8), INTENT(IN) :: XPT(NPT,N), HQ(N*(N+1)/2),PQ(NPT)
        REAL(KIND = 8), INTENT(OUT) :: H(N*(N+1)/2)
        REAL(KIND = 8) :: TEMP 
        INTEGER(KIND = 4) :: I, J, K, IH
        H = HQ
        DO K = 1, NPT
            IH=0
            DO I=1,N
                TEMP=PQ(K)*XPT(K,I)
                DO J=1,I
                    IH=IH+1
                    H(IH)=H(IH)+TEMP*XPT(K,J)
                END DO
            END DO
        END DO
      END


      SUBROUTINE HESSIANSQ(HSQ, HQ, PQ, XPT, N, NPT)
        IMPLICIT NONE
        INTEGER(KIND = 4), INTENT(IN) :: N, NPT
        REAL(KIND = 8), INTENT(IN) :: XPT(NPT,N), HQ(N*(N+1)/2),PQ(NPT)
        REAL(KIND = 8), INTENT(OUT) :: HSQ
        REAL(KIND = 8) :: TEMP, H(N*(N+1)/2)
        INTEGER(KIND = 4) :: I, J, K, IH
        H = HQ
        DO K = 1, NPT
            IH=0
            DO I=1,N
                TEMP=PQ(K)*XPT(K,I)
                DO J=1,I
                    IH=IH+1
                    H(IH)=H(IH)+TEMP*XPT(K,J)
                END DO
            END DO
        END DO
        HSQ = 0.0D0
        IH = 0
        DO I=1,N
            DO J=1,I
                IH=IH+1
                TEMP = H(IH)*H(IH)
                IF (J < I) TEMP = TEMP + TEMP
                HSQ = HSQ + TEMP
            END DO
        END DO
      END



!      SUBROUTINE EIGVALUE(A,N,E)
!C E = EIG(A), in ascending order, where A is a symmetric matrix of 
!C order N. Only the upper triangular part of A will be used. 
!C A will not be disdoyed.
!         IMPLICIT NONE
!         INTEGER(KIND=4), INTENT(IN) :: N
!         REAL(KIND=8), INTENT(IN) :: A(N,N)
!         REAL(KIND=8), INTENT(OUT) :: E(N)
!
!         REAL(KIND=8) :: TA(N,N), WORK(1+6*N+2*N*N)
!         INTEGER(KIND=4) :: IWORK(3+5*N), INFO,LDA,LWORK,LIWORK
!         CHARACTER :: JOBZ='N'
!         CHARACTER :: UPLO='U'
!         LDA=N
!         LWORK=1+6*N+2*N*N
!         LIWORK=3+5*N
!         TA = A
!         CALL DSYEVD(JOBZ,UPLO, N, TA, LDA, E,WORK, LWORK, IWORK, 
!     1   LIWORK, INFO)
!         IF (INFO /= 0) WRITE(*,*) "ERROR IN EIGVALUE. INFO=", INFO
!      END
!
!
!      SUBROUTINE EIG(A, N, K, LAMBDA, V)
!C LAMBDA are the largest K eigenvalues of A, in descending order,
!C and the columns of V are corresponding eigenvectors, where A 
!C is a symmetric matrix of order N. Only the upper triangular 
!C part of A will be used. A will not be disdoyed.
!         IMPLICIT NONE
!         INTEGER(KIND=4), INTENT(IN) :: N, K
!         REAL(KIND=8), INTENT(IN) :: A(N,N) 
!         REAL(KIND=8), INTENT(OUT) :: LAMBDA(K), V(N,K) 
!
!         REAL(KIND=8) :: WORK(1+6*N+2*N**2),TEMPA(N,N), W(N)
!         INTEGER(KIND=4) :: IWORK(3+5*N) 
!         INTEGER(KIND=4) :: LWORK, LIWORK, INFO, I
!         CHARACTER :: JOBZ, UPLO
!         JOBZ = 'V'
!         UPLO = 'U'
!         LWORK = 1+6*N+2*N**2
!         LIWORK = 3+5*N
!         TEMPA = A
!         CALL DSYEVD( JOBZ, UPLO, N, TEMPA, N, W, WORK, LWORK, IWORK,
!     1   LIWORK, INFO )
!         IF (INFO /= 0) WRITE(*,*) "ERROR IN EIG. INFO=", INFO
!         DO I = 1, K
!            LAMBDA(I) = W(N-I+1)
!            V(:,I) = TEMPA(:, N-I+1)
!         END DO
!      END
!      
!c      SUBROUTINE DNOISE(X, N, F)
!c        INTEGER*4 N, I
!c        REAL*8 X(N)
!c        REAL*8 F, N1, N2, NI
!c        
!c
!c        N1 = 0.0D0
!c        N2 = 0.0D0
!c        NI = 0.0D0
!c        DO I = 1, N
!c            N1 = N1 + DABS(X(I))
!c            N2 = N2 + X(I)*X(I)
!c            NI = DMAX1(NI, DABS(X(I)))
!c        END DO
!c        N2 = SQRT(N2)
!c        F = 0.9D0*SIN(1.0D2*N1)*COS(1.0D2*NI)+0.1D0*COS(N2)
!c        F = F*(4.0D0*F*F-3.0D0)
!c        RETURN
!c      END SUBROUTINE
!      
      SUBROUTINE randperm(PERM, IPERM, N)
        INTEGER (KIND=4), INTENT(IN) :: N
        INTEGER (KIND=4), INTENT(OUT) :: PERM(N), IPERM(N)
        INTEGER (KIND=4) :: I, J, K, ITEMP, JTEMP(N)
        REAL (KIND=8) :: TEMP

        JTEMP = 0
        DO I = 1, N
            DO WHILE (1==1)
                CALL RANDOM_NUMBER(TEMP)
                ITEMP = CEILING(TEMP*(N-I+1))
                IF (ITEMP /= 0) EXIT
            END DO
            K = 0
            DO J = 1, N
                IF (JTEMP(J) == 0) K = K+1
                IF (K == ITEMP) EXIT
            END DO
            PERM(I) = J
            IPERM(J) = I
            JTEMP(J) = 1
        END DO
      END SUBROUTINE

      SUBROUTINE init_random_seed()
      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
                                          
      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))
                  
      CALL SYSTEM_CLOCK(COUNT=clock)
                  
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)
             
      DEALLOCATE(seed)
      END SUBROUTINE


!      SUBROUTINE SOLVEPLE(A, b, N)
!        INTEGER(KIND=4), INTENT(IN) :: N
!        REAL(KIND=8), INTENT(INOUT) :: A(N,N), b(N)
!
!        CHARACTER UPLO 
!        INTEGER(KIND=4) :: LDA, LDB, NRHS
!        INTEGER(KIND=4) :: INFO
!
!        UPLO = 'U'
!        LDA = N
!        LDB = N
!        NRHS = 1
!        
!        CALL DPOTRF( UPLO, N, A, LDA, INFO )
!        IF(INFO /= 0) THEN
!            WRITE(*,*) "DPOTRF exites with INFO = ", INFO
!        END IF
!
!        
!        CALL DPOTRS( UPLO, N, NRHS, A, LDA, b, LDB, INFO )
!        IF(INFO /= 0) THEN
!            WRITE(*,*) "DPOTRS exites with INFO = ", INFO
!        END IF
!
!      END SUBROUTINE






