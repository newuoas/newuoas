      SUBROUTINE SETUP (PROB, N, X, RHOBEG)
        IMPLICIT NONE
        CHARACTER(LEN=20), INTENT(IN) :: PROB
        INTEGER(KIND=4), INTENT(IN) :: N
        REAL(KIND=8), INTENT(OUT) :: X(N), RHOBEG
        INTEGER(KIND=4) :: I
        REAL(KIND=8) :: SCALE
        CHARACTER(LEN=7) :: NAME7 ! For Problems DIXMAAN*.
        NAME7 = PROB(1:7)
        RHOBEG = 1.0D0

        IF (TRIM(PROB) .EQ. 'ARGLINA') THEN
            WRITE(*,*) "ARGLINA" 
            DO I = 1, N
              X(I) = 1.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'ARGLINA4') THEN
            WRITE(*,*) "ARGLINA4 (DIFFERENT FROM ARGLINA. Quartic 
     1instead of  Quadratic)"
            DO I = 1, N
              X(I) = 1.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'ARGLINB') THEN
            WRITE(*,*) "ARGLINB"
            DO I = 1, N
              X(I) = 1.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'ARGLINC') THEN
            WRITE(*,*) "ARGLINC"
            DO I = 1, N
              X(I) = 1.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'ARGTRIG') THEN
            WRITE(*,*) "ARGTRIG"
            DO I = 1, N
              X(I) = 1.0D0/DFLOAT(N)
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'ARWHEAD') THEN
            write(*,*) "ARWHEAD."
            DO I = 1, N
              !X(I) = 10.0D0
              X(I) = 1.0D0
            END DO
            RHOBEG = 0.5D0
        ELSE IF (TRIM(PROB) .EQ. 'BDQRTIC') THEN
            write(*,*) "BDQRTIC."
            DO I = 1, N
              X(I) = 1.0D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'BDQRTICP') THEN
            write(*,*) "BDQRTICP."
            DO I = 1, N
              X(I) = 1.0D0
            END DO
            RHOBEG = 0.5D0
        ELSE IF (TRIM(PROB) .EQ. 'BDVALUE') THEN
            WRITE(*,*) "BDVALUE"
            DO I = 1, N
              X(I) = DFLOAT(I*(I-N-1))/DFLOAT((N+1)*(N+1)) 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'BROWNAL') THEN
            WRITE(*,*) "BROWNAL"
            DO I = 1, N
              X(I) = 0.5D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'BROYDN3D') THEN
            WRITE(*,*) "BROYDN3D"
            DO I = 1, N
              X(I) = -1.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'BROYDN7D') THEN
            WRITE(*,*) "BROYDN7D"
            DO I = 1, N
              X(I) = -1.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'BRYBND') THEN
            WRITE(*,*) "BRYBND"
            DO I = 1, N
              X(I) = -1.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'CHAINWOO') THEN
            WRITE(*,*) "CHAINWOO"
            DO I = 1, N
                IF (I ==1 .OR. I==3) THEN 
                    X(I) = -3.0D0
                ELSE IF (I ==2 .OR. I==4) THEN 
                    X(I) = -1.0D0
                ELSE
                    X(I) = -2.0D0 
                END IF
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'CHEBQUAD') THEN
            write(*,*) "CHEBQUAD."
            DO I = 1, N
              X(I)=DFLOAT(I)/DFLOAT(N+1)
            END DO
            RHOBEG=0.2D0*X(1)
        ELSE IF (TRIM(PROB) .EQ. 'CHPOWELLB') THEN
            WRITE(*,*) "CHPOWELLB"
            DO I = 1, N
              X(I) = 1.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'CHPOWELLS') THEN
            WRITE(*,*) "CHPOWELLS"
            DO I = 1, N
              X(I) = 1.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'CHNROSNB') THEN
            write(*,*) "CHNROSNB."
            DO I = 1, N
              X(I) = -1.0D0
            END DO
            !RHOBEG = 0.5D0
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'CHROSEN') THEN
            write(*,*) "CHROSEN."
            DO I = 1, N
C              X(I) = -10.0D0
              X(I) = -1.0D0
            END DO
            RHOBEG = 0.5D0
        ELSE IF (TRIM(PROB) .EQ. 'COSINE') THEN
            write(*,*) "COSINE."
            DO I = 1, N
              X(I) = 1.0D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'CRAGGLVY') THEN
            WRITE(*,*) "CRAGGLVY"
            DO I = 2, N
              X(I) = 2.0D0 
            END DO
            X(1) = 1.0D0
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'CURLY10') THEN
            WRITE(*,*) "CURLY10"
            DO I = 1, N
                X(I) = (1.0D-4)/(DFLOAT(N+1)) 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'CURLY20') THEN
            WRITE(*,*) "CURLY20"
            DO I = 1, N
                X(I) = (1.0D-4)/(DFLOAT(N+1)) 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'CURLY30') THEN
            WRITE(*,*) "CURLY30"
            DO I = 1, N
                X(I) = (1.0D-4)/(DFLOAT(N+1)) 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'CUBE') THEN
            WRITE(*,*) "CUBE"
            X(1) = -1.2D0
            IF (N >= 2) X(2) = 1.0D0
            DO I = 3, N
                X(I) = X(I-2)
            END DO
            RHOBEG = 1.0D0
        ELSE IF (NAME7 .EQ. 'DIXMAAN') THEN
            WRITE(*,*) TRIM(PROB)
            DO I = 1, N
              X(I) = 2.0D0 
            END DO
            RHOBEG = 1.0D0

        ELSE IF (TRIM(PROB) .EQ. 'DQRTIC') THEN
            WRITE(*,*) "DQRTIC"
            DO I = 1, N
              X(I) = 2.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'EDENSCH') THEN
            WRITE(*,*) "EDENSCH"
            DO I = 1, N
              X(I) = 0.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'EG2') THEN
            WRITE(*,*) "EG2"
            DO I = 1, N
              X(I) = 1.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'ENGVAL1') THEN
            WRITE(*,*) "ENGVAL1"
            DO I = 1, N
              X(I) = 2.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'ERRINROS') THEN
            write(*,*) "ERRINROS."
            DO I = 1, N
              X(I) = -1.0D0
            END DO
            !RHOBEG = 0.5D0
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'EXTROSNB') THEN
            write(*,*) "EXTROSNB."
            DO I = 1, N
              X(I) = -1.0D0
            END DO
            !RHOBEG = 0.5D0
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'EXTTET') THEN
            WRITE(*,*) "EXTTET"
            DO I = 1, N
              X(I) = 0.1D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'FIROSE') THEN
            write(*,*) "FIROSE."
            DO I = 1, N
C              X(I) = -10.0D0
              X(I) = -1.0D0
            END DO
            RHOBEG = 0.5D0
        ELSE IF (TRIM(PROB) .EQ. 'FLETCBV2') THEN
            WRITE(*,*) "FLETCBV2"
            DO I = 1, N
                X(I) = DFLOAT(I)/DFLOAT(N+1)
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'FLETCBV3') THEN
            WRITE(*,*) "FLETCBV3"
            DO I = 1, N
                X(I) = DFLOAT(I)/DFLOAT(N+1)
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'FLETCHCR') THEN
            WRITE(*,*) "FLETCHCR"
            DO I = 1, N
                X(I) = 0.0D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'FREUROTH') THEN
            WRITE(*,*) "FREUROTH"
            X(1) = 0.5D0
            IF (N>=2) THEN
                X(2) = -2.0D0
            END IF
            DO I = 3, N
                X(I) = 0.0D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'GENBROWN') THEN
            WRITE(*,*) "GENBROWN"
            DO I = 1, N
              X(I) = 1.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'GENHUMPS') THEN
            WRITE(*,*) "GENHUMPS"
            X(1) = -506.0D0
            DO I = 2, N
              X(I) = -506.2D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'GENROSE') THEN
            WRITE(*,*) "GENROSE"
            DO I = 1, N
                X(I) = DFLOAT(I)/DFLOAT(N+1) 
            END DO
C            RHOBEG = 0.5D0
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'INDEF') THEN
            WRITE(*,*) "INDEF"
            DO I = 1, N
                X(I) = DFLOAT(I)/DFLOAT(N+1) 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'INTEGREQ') THEN
            WRITE(*,*) "INTEGREQ"
            DO I = 1, N
              X(I) = DFLOAT(I*(I-N-1))/DFLOAT((N+1)*(N+1)) 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'LIARWHD') THEN
            WRITE(*,*) "LIARWHD"
            DO I = 1, N
              X(I) = 4.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'LILIFUN3') THEN
            WRITE(*,*) "LILIFUN3"
            DO I = 1, N
              X(I) = 1.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'LILIFUN4') THEN
            WRITE(*,*) "LILIFUN4"
            DO I = 1, N
              X(I) = 1.0D0 
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'MOREBV') THEN
            WRITE(*,*) "MOREBV"
            DO I = 1, N
              X(I) = DFLOAT(I*(I-N-1))/DFLOAT((N+1)*(N+1)) 
            END DO
        ELSE IF (TRIM(PROB) .EQ. 'MOREBVL') THEN
            RHOBEG = 1.0D0
            WRITE(*,*) "MOREBVL"
            DO I = 1, N
              X(I) = 0.0D0 ! Suggested by 
                            ! Luksan, Matonoha, Vlcek, "Modified CUTE 
                            ! Probelems for Unconstrained Optimization.
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'NCB20') THEN
            WRITE(*,*) "NCB20"
            IF (N<20) THEN
                DO I = 1, N
                    X(I) = 0.0D0
                END DO
            ELSE 
                DO I = 1, N-10
                    X(I) = 0.0D0
                END DO
                DO I = N-9, N
                    X(I) = 1.0D0
                END DO
                    
            END IF
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'NCB20B') THEN
            WRITE(*,*) "NCB20."
            DO I = 1, N
                X(I) = 0.0D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'NONCVXU2') THEN
            WRITE(*,*) "NONCVXU2"
            DO I = 1,N
                X(I) = DFLOAT(I)
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'NONCVXUN') THEN
            WRITE(*,*) "NONCVXUN"
            DO I = 1,N
                X(I) = DFLOAT(I)
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'NONDIA') THEN
            WRITE(*,*) "NONDIA"
            DO I = 1, N
                X(I) = -1.0D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'NONDQUAR') THEN
            WRITE(*,*) "NONDQUAR"
            DO I = 1, N
                IF (MOD(I,2) == 1) X(I) = 1.0D0    
                IF (MOD(I,2) == 0) X(I) = -1.0D0    
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'PENALTY1') THEN
            write(*,*) "PENALTY1."
            DO I = 1, N
              X(I) = DFLOAT(I)
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'PENALTY2') THEN
            write(*,*) "PENALTY2."
            DO I = 1, N
              X(I) = 0.5D0
            END DO
            RHOBEG = 0.1D0
        ELSE IF (TRIM(PROB) .EQ. 'PENALTY3') THEN
            write(*,*) "PENALTY3."
            DO I = 1, N
              X(I) = DFLOAT(I)/DFLOAT(N+1)
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'PENALTY3P') THEN
            write(*,*) "PENALTY3P."
            DO I = 1, N
              X(I) = 0.0D0  ! Suggested by Powell, "The NEWUOA Software 
C            for unconstrained optimization without
C             derivatives", 2004, section 8.                       
            END DO
            RHOBEG = 0.1D0
        ELSE IF(TRIM(PROB) .EQ. 'POWELLSG') THEN
            WRITE(*,*) "POWELLSG"
            DO I = 1, N
                IF (MOD(I,4) == 1) X(I) = 3.0D0
                IF (MOD(I,4) == 2) X(I) = -1.0D0
                IF (MOD(I,4) == 3) X(I) = 0.0D0
                IF (MOD(I,4) == 4) X(I) = 1.0D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'POWER') THEN
            write(*,*) "POWER"
            DO I = 1, N
              X(I) = 1.0D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'ROSENBROCK') THEN
            write(*,*) "ROSENBROCK."
            DO I = 1, N
              X(I) = -1.0D0
            END DO
            RHOBEG = 0.5D0
        ELSE IF (TRIM(PROB) .EQ. 'SBRYBND') THEN
            WRITE(*,*) "SBRYBND"
            SCALE = 1.2D1
            DO I = 1, N
                X(I) = EXP(SCALE*DFLOAT(1-I)/DFLOAT(N-1))
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'SBRYBNDL') THEN
            WRITE(*,*) "SBRYBNDL"
            SCALE = 1.2D1
            DO I = 1, N
                X(I) = EXP(SCALE*DFLOAT(1-I)/DFLOAT(N-1))
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'SCHMVETT') THEN
            WRITE(*,*) "SCHMVETT"
            DO I = 1, N
                X(I) = 3.0D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'SCOSINE') THEN
            WRITE(*,*) "SCOSINE"
            SCALE = 1.2D1
            DO I = 1, N
                X(I) = EXP(SCALE*DFLOAT(1-I)/DFLOAT(N-1))
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'SCOSINEL') THEN
            WRITE(*,*) "SCOSINEL"
            SCALE = 1.2D1
            DO I = 1, N
                X(I) = EXP(SCALE*DFLOAT(1-I)/DFLOAT(N-1))
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'SEROSE') THEN
            write(*,*) "SEROSE."
            DO I = 1, N
C             X(I) = -10.0D0
              X(I) = -1.0D0
            END DO
            RHOBEG = 0.5D0
        ELSE IF (TRIM(PROB) .EQ. 'SINQUAD') THEN
            write(*,*) "SINQUAD."
            DO I = 1, N
              X(I) = 0.1D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'SPARSINE') THEN
            write(*,*) "SPARSINE."
            DO I = 1, N
              X(I) = 0.5D0
            END DO
            RHOBEG = 1.0D0

        ELSE IF (TRIM(PROB) .EQ. 'SPARSQUR') THEN
            write(*,*) "SPARSQUR."
            DO I = 1, N
              X(I) = 0.5D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'SPMSRTLS') THEN
            WRITE(*,*) "SPMSRTLS"
            DO I = 1, N
                X(I) = 0.2D0*SIN(DFLOAT(I)**2)
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'SPHRPTS') THEN
            write(*,*) "SPHRPTS."
            do I=1,n/2
              X(2*I) = 0.0D0
              X(2*I-1) = 8.0D0*ACOS(0.0D0)*DFLOAT(I)/DFLOAT(n)
            end do
            RHOBEG = 1.0D0/DFLOAT(n)
        ELSE IF (TRIM(PROB) .EQ. 'SROSENBR') THEN
            WRITE(*,*) "SROSENBR"
            DO I = 1, N
                IF (MOD(I,2) == 1) X(I) = -1.2D0
                IF (MOD(I,2) ==0) X(I) = 1.0D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'TOINTGSS') THEN
            WRITE(*,*) "TOINTGSS"
            DO I = 1, N
                X(I) = 3.0D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'TOINTTRIG') THEN
            WRITE(*,*) "TOINTTRIG"
            DO I = 1, N
                X(I) = 1.0D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'TQUARTIC') THEN
            WRITE(*,*) "TQUARTIC"
            DO I = 1, N
                X(I) = 0.1D0
            END DO
            RHOBEG = 1.0D0
        ELSE IF (TRIM(PROB) .EQ. 'TRIROSE') THEN
            write(*,*) "TRIROSE."
            DO I = 1, N
C             X(I) = -10.0D0
              X(I) = -1.0D0
            END DO
            RHOBEG = 0.5D0
        ELSE IF (TRIM(PROB) .EQ. 'VARDIM') THEN
            write(*,*) "VARDIM."
            DO I = 1, N
               X(I) = DFLOAT(N-I)/DFLOAT(N)
            END DO
            RHOBEG = 1.0D0/DFLOAT(2*N)
        ELSE IF (TRIM(PROB) .EQ. 'WOODS') THEN
            WRITE(*,*) "WOODS"
            DO I = 1, N
                IF (MOD(I,2) == 1) X(I) = -3.0D0
                IF (MOD(I,2) ==0) X(I) = -1.0D0
            END DO
            RHOBEG = 1.0D0
        ELSE
            WRITE(*,*)"ERROR in MAIN: Unknown problem name:",TRIM(PROB)
            STOP
        END IF 
      END
