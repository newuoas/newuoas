      PROGRAM MAIN
      USE PROBLEM
      USE PERMUTE
      USE STATISTICS

      IMPLICIT NONE
      REAL(KIND=8) :: X(MAXDIM),XTEMP(MAXDIM),W(LWORK) 
      REAL(KIND=8) :: RHOBEG,RHOEND
      INTEGER(KIND=4) :: I, IPRINT, N, NPT, MAXFUN
      INTEGER(KIND=4) :: FN,FF,FFOPTREC, FDIM
      INTEGER(KIND=4) :: PERMS(MAXDIM,100), IPERMS(MAXDIM,100)
      INTEGER(KIND=4) :: NB, NE, ND, NSOLVER
      CHARACTER(LEN=20) :: CSOLVER, CNB, CNE, CND, CTEMP
      INTEGER(KIND=4) :: ARGC
      INTEGER(KIND=4) :: SOLVER

      IPRINT=3
      RHOEND=1.0D-6
      MAXFUN = MAXFUNEVL

      FN = 20
      FF = 21
      FFOPTREC = 22
      FDIM = 23
      
      ARGC = IARGC()
      IF (ARGC == 3) THEN 
          CALL GETARG(1,PROB)
          CALL GETARG(2, CNB)
          CALL GETARG(3,CSOLVER)
          READ(CNB,*) NB
          READ(CSOLVER,*) NSOLVER
          NE = NB
          ND = 1
      ELSE IF (ARGC == 5) THEN
          CALL GETARG(1, PROB)
          CALL GETARG(2, CNB)
          CALL GETARG(3, CNE)
          CALL GETARG(4, CND)
          CALL GETARG(5,CSOLVER)
          READ(CNB,*) NB
          READ(CNE,*) NE
          READ(CND,*) ND 
          READ(CSOLVER,*) NSOLVER
      ELSE
          WRITE(*,*) "USAGE: test PROBLEM LOWERDIM UPPERDIM DIMSTEP SOLV
     1ERNUM, or test PROBLEM DIM SOLVERNUM."
          STOP
      END IF
      IF (NE > MAXDIM ) THEN
          WRITE(*,*) "The largest dimension supported is", MAXDIM
          STOP
      END IF
      IF (NSOLVER > 99 ) THEN
          WRITE(*,*) "The largest dimension supported is 99."
          STOP
      END IF

      CALL init_random_seed()
      DO I = NB,NE,ND
        CALL RANDPERM(PERMS(:,(I-NB)/ND+1),IPERMS(:,(I-NB)/ND+1),I)
      END DO


      DO SOLVER = 1, NSOLVER
        DO N=NB, NE, ND
          NPT=2*N+1
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          MAXFUN = 100*N
          !MAXFUN = 80*N
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF ((NPT+13)*(NPT+N)+3*N*(N+3)/2 > LWORK) THEN
              WRITE(*,*) "The work space is too small."
              WRITE(*,*) "The smallest work space you need is an real*8
     1        array of size", (NPT+13)*(NPT+N)+3*N*(N+3)/2
              ! See modules.f for the definition of LWORK.
          END IF

          PERM(1:N) = PERMS(1:N,(N-NB)/ND+1)
    
          CALL SETUP(PROB, N, XTEMP, RHOBEG)
          DO I = 1, N
              X(PERM(I)) = XTEMP(I)+
     1         10.0D0*SIN(DFLOAT(I))*DMAX1(1.0D0,DABS(XTEMP(I)))
          END DO
          
          IF (N == NB) THEN
            WRITE( CTEMP,'(I2)' ) SOLVER
            OPEN(FN, FILE = TRIM(PROB) //'.N'// TRIM(ADJUSTL(CTEMP)))
            OPEN(FF, FILE = TRIM(PROB) //'.F'// TRIM(ADJUSTL(CTEMP)))
            OPEN(FFOPTREC, FILE = TRIM(PROB) //'.FOPTREC'// 
     1       TRIM(ADJUSTL(CTEMP)))
            OPEN(FDIM, FILE = TRIM(PROB) //'.DIM')
          END IF
    
          WRITE(*,*) 'SOLVER'//TRIM(ADJUSTL(CTEMP))
          WRITE(*,"(//4X,'Results with N =',I2,' and NPT =',I3)") N,NPT
C          PRINT 20, N,NPT
C   20   FORMAT (//4X,'Results with N =',I2,' and NPT =',I3)
          NFEVL = 0
          FOPTIMAL = 1.0D300
          IF (SOLVER == 1) THEN
              CALL NEWUOA1 (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
              WRITE (FDIM,*) N
          ELSE IF (SOLVER == 2) THEN
              CALL NEWUOA2 (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
          ELSE IF (SOLVER == 3) THEN
              CALL NEWUOA3 (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
          ELSE IF (SOLVER == 4) THEN
              CALL NEWUOA4 (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
C          ELSE IF (SOLVER == 5) THEN
C              CALL NEWUOA5 (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
C          ELSE IF (SOLVER == 6) THEN
C              CALL NEWUOA6 (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
          END IF
          WRITE (FN,*) NFEVL
          WRITE (FF,*) FOPTIMAL
          WRITE (FFOPTREC,*) FOPTREC(1:NFEVL)
          IF (N+ND > NE) THEN
             CLOSE(UNIT = FN, STATUS = 'KEEP')
             CLOSE(UNIT = FF, STATUS = 'KEEP')
             CLOSE(UNIT = FFOPTREC, STATUS = 'KEEP')
             CLOSE(UNIT = FDIM, STATUS = 'KEEP')
          END IF
        END DO 
      
      END DO
      
      STOP
      END 


