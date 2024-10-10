      SUBROUTINE NEWUOA1(N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),W(*)
C
C     This subroutine seeks the least value of a function of many variables,
C     by a trust region method that forms quadratic models by interpolation.
C     There can be some freedom in the interpolation conditions, which is
C     taken up by minimizing the Frobenius norm of the change to the second
C     derivative of the quadratic model, beginning with a zero matrix. The
C     arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
C     NPT is the number of interpolation conditions. Its value must be in the
C       interval [N+2,(N+1)(N+2)/2].
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C       will be changed to the values that give the least calculated F.
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
C       RHOBEG should be about one tenth of the greatest expected change to a
C       variable, and RHOEND should indicate the accuracy that is required in
C       the final values of the variables.
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C       amount of printing. Specifically, there is no output if IPRINT=0 and
C       there is output only at the return if IPRINT=1. Otherwise, each new
C       value of RHO is printed, with the best vector of variables so far and
C       the corresponding value of the objective function. Further, each new
C       value of F with its variables are output if IPRINT=3.
C     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
C     The array W will be used for working space. Its length must be at least
C     (NPT+13)*(NPT+N)+3*N*(N+3)/2.
C
C     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
C     the value of the objective function for the variables X(1),X(2),...,X(N).
C
C     Partition the working space array, so that different parts of it can be
C     treated separately by the subroutine that performs the main calculation.
C
      NP=N+1
      NPTM=NPT-NP
      IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
          PRINT 10
   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in',
     1      ' the required interval')
          GO TO 20
      END IF
      NDIM=NPT+N
      IXB=1
      IXO=IXB+N
      IXN=IXO+N
      IXP=IXN+N
      IFV=IXP+N*NPT
      IGQ=IFV+NPT
      IHQ=IGQ+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ID=IZMAT+NPT*NPTM
      IVL=ID+N
      IW=IVL+NDIM
C
C     The above settings provide a partition of W for subroutine NEWUOB.
C     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
C     W plus the space that is needed by the last array of NEWUOB.
C
      CALL NEWUOB1 (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
     1  W(IXO),W(IXN),W(IXP),W(IFV),W(IGQ),W(IHQ),W(IPQ),W(IBMAT),
     2  W(IZMAT),NDIM,W(ID),W(IVL),W(IW))
   20 RETURN
      END 



      SUBROUTINE NEWUOA2(N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),W(*)
C
C     This subroutine seeks the least value of a function of many variables,
C     by a trust region method that forms quadratic models by interpolation.
C     There can be some freedom in the interpolation conditions, which is
C     taken up by minimizing the Frobenius norm of the change to the second
C     derivative of the quadratic model, beginning with a zero matrix. The
C     arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
C     NPT is the number of interpolation conditions. Its value must be in the
C       interval [N+2,(N+1)(N+2)/2].
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C       will be changed to the values that give the least calculated F.
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
C       RHOBEG should be about one tenth of the greatest expected change to a
C       variable, and RHOEND should indicate the accuracy that is required in
C       the final values of the variables.
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C       amount of printing. Specifically, there is no output if IPRINT=0 and
C       there is output only at the return if IPRINT=1. Otherwise, each new
C       value of RHO is printed, with the best vector of variables so far and
C       the corresponding value of the objective function. Further, each new
C       value of F with its variables are output if IPRINT=3.
C     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
C     The array W will be used for working space. Its length must be at least
C     (NPT+13)*(NPT+N)+3*N*(N+3)/2.
C
C     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
C     the value of the objective function for the variables X(1),X(2),...,X(N).
C
C     Partition the working space array, so that different parts of it can be
C     treated separately by the subroutine that performs the main calculation.
C
      NP=N+1
      NPTM=NPT-NP
      IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
          PRINT 10
   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in',
     1      ' the required interval')
          GO TO 20
      END IF
      NDIM=NPT+N
      IXB=1
      IXO=IXB+N
      IXN=IXO+N
      IXP=IXN+N
      IFV=IXP+N*NPT
      IGQ=IFV+NPT
      IHQ=IGQ+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ID=IZMAT+NPT*NPTM
      IVL=ID+N
      IW=IVL+NDIM
C
C     The above settings provide a partition of W for subroutine NEWUOB.
C     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
C     W plus the space that is needed by the last array of NEWUOB.
C
      CALL NEWUOB2 (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
     1  W(IXO),W(IXN),W(IXP),W(IFV),W(IGQ),W(IHQ),W(IPQ),W(IBMAT),
     2  W(IZMAT),NDIM,W(ID),W(IVL),W(IW))
   20 RETURN
      END 



      SUBROUTINE NEWUOA3(N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),W(*)
C
C     This subroutine seeks the least value of a function of many variables,
C     by a trust region method that forms quadratic models by interpolation.
C     There can be some freedom in the interpolation conditions, which is
C     taken up by minimizing the Frobenius norm of the change to the second
C     derivative of the quadratic model, beginning with a zero matrix. The
C     arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
C     NPT is the number of interpolation conditions. Its value must be in the
C       interval [N+2,(N+1)(N+2)/2].
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C       will be changed to the values that give the least calculated F.
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
C       RHOBEG should be about one tenth of the greatest expected change to a
C       variable, and RHOEND should indicate the accuracy that is required in
C       the final values of the variables.
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C       amount of printing. Specifically, there is no output if IPRINT=0 and
C       there is output only at the return if IPRINT=1. Otherwise, each new
C       value of RHO is printed, with the best vector of variables so far and
C       the corresponding value of the objective function. Further, each new
C       value of F with its variables are output if IPRINT=3.
C     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
C     The array W will be used for working space. Its length must be at least
C     (NPT+13)*(NPT+N)+3*N*(N+3)/2.
C
C     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
C     the value of the objective function for the variables X(1),X(2),...,X(N).
C
C     Partition the working space array, so that different parts of it can be
C     treated separately by the subroutine that performs the main calculation.
C
      NP=N+1
      NPTM=NPT-NP
      IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
          PRINT 10
   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in',
     1      ' the required interval')
          GO TO 20
      END IF
      NDIM=NPT+N
      IXB=1
      IXO=IXB+N
      IXN=IXO+N
      IXP=IXN+N
      IFV=IXP+N*NPT
      IGQ=IFV+NPT
      IHQ=IGQ+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ID=IZMAT+NPT*NPTM
      IVL=ID+N
      IW=IVL+NDIM
C
C     The above settings provide a partition of W for subroutine NEWUOB.
C     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
C     W plus the space that is needed by the last array of NEWUOB.
C
      CALL NEWUOB3 (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
     1  W(IXO),W(IXN),W(IXP),W(IFV),W(IGQ),W(IHQ),W(IPQ),W(IBMAT),
     2  W(IZMAT),NDIM,W(ID),W(IVL),W(IW))
   20 RETURN
      END 

      SUBROUTINE NEWUOA4(N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),W(*)
C
C     This subroutine seeks the least value of a function of many variables,
C     by a trust region method that forms quadratic models by interpolation.
C     There can be some freedom in the interpolation conditions, which is
C     taken up by minimizing the Frobenius norm of the change to the second
C     derivative of the quadratic model, beginning with a zero matrix. The
C     arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
C     NPT is the number of interpolation conditions. Its value must be in the
C       interval [N+2,(N+1)(N+2)/2].
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C       will be changed to the values that give the least calculated F.
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
C       RHOBEG should be about one tenth of the greatest expected change to a
C       variable, and RHOEND should indicate the accuracy that is required in
C       the final values of the variables.
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C       amount of printing. Specifically, there is no output if IPRINT=0 and
C       there is output only at the return if IPRINT=1. Otherwise, each new
C       value of RHO is printed, with the best vector of variables so far and
C       the corresponding value of the objective function. Further, each new
C       value of F with its variables are output if IPRINT=3.
C     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
C     The array W will be used for working space. Its length must be at least
C     (NPT+13)*(NPT+N)+3*N*(N+3)/2.
C
C     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
C     the value of the objective function for the variables X(1),X(2),...,X(N).
C
C     Partition the working space array, so that different parts of it can be
C     treated separately by the subroutine that performs the main calculation.
C
      NP=N+1
      NPTM=NPT-NP
      IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
          PRINT 10
   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in',
     1      ' the required interval')
          GO TO 20
      END IF
      NDIM=NPT+N
      IXB=1
      IXO=IXB+N
      IXN=IXO+N
      IXP=IXN+N
      IFV=IXP+N*NPT
      IGQ=IFV+NPT
      IHQ=IGQ+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ID=IZMAT+NPT*NPTM
      IVL=ID+N
      IW=IVL+NDIM
C
C     The above settings provide a partition of W for subroutine NEWUOB.
C     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
C     W plus the space that is needed by the last array of NEWUOB.
C
      CALL NEWUOB4 (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
     1  W(IXO),W(IXN),W(IXP),W(IFV),W(IGQ),W(IHQ),W(IPQ),W(IBMAT),
     2  W(IZMAT),NDIM,W(ID),W(IVL),W(IW))
   20 RETURN
      END 


C      SUBROUTINE NEWUOA5(N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
C      IMPLICIT REAL*8 (A-H,O-Z)
C      DIMENSION X(*),W(*)
CC
CC     This subroutine seeks the least value of a function of many variables,
CC     by a trust region method that forms quadratic models by interpolation.
CC     There can be some freedom in the interpolation conditions, which is
CC     taken up by minimizing the Frobenius norm of the change to the second
CC     derivative of the quadratic model, beginning with a zero matrix. The
CC     arguments of the subroutine are as follows.
CC
CC     N must be set to the number of variables and must be at least two.
CC     NPT is the number of interpolation conditions. Its value must be in the
CC       interval [N+2,(N+1)(N+2)/2].
CC     Initial values of the variables must be set in X(1),X(2),...,X(N). They
CC       will be changed to the values that give the least calculated F.
CC     RHOBEG and RHOEND must be set to the initial and final values of a trust
CC       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
CC       RHOBEG should be about one tenth of the greatest expected change to a
CC       variable, and RHOEND should indicate the accuracy that is required in
CC       the final values of the variables.
CC     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
CC       amount of printing. Specifically, there is no output if IPRINT=0 and
CC       there is output only at the return if IPRINT=1. Otherwise, each new
CC       value of RHO is printed, with the best vector of variables so far and
CC       the corresponding value of the objective function. Further, each new
CC       value of F with its variables are output if IPRINT=3.
CC     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
CC     The array W will be used for working space. Its length must be at least
CC     (NPT+13)*(NPT+N)+3*N*(N+3)/2.
CC
CC     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
CC     the value of the objective function for the variables X(1),X(2),...,X(N).
CC
CC     Partition the working space array, so that different parts of it can be
CC     treated separately by the subroutine that performs the main calculation.
CC
C      NP=N+1
C      NPTM=NPT-NP
C      IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
C          PRINT 10
C   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in',
C     1      ' the required interval')
C          GO TO 20
C      END IF
C      NDIM=NPT+N
C      IXB=1
C      IXO=IXB+N
C      IXN=IXO+N
C      IXP=IXN+N
C      IFV=IXP+N*NPT
C      IGQ=IFV+NPT
C      IHQ=IGQ+N
C      IPQ=IHQ+(N*NP)/2
C      IBMAT=IPQ+NPT
C      IZMAT=IBMAT+NDIM*N
C      ID=IZMAT+NPT*NPTM
C      IVL=ID+N
C      IW=IVL+NDIM
CC
CC     The above settings provide a partition of W for subroutine NEWUOB.
CC     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
CC     W plus the space that is needed by the last array of NEWUOB.
CC
C      CALL NEWUOB5 (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
C     1  W(IXO),W(IXN),W(IXP),W(IFV),W(IGQ),W(IHQ),W(IPQ),W(IBMAT),
C     2  W(IZMAT),NDIM,W(ID),W(IVL),W(IW))
C   20 RETURN
C      END 
C
