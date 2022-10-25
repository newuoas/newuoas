      SUBROUTINE NEWUOB4 (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,XBASE,
     1  XOPT,XNEW,XPT,FVAL,GQ,HQ,PQ,BMAT,ZMAT,NDIM,D,VLAG,W)
      IMPLICIT NONE
CCCCC DUMMIES CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER(KIND=4), INTENT(IN) :: N, NDIM, IPRINT,
     1 MAXFUN
      INTEGER(KIND=4), INTENT(INOUT) :: NPT
      REAL(KIND=8), INTENT(IN) :: RHOBEG, RHOEND 
      REAL(KIND=8), INTENT(INOUT) :: X(N), XBASE(N), XOPT(N), 
     1 XNEW(N), XPT(NPT,N), FVAL(NPT), GQ(N), HQ((N*(N+1))/2), 
     1 PQ(NPT), BMAT(NDIM,N), ZMAT(NPT,NPT-N-1),D(N), VLAG(*),W(*)
C
C     The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical
C       to the corresponding arguments in SUBROUTINE NEWUOA.
C     XBASE will hold a shift of origin that should reduce the contributions
C       from rounding errors to values of the model and Lagrange functions.
C     XOPT will be set to the displacement from XBASE of the vector of
C       variables that provides the least calculated F so far.
C     XNEW will be set to the displacement from XBASE of the vector of
C       variables for the current calculation of F.
C     XPT will contain the interpolation point coordinates relative to XBASE.
C     FVAL will hold the values of F at the interpolation points.
C     GQ will hold the gradient of the quadratic model at XBASE.
C     HQ will hold the explicit second derivatives of the quadratic model.
C     PQ will contain the parameters of the implicit second derivatives of
C       the quadratic model.
C     BMAT will hold the last N columns of H.
C     ZMAT will hold the factorization of the leading NPT by NPT submatrix of
C       H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, where
C       the elements of DZ are plus or minus one, as specified by IDZ.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     D is reserved for trial steps from XOPT.
C     VLAG will contain the values of the Lagrange functions at a new point X.
C       They are part of a product that requires VLAG to be of length NDIM.
C     The array W will be used for working space. Its length must be at least
C       10*NDIM = 10*(NPT+N).
C
CCCCC VARIABLES USED BY NEWUOA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER(KIND=4) :: NP, NH, NPTM, NFTEST
      REAL(KIND=8) :: HALF, ONE, ZERO, TENTH 
      INTEGER(KIND=4) :: I, J, K, IH, IP, IPT, JP, JPT, ITEMP, IDZ, 
     1 ITEST, KNEW, KOPT, KSAVE, NF, NFM, NFMM, NFSAV, KTEMP 
      REAL(KIND=8) :: ALPHA, BETA, BSUM, DELTA, DETRAT, DIFF, DIFFA, 
     1 DIFFB, DIFFC, DISTSQ, DNORM, DSQ, DSTEP, DX, F, FBEG, FOPT, 
     1 FSAVE, GISQ, GQSQ, HDIAG, VQUAD, RATIO, CRVMIN, RECIP, RHO, 
     1 RHOSQ, SUM, SUMA, SUMB, SUMZ, TEMP, TEMPQ, XIPT, XJPT, XOPTSQ,
     1 RECIQ 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC  NEW VARIABLES CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER(KIND=4) :: ITERNUM, SUBSIZE, SUBN, SUBNPT, CALCULATED,
     1  NFACT, SUBMAXFUN, SUBSUCCESS, ISUBUNSUC, IMAX, NEWITER
      REAL(KIND=8) :: SUBDIR(N,50), SUBBS(N,50), SUBD(50), G(N), 
     1  GNORM, FCALCULATED, SUBRHOBEG, SUBRHOEND,SUBDIRN(50),
     1  DIR(N), DN, HXOPT(N), DLAST(N), SUBW(100000), FMAX, DIAG(N),
     1  MDIAG
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      NPT = 2*N+1 !!!???? 
      ITERNUM = 1
      ISUBUNSUC = 0
      NFACT = 0
      DNORM = RHOBEG
      NEWITER = 1
      CALCULATED = 0
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      TENTH=0.1D0
      ZERO=0.0D0
      NP=N+1
      NH=(N*NP)/2
      NPTM=NPT-NP
      NFTEST=MAX0(MAXFUN,1)
      RHO = RHOBEG
C
C     Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
C
  001 WRITE(*,*) "Iter: ", ITERNUM
      DO 20 J=1,N
      XBASE(J)=X(J)
      DO 10 K=1,NPT
   10 XPT(K,J)=ZERO
      DO 20 I=1,NDIM
   20 BMAT(I,J)=ZERO
      DO 30 IH=1,NH
   30 HQ(IH)=ZERO
      DO 40 K=1,NPT
      PQ(K)=ZERO
      DO 40 J=1,NPTM
   40 ZMAT(K,J)=ZERO
C
C     Begin the initialization procedure. NF becomes one more than the number
C     of function values so far. The coordinates of the displacement of the
C     next initial interpolation point from XBASE are set in XPT(NF,.).
C
      RHOSQ=RHO*RHO
      RECIP=ONE/RHOSQ
      RECIQ=DSQRT(HALF)/RHOSQ
      NF=0
   50 NFM=NF
      NFMM=NF-N
      NF=NF+1
      IF (NFM .LE. 2*N) THEN
          IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
              XPT(NF,NFM)=RHO
          ELSE IF (NFM .GT. N) THEN
              XPT(NF,NFMM)=-RHO
          END IF
      ELSE
          ITEMP=(NFMM-1)/N
          JPT=NFM-ITEMP*N-N
          IPT=JPT+ITEMP
          IF (IPT .GT. N) THEN
              ITEMP=JPT
              JPT=IPT-N
              IPT=ITEMP
          END IF
          XIPT=RHO
          IF (FVAL(IPT+NP) .LT. FVAL(IPT+1)) XIPT=-XIPT
          XJPT=RHO
          IF (FVAL(JPT+NP) .LT. FVAL(JPT+1)) XJPT=-XJPT
          XPT(NF,IPT)=XIPT
          XPT(NF,JPT)=XJPT
      END IF
C
C     Calculate the next value of F, label 70 being reached immediately
C     after this calculation. The least function value so far and its index
C     are required.
C
      DO 60 J=1,N
   60 X(J)=XPT(NF,J)+XBASE(J)
      GOTO 310
   70 FVAL(NF)=F
      IF (NF .EQ. 1) THEN
          FBEG=F
          FOPT=F
          KOPT=1
      ELSE IF (F .LT. FOPT) THEN
          FOPT=F
          KOPT=NF
      END IF
C
C     Set the nonzero initial elements of BMAT and the quadratic model in
C     the cases when NF is at most 2*N+1.
C
      IF (NFM .LE. 2*N) THEN
          IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
              GQ(NFM)=(F-FBEG)/RHO
              IF (NPT .LT. NF+N) THEN
                  BMAT(1,NFM)=-ONE/RHO
                  BMAT(NF,NFM)=ONE/RHO
                  BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
              END IF
          ELSE IF (NFM .GT. N) THEN
              BMAT(NF-N,NFMM)=HALF/RHO
              BMAT(NF,NFMM)=-HALF/RHO
              ZMAT(1,NFMM)=-RECIQ-RECIQ
              ZMAT(NF-N,NFMM)=RECIQ
              ZMAT(NF,NFMM)=RECIQ
              IH=(NFMM*(NFMM+1))/2
              TEMP=(FBEG-F)/RHO
              HQ(IH)=(GQ(NFMM)-TEMP)/RHO
              GQ(NFMM)=HALF*(GQ(NFMM)+TEMP)
          END IF
C
C     Set the off-diagonal second derivatives of the Lagrange functions and
C     the initial quadratic model.
C
      ELSE
          IH=(IPT*(IPT-1))/2+JPT
          IF (XIPT .LT. ZERO) IPT=IPT+N
          IF (XJPT .LT. ZERO) JPT=JPT+N
          ZMAT(1,NFMM)=RECIP
          ZMAT(NF,NFMM)=RECIP
          ZMAT(IPT+1,NFMM)=-RECIP
          ZMAT(JPT+1,NFMM)=-RECIP
          HQ(IH)=(FBEG-FVAL(IPT+1)-FVAL(JPT+1)+F)/(XIPT*XJPT)
      END IF
      IF (NF .LT. NPT) GOTO 50
C
C     Begin the iterative procedure, because the initial model is complete.
C
      DELTA=RHO
      IDZ=1
      DIFFA=ZERO
      DIFFB=ZERO
      ITEST=0
      XOPTSQ=ZERO
      DO 80 I=1,N
      XOPT(I)=XPT(KOPT,I)
   80 XOPTSQ=XOPTSQ+XOPT(I)**2
   90 NFSAV=NF
C
C     Generate the next trust region step and test its length. Set KNEW
C     to -1 if the purpose of the next F will be to improve the model.
C
  100 KNEW=0
      IF (ISUBUNSUC >= 1) THEN 
          IF (NEWITER == 1) THEN
             !!!!!!!!!!!!!!!! ?????????????????????!!!!!!!!!!!!!!!!!
             DELTA = DMAX1(RHO,1.0D-1*SUBRHOBEG) !!!! What should be the
                                                 !!!! first DELTA of the
                                                 !!!! "new" iterations?
             NEWITER = 0
          END IF
          GOTO 101
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBSIZE = 4 
      CALL TRSAPP (N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,D,W,W(NP),
     1  W(NP+N),W(NP+2*N),CRVMIN)
    
      SUBDIR = 0.0D0

      CALL HSSD(N, NPT, HXOPT, HQ, PQ, XPT, XOPT)
      G = GQ+HXOPT  ! This is necessary, since XOPT may not be zero
                    ! (i.e., the current best point may not be XBASE),
                    ! even if no step has been taken, as the function
                    ! might get a descent during the sampling. 
C      SUBDIR(:,1) = -G
      GNORM = SQRT(DOT_PRODUCT(G,G))
      DO I = 1,N 
          DIAG(I) = HQ(I*(I+1)/2)
          DO J = 1, NPT  
              DIAG(I) = DIAG(I)+PQ(J)*XPT(J,I)**2
          END DO  
      END DO                                                            
      MDIAG = MAXVAL(DABS(DIAG(1:N)))
      WRITE(*,*) "PRECONDITION."
C      TEMP = 1.0D-6*MDIAG
C      DO I = 1, N
C          IF (DIAG(I) < TEMP) THEN
C             WRITE(*,*) "SMALL DIAGONAL."
C             G(I)=G(I)*((DIAG(I)-TEMP)**2/TEMP**3
C     1        -(DIAG(I)-TEMP)/TEMP**2+1.0D0/TEMP)
C          ELSE
C             G(I) = G(I)/DIAG(I)
C          END IF
C      END DO
      TEMP = 1.0D-6*MDIAG
      DO I = 1, N
          IF (DIAG(I) < TEMP) THEN
              WRITE(*,*) "SMALL DIAGONAL."
             G(I)=G(I)*(2.0D0/TEMP-DIAG(I)/TEMP**2)
          ELSE
             G(I) = G(I)/DIAG(I)
          END IF
      END DO
C      TEMP = 1.0D-6*MDIAG
C      DO I = 1, N
C          IF (ABS(DIAG(I)) < TEMP) THEN
C              G(I) = G(I)*DIAG(I)/TEMP
C          ELSE
C              G(I) = G(I)/DIAG(I)
C          END IF
C      END DO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DELTA = DMAX1(DNORM,RHO) !!!!!! IMPORTANT 
                               !!!!!! What should be the DELTA for the 
                               !!!!!! trial step?
                               !!!!!! Maybe the angle between the real
                               !!!!!! step and this trial step is a 
                               !!!!!! good reference!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C      SUBDIR(:,1) = D  !!! NOTICE: The order and sign of the directions
C                       !!! are important; they make big difference,
C                       !!! because SUBNEWUOA will sample
C                       !!! along the coordinate directions, the first
C                       !!! point being x_0 + rhobeg*e_1, and the second
C                       !!! one being x_0 + rhobeg*e_1 or x_0-rhobeg*e_1,
C                       !!! according to the behavior of x_0 +rhobeg*e_1.
C                       !!! The order and sign of the directions will 
C                       !!! determine who is e_1.
C      SUBDIR(:,2) = -G
C      SUBDIR(:,3) = DLAST
C      SUBDIR(:,2) = D  !!! NOTICE: The order and sign of the directions
                       !!! are important; they make big difference,
                       !!! because SUBNEWUOA will sample along the 
                       !!! coordinate directions. 
      SUBDIR(:,3) = -G
C      SUBDIR(:,2) = 1.0D0
      SUBDIR(:,4) = DLAST
C Nemiroveski Direction 3 
C      FMAX = FOPT
C      IMAX = 0
C      DO I = 1, NPT
C          IF (FVAL(I) > FMAX) THEN
C            FMAX = FVAL(I)
C            IMAX = I
C          END IF
C      END DO
C      SUBDIR(:,4) = XOPT - XPT(IMAX,:)
C CAUTION!!!! XPT is NPT*N !!!!!! Thus XPT(:,IMAX) is WRONG!!!!!!


      DO I = 1, SUBSIZE
         SUBDIRN(I) = SQRT(DOT_PRODUCT(SUBDIR(:,I),SUBDIR(:,I)))
         If (SUBDIRN(I) >= 1.0D-6*DMAX1(RHO,RHOEND)) THEN
             SUBDIR(:,I) = SUBDIR(:,I)/SUBDIRN(I)
         ELSE
             SUBDIR(:,I) = ZERO
             SUBDIRN(I) = ZERO
         END IF
      END DO

      SUBN = 0
      DO I = 1, SUBSIZE
         DIR = SUBDIR(:,I)
         DN = SQRT(DOT_PRODUCT(DIR,DIR))
C         IF (DN > 1.0D-6*DMAX1(RHO,RHOEND)) THEN
         IF (DN > 1.0D-6) THEN  !!!! Should be a relative value, since
                                !!!! SUBDIR has been normalized in last
                                !!!! step.

C  CAUTION!!!! Should be ">" instead of ">=", as SUBDIR(:,I) might be
C  ZERO, which should not be included into the basis.
            SUBN = SUBN+1
            DIR = DIR/DN
            SUBBS(:,SUBN) = DIR 
            DO J = I+1, SUBSIZE
               SUBDIR(:,J) = SUBDIR(:,J)-
     1         DOT_PRODUCT(SUBDIR(:,J),DIR)*DIR/DOT_PRODUCT(DIR,DIR)
            END DO
         END IF
      END DO
     
      !!!!! WHAT IF SUBN=0 ????
      SUBNPT = 2*SUBN + 1
      !SUBNPT = (SUBN+1)*(SUBN+2)/2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBRHOBEG = MAX(RHO,DNORM,1.0D-1*SUBRHOBEG)  !!!!! IMPORTANT
      !SUBRHOBEG = MAX(RHO,DNORM)  !!!!! IMPORTANT
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !SUBRHOEND = 0.1D0*SUBRHOBEG
      SUBRHOEND = RHOEND
      !SUBRHOEND = DMAX1(1.0D-2*SUBRHOBEG, 1.0D-1*RHO)
      !SUBMAXFUN = MIN(MAXFUN-NF, 20*SUBN)
      SUBMAXFUN = MAXFUN-NFACT

      SUBD = 0.0D0
      X = XOPT + XBASE
      F = FOPT

      WRITE(*,*) "RHO", RHO
      WRITE(*,*) "DELTA", DELTA
      WRITE(*,*) "SUBRHOBEG", SUBRHOBEG
      WRITE(*,*) "SUBN",SUBN
      WRITE(*,*) "SUBDIRN"
      WRITE(*,*) SUBDIRN(1:SUBSIZE)
    
      CALL SUBNEWUOA(SUBN,SUBNPT,SUBD,SUBRHOBEG,SUBRHOEND,IPRINT,
     1 SUBMAXFUN,SUBW,N,X,F,NFACT,SUBSUCCESS,SUBBS(:,1:SUBN))
      D = ZERO
      DO I = 1, SUBN
          D = D + SUBD(I)*SUBBS(:,I)
      END DO
      DLAST = D  !?????? When should DLAST be updated?? When D is long
                 !!! enough??
      DNORM = SQRT(DOT_PRODUCT(D,D))

      WRITE(*,*) "DNORM", DNORM

      XOPT = XOPT + D
      FOPT = F
      X = X + D
      CALCULATED = 1
      FCALCULATED = F
      IF (DNORM <= RHOEND) THEN   !!!!???????
          ISUBUNSUC = ISUBUNSUC + 1
          SUBSUCCESS = 0 !!!!!!!!!!!!??????????
          GOTO 541
      END IF
C      IF (SUBSUCCESS == 0) THEN
C          ISUBUNSUC = ISUBUNSUC + 1
C      END IF
C      IF (DMAX1(DNORM, GNORM, RHO) <= RHOEND) THEN 
      IF (DMAX1(GNORM, RHO) <= RHOEND) THEN 
      !!!!!!!!!!!!!????????????????????????? 
          GOTO 541
      ELSE 
          ITERNUM = ITERNUM + 1
          IF (SUBSUCCESS == 0) THEN
              RHO = 0.1D0*RHO !!!!!!!!!??????????????
              !!!! What if when DNORM <= 0.5*RHO for 3 times, then RHO =
              !!!! 0.1*RHO??
              !!!! What if when 0.1D0*RHO <= RHOEND, then break to "New"
              !!!! iterations?
              WRITE(*,*) "SUBFAIL", DNORM
          ELSE
              !RHO = DMAX1(0.5D0*RHO,1.0D-1*RHOEND) !!!!!!!!!!????
               RHO = DMAX1(0.5D0*RHO,RHOEND) !!!!!!!!!!????
              !!!! What if when 0.1D0*RHO <= RHOEND, then break to "New"
              !!!! iterations?
          END IF
          GOTO 001
      END IF


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC





  101 CALL TRSAPP (N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,D,W,W(NP),
     1  W(NP+N),W(NP+2*N),CRVMIN)
      DSQ=ZERO
      DO 110 I=1,N
  110 DSQ=DSQ+D(I)**2
      DNORM=DMIN1(DELTA,DSQRT(DSQ))
      IF (DNORM .LT. HALF*RHO) THEN
          KNEW=-1
          DELTA=TENTH*DELTA
          RATIO=-1.0D0
          IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
          IF (NF .LE. NFSAV+2) GOTO 460
          TEMP=0.125D0*CRVMIN*RHO*RHO
          IF (TEMP .LE. DMAX1(DIFFA,DIFFB,DIFFC)) GOTO 460
          GOTO 490
      END IF
C
C     Shift XBASE if XOPT may be too far from XBASE. First make the changes
C     to BMAT that do not depend on ZMAT.
C
  120 IF (DSQ .LE. 1.0D-3*XOPTSQ) THEN
          TEMPQ=0.25D0*XOPTSQ
          DO 140 K=1,NPT
          SUM=ZERO
          DO 130 I=1,N
  130     SUM=SUM+XPT(K,I)*XOPT(I)
          TEMP=PQ(K)*SUM
          SUM=SUM-HALF*XOPTSQ
          W(NPT+K)=SUM
          DO 140 I=1,N
          GQ(I)=GQ(I)+TEMP*XPT(K,I)
          XPT(K,I)=XPT(K,I)-HALF*XOPT(I)
          VLAG(I)=BMAT(K,I)
          W(I)=SUM*XPT(K,I)+TEMPQ*XOPT(I)
          IP=NPT+I
          DO 140 J=1,I
  140     BMAT(IP,J)=BMAT(IP,J)+VLAG(I)*W(J)+W(I)*VLAG(J)
C
C     Then the revisions of BMAT that depend on ZMAT are calculated.
C
          DO 180 K=1,NPTM
          SUMZ=ZERO
          DO 150 I=1,NPT
          SUMZ=SUMZ+ZMAT(I,K)
  150     W(I)=W(NPT+I)*ZMAT(I,K)
          DO 170 J=1,N
          SUM=TEMPQ*SUMZ*XOPT(J)
          DO 160 I=1,NPT
  160     SUM=SUM+W(I)*XPT(I,J)
          VLAG(J)=SUM
          IF (K .LT. IDZ) SUM=-SUM
          DO 170 I=1,NPT
  170     BMAT(I,J)=BMAT(I,J)+SUM*ZMAT(I,K)
          DO 180 I=1,N
          IP=I+NPT
          TEMP=VLAG(I)
          IF (K .LT. IDZ) TEMP=-TEMP
          DO 180 J=1,I
  180     BMAT(IP,J)=BMAT(IP,J)+TEMP*VLAG(J)
C
C     The following instructions complete the shift of XBASE, including
C     the changes to the parameters of the quadratic model.
C
          IH=0
          DO 200 J=1,N
          W(J)=ZERO
          DO 190 K=1,NPT
          W(J)=W(J)+PQ(K)*XPT(K,J)
  190     XPT(K,J)=XPT(K,J)-HALF*XOPT(J)
          DO 200 I=1,J
          IH=IH+1
          IF (I .LT. J) GQ(J)=GQ(J)+HQ(IH)*XOPT(I)
          GQ(I)=GQ(I)+HQ(IH)*XOPT(J)
          HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
  200     BMAT(NPT+I,J)=BMAT(NPT+J,I)
          DO 210 J=1,N
          XBASE(J)=XBASE(J)+XOPT(J)
  210     XOPT(J)=ZERO
          XOPTSQ=ZERO
      END IF
C
C     Pick the model step if KNEW is positive. A different choice of D
C     may be made later, if the choice of D by BIGLAG causes substantial
C     cancellation in DENOM.
C
      IF (KNEW .GT. 0) THEN
          CALL BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KNEW,DSTEP,
     1      D,ALPHA,VLAG,VLAG(NPT+1),W,W(NP),W(NP+N))
      END IF
C
C     Calculate VLAG and BETA for the current choice of D. The first NPT
C     components of W_check will be held in W.
C
      DO 230 K=1,NPT
      SUMA=ZERO
      SUMB=ZERO
      SUM=ZERO
      DO 220 J=1,N
      SUMA=SUMA+XPT(K,J)*D(J)
      SUMB=SUMB+XPT(K,J)*XOPT(J)
  220 SUM=SUM+BMAT(K,J)*D(J)
      W(K)=SUMA*(HALF*SUMA+SUMB)
  230 VLAG(K)=SUM
      BETA=ZERO
      DO 250 K=1,NPTM
      SUM=ZERO
      DO 240 I=1,NPT
  240 SUM=SUM+ZMAT(I,K)*W(I)
      IF (K .LT. IDZ) THEN
          BETA=BETA+SUM*SUM
          SUM=-SUM
      ELSE
          BETA=BETA-SUM*SUM
      END IF
      DO 250 I=1,NPT
  250 VLAG(I)=VLAG(I)+SUM*ZMAT(I,K)
      BSUM=ZERO
      DX=ZERO
      DO 280 J=1,N
      SUM=ZERO
      DO 260 I=1,NPT
  260 SUM=SUM+W(I)*BMAT(I,J)
      BSUM=BSUM+SUM*D(J)
      JP=NPT+J
      DO 270 K=1,N
  270 SUM=SUM+BMAT(JP,K)*D(K)
      VLAG(JP)=SUM
      BSUM=BSUM+SUM*D(J)
  280 DX=DX+D(J)*XOPT(J)
      BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     If KNEW is positive and if the cancellation in DENOM is unacceptable,
C     then BIGDEN calculates an alternative model step, XNEW being used for
C     working space.
C
      IF (KNEW .GT. 0) THEN
          TEMP=ONE+ALPHA*BETA/VLAG(KNEW)**2
          IF (DABS(TEMP) .LE. 0.8D0) THEN
              CALL BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
     1          KNEW,D,W,VLAG,BETA,XNEW,W(NDIM+1),W(6*NDIM+1))
          END IF
      END IF
C
C     Calculate the next value of the objective function.
C
  290 DO 300 I=1,N
      XNEW(I)=XOPT(I)+D(I)
  300 X(I)=XBASE(I)+XNEW(I)
      NF=NF+1

  310 IF (CALCULATED == 1) THEN
          F = FCALCULATED
          CALCULATED = 0
C          NF = NF - 1
      ELSE
          NFACT = NFACT + 1
          IF (NFACT .GT. NFTEST) THEN
              NFACT=NFACT-1
              IF (IPRINT .GT. 0) PRINT 320
  320         FORMAT (/4X,'Return from NEWUOA because CALFUN has been',
     1      ' called MAXFUN times.')
              GOTO 530
          END IF
          CALL CALFUN (N,X,F)
      END IF
      IF (IPRINT .EQ. 3) THEN
          PRINT 330, NFACT,F,(X(I),I=1,N)
  330      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (NF .LE. NPT) GOTO 70
      IF (KNEW .EQ. -1) GOTO 530
C
C     Use the quadratic model to predict the change in F due to the step D,
C     and set DIFF to the error of this prediction.
C
      VQUAD=ZERO
      IH=0
      DO 340 J=1,N
      VQUAD=VQUAD+D(J)*GQ(J)
      DO 340 I=1,J
      IH=IH+1
      TEMP=D(I)*XNEW(J)+D(J)*XOPT(I)
      IF (I .EQ. J) TEMP=HALF*TEMP
  340 VQUAD=VQUAD+TEMP*HQ(IH)
      DO 350 K=1,NPT
  350 VQUAD=VQUAD+PQ(K)*W(K)
      DIFF=F-FOPT-VQUAD
      DIFFC=DIFFB
      DIFFB=DIFFA
      DIFFA=DABS(DIFF)
      IF (DNORM .GT. RHO) NFSAV=NF
C
C     Update FOPT and XOPT if the new F is the least value of the objective
C     function so far. The branch when KNEW is positive occurs if D is not
C     a trust region step.
C
      FSAVE=FOPT
      IF (F .LT. FOPT) THEN
          FOPT=F
          XOPTSQ=ZERO
          DO 360 I=1,N
          XOPT(I)=XNEW(I)
  360     XOPTSQ=XOPTSQ+XOPT(I)**2
      END IF
      KSAVE=KNEW
      IF (KNEW .GT. 0) GOTO 410
C
C     Pick the next value of DELTA after a trust region step.
C
      IF (VQUAD .GE. ZERO) THEN
          IF (IPRINT .GT. 0) PRINT 370
  370     FORMAT (/4X,'Return from NEWUOA because a trust',
     1      ' region step has failed to reduce Q.')
          GOTO 530
      END IF
      RATIO=(F-FSAVE)/VQUAD
      IF (RATIO .LE. TENTH) THEN
          DELTA=HALF*DNORM
      ELSE IF (RATIO. LE. 0.7D0) THEN
          DELTA=DMAX1(HALF*DELTA,DNORM)
      ELSE
          DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
      END IF
      IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
C
C     Set KNEW to the index of the next interpolation point to be deleted.
C
      RHOSQ=DMAX1(TENTH*DELTA,RHO)**2
      KTEMP=0
      DETRAT=ZERO
      IF (F .GE. FSAVE) THEN
          KTEMP=KOPT
          DETRAT=ONE
      END IF
      DO 400 K=1,NPT
      HDIAG=ZERO
      DO 380 J=1,NPTM
      TEMP=ONE
      IF (J .LT. IDZ) TEMP=-ONE
  380 HDIAG=HDIAG+TEMP*ZMAT(K,J)**2
      TEMP=DABS(BETA*HDIAG+VLAG(K)**2)
      DISTSQ=ZERO
      DO 390 J=1,N
  390 DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
      IF (DISTSQ .GT. RHOSQ) TEMP=TEMP*(DISTSQ/RHOSQ)**3
      IF (TEMP .GT. DETRAT .AND. K .NE. KTEMP) THEN
          DETRAT=TEMP
          KNEW=K
      END IF
  400 CONTINUE
      IF (KNEW .EQ. 0) GOTO 460
C
C     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
C     can be moved. Begin the updating of the quadratic model, starting
C     with the explicit second derivative term.
C
  410 CALL UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
      FVAL(KNEW)=F
      IH=0
      DO 420 I=1,N
      TEMP=PQ(KNEW)*XPT(KNEW,I)
      DO 420 J=1,I
      IH=IH+1
  420 HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
      PQ(KNEW)=ZERO
C
C     Update the other second derivative parameters, and then the gradient
C     vector of the model. Also include the new interpolation point.
C
      DO 440 J=1,NPTM
      TEMP=DIFF*ZMAT(KNEW,J)
      IF (J .LT. IDZ) TEMP=-TEMP
      DO 440 K=1,NPT
  440 PQ(K)=PQ(K)+TEMP*ZMAT(K,J)
      GQSQ=ZERO
      DO 450 I=1,N
      GQ(I)=GQ(I)+DIFF*BMAT(KNEW,I)
      GQSQ=GQSQ+GQ(I)**2
  450 XPT(KNEW,I)=XNEW(I)
C
C     If a trust region step makes a small change to the objective function,
C     then calculate the gradient of the least Frobenius norm interpolant at
C     XBASE, and store it in W, using VLAG for a vector of right hand sides.
C
      IF (KSAVE .EQ. 0 .AND. DELTA .EQ. RHO) THEN
          IF (RATIO .GT. 1.0D-2) THEN
C          IF (DABS(RATIO) .GT. 1.0D-2) THEN
              ITEST=0
          ELSE
              DO 700 K=1,NPT
  700         VLAG(K)=FVAL(K)-FVAL(KOPT)
              GISQ=ZERO
              DO 720 I=1,N
              SUM=ZERO
              DO 710 K=1,NPT
  710         SUM=SUM+BMAT(K,I)*VLAG(K)
              GISQ=GISQ+SUM*SUM
  720         W(I)=SUM
C
C     Test whether to replace the new quadratic model by the least Frobenius
C     norm interpolant, making the replacement if the test is satisfied.
C
              ITEST=ITEST+1
              IF (GQSQ .LT. 1.0D2*GISQ) ITEST=0
              IF (ITEST .GE. 3) THEN
                write(*,*)"NF = "
                write(*,*)NF
                write(*,*)"RESTART"
                  
                  DO 730 I=1,N
  730             GQ(I)=W(I)
                  DO 740 IH=1,NH
  740             HQ(IH)=ZERO
                  DO 760 J=1,NPTM
                  W(J)=ZERO
                  DO 750 K=1,NPT
  750             W(J)=W(J)+VLAG(K)*ZMAT(K,J)
  760             IF (J .LT. IDZ) W(J)=-W(J)
                  DO 770 K=1,NPT
                  PQ(K)=ZERO
                  DO 770 J=1,NPTM
  770             PQ(K)=PQ(K)+ZMAT(K,J)*W(J)
                  ITEST=0
              END IF
          END IF
      END IF
      IF (F .LT. FSAVE) KOPT=KNEW
C
C     If a trust region step has provided a sufficient decrease in F, then
C     branch for another trust region calculation. The case KSAVE>0 occurs
C     when the new function value was calculated by a model step.
C
      IF (F .LE. FSAVE+TENTH*VQUAD) GOTO 100
      IF (KSAVE .GT. 0) GOTO 100
C
C     Alternatively, find out if the interpolation points are close enough
C     to the best point so far.
C
      KNEW=0
  460 DISTSQ=4.0D0*DELTA*DELTA
      DO 480 K=1,NPT
      SUM=ZERO
      DO 470 J=1,N
  470 SUM=SUM+(XPT(K,J)-XOPT(J))**2
      IF (SUM .GT. DISTSQ) THEN
          KNEW=K
          DISTSQ=SUM
      END IF
  480 CONTINUE
C
C     If KNEW is positive, then set DSTEP, and branch back for the next
C     iteration, which will generate a "model step".
C
      IF (KNEW .GT. 0) THEN
          DSTEP=DMAX1(DMIN1(TENTH*DSQRT(DISTSQ),HALF*DELTA),RHO)
          DSQ=DSTEP*DSTEP
          GOTO 120
      END IF
      IF (RATIO .GT. ZERO) GOTO 100
      IF (DMAX1(DELTA,DNORM) .GT. RHO) GOTO 100
C
C     The calculations with the current value of RHO are complete. Pick the
C     next values of RHO and DELTA.
C
  490 IF (RHO .GT. RHOEND) THEN
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          IF (RATIO .LE. 16.0D0) THEN
              RHO=RHOEND
          ELSE IF (RATIO .LE. 250.0D0) THEN
              RHO=DSQRT(RATIO)*RHOEND
          ELSE
              RHO=TENTH*RHO
          END IF
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT .GE. 2) THEN
              IF (IPRINT .GE. 3) PRINT 500
  500         FORMAT (5X)
              PRINT 510, RHO,NF
  510         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     1          ' function values =',I6)
              PRINT 520, FOPT,(XBASE(I)+XOPT(I),I=1,N)
  520         FORMAT (4X,'Least value of F =',1PD23.15,9X,
     1          'The corresponding X is:'/(2X,5D15.6))
          END IF
          GOTO 90
      END IF
C
C     Return from the calculation, after another Newton-Raphson step, if
C     it is too short to have been tried before.
C
      IF (KNEW .EQ. -1) GOTO 290
  530 IF (FOPT .LE. F) THEN
          DO 540 I=1,N
  540     X(I)=XBASE(I)+XOPT(I)
          F=FOPT
      END IF
  541 IF (IPRINT .GE. 1) THEN
          PRINT 550, NFACT
  550     FORMAT (/4X,'At the return from NEWUOA',5X,
     1      'Number of function values =',I6)
          PRINT 520, F,(X(I),I=1,N)
      END IF
      RETURN
      END 
