      SUBROUTINE  SGBCO( ABD, LDA, N, ML, MU, IPVT, RCOND, Z )
C
C         FACTORS A REAL BAND MATRIX BY GAUSSIAN ELIMINATION
C         AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     IF  RCOND  IS NOT NEEDED, SGBFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW SBGCO BY SGBSL.
C
C     INPUT:
C
C        ABD     REAL(LDA, N)
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
C                SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C                0 .LE. ML .LT. N .
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. MU .LT. N .
C                MORE EFFICIENT IF  ML .LE. MU .
C
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       REAL(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     BAND STORAGE
C
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
C           WILL SET UP THE INPUT.
C
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-MU)
C                      I2 = MIN0(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
C
C     EXAMPLE:  IF THE ORIGINAL MATRIX IS
C
C           11 12 13  0  0  0
C           21 22 23 24  0  0
C            0 32 33 34 35  0
C            0  0 43 44 45 46
C            0  0  0 54 55 56
C            0  0  0  0 65 66
C
C      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN
C
C            *  *  *  +  +  +  , * = NOT USED
C            *  * 13 24 35 46  , + = USED FOR PIVOTING
C            * 12 23 34 45 56
C           11 22 33 44 55 66
C           21 32 43 54 65  *
C
C
C     ROUTINES CALLED:  FROM LINPACK: SGBFA
C                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
C                       FROM FORTRAN: ABS, AMAX1, MAX0, MIN0, SIGN
C
      INTEGER  LDA, N, ML, MU, IPVT(*)
      REAL     ABD(LDA,*), Z(*)
      REAL     RCOND
C
      REAL     SDOT, EK, T, WK, WKM
      REAL     ANORM, S, SASUM, SM, YNORM
      INTEGER  IS, INFO, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM
C
C
C                       ** COMPUTE 1-NORM OF A
      ANORM = 0.0E0
      L = ML + 1
      IS = L + MU
      DO 10 J = 1, N
         ANORM = AMAX1(ANORM, SASUM(L,ABD(IS,J), 1))
         IF (IS .GT. ML + 1) IS = IS - 1
         IF (J .LE. MU) L = L + 1
         IF (J .GE. N - ML) L = L - 1
   10 CONTINUE
C                                               ** FACTOR
      CALL SGBFA(ABD, LDA, N, ML, MU, IPVT, INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C                     ** SOLVE TRANS(U)*W = E
      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE
C
      M = ML + MU + 1
      JU = 0
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK, -Z(K))
         IF (ABS(EK-Z(K)) .GT. ABS(ABD(M,K))) THEN
            S = ABS(ABD(M,K))/ABS(EK-Z(K))
            CALL SSCAL(N, S, Z, 1)
            EK = S*EK
         ENDIF
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (ABD(M,K) .NE. 0.0E0) THEN
            WK  = WK /ABD(M,K)
            WKM = WKM/ABD(M,K)
         ELSE
            WK  = 1.0E0
            WKM = 1.0E0
         ENDIF
         KP1 = K + 1
         JU = MIN0(MAX0(JU, MU+IPVT(K)), N)
         MM = M
         IF (KP1 .LE. JU) THEN
            DO 60 J = KP1, JU
               MM = MM - 1
               SM = SM + ABS(Z(J)+WKM*ABD(MM,J))
               Z(J) = Z(J) + WK*ABD(MM,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .LT. SM) THEN
               T = WKM - WK
               WK = WKM
               MM = M
               DO 70 J = KP1, JU
                  MM = MM - 1
                  Z(J) = Z(J) + T*ABD(MM,J)
   70          CONTINUE
            ENDIF
         ENDIF
         Z(K) = WK
  100 CONTINUE
C
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
C
C                         ** SOLVE TRANS(L)*Y = W
      DO 120 KB = 1, N
         K = N + 1 - KB
         LM = MIN0(ML, N-K)
         IF (K .LT. N) Z(K) = Z(K) + SDOT(LM, ABD(M+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0 / ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
         ENDIF
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
C
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
C
      YNORM = 1.0E0
C                         ** SOLVE L*V = Y
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         LM = MIN0(ML, N-K)
         IF (K .LT. N) CALL SAXPY(LM, T, ABD(M+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0 / ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
  140 CONTINUE
C
      S = 1.0E0/SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
      YNORM = S*YNORM
C                           ** SOLVE  U*Z = W
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .GT. ABS(ABD(M,K))) THEN
            S = ABS(ABD(M,K)) / ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
         IF (ABD(M,K) .NE. 0.0E0) Z(K) = Z(K)/ABD(M,K)
         IF (ABD(M,K) .EQ. 0.0E0) Z(K) = 1.0E0
         LM = MIN0(K, M) - 1
         LA = M - LM
         LZ = K - LM
         T = -Z(K)
         CALL SAXPY(LM, T, ABD(LA,K), 1, Z(LZ), 1)
  160 CONTINUE
C                              ** MAKE ZNORM = 1.0
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
      SUBROUTINE  SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )
C
C         FACTORS A REAL BAND MATRIX BY ELIMINATION.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     SGBFA IS USUALLY CALLED BY SBGCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C
C     INPUT:  SAME AS 'SGBCO'
C
C     ON RETURN:
C
C        ABD,IPVT    SAME AS 'SGBCO'
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGBSL WILL DIVIDE BY ZERO IF
C                     CALLED.  USE  RCOND  IN SBGCO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     (SEE 'SGBCO' FOR DESCRIPTION OF BAND STORAGE MODE)
C
C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
C                       FROM FORTRAN: MAX0, MIN0
C
      INTEGER  LDA, N, ML, MU, IPVT(*), INFO
      REAL     ABD(LDA,*)
C
      REAL     T
      INTEGER  I,ISAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
C
C
      M = ML + MU + 1
      INFO = 0
C                        ** ZERO INITIAL FILL-IN COLUMNS
      J0 = MU + 2
      J1 = MIN0(N, M) - 1
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0E0
   10    CONTINUE
   20 CONTINUE
      JZ = J1
      JU = 0
C
C                       ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
      NM1 = N - 1
      DO 120 K = 1, NM1
         KP1 = K + 1
C                                  ** ZERO NEXT FILL-IN COLUMN
         JZ = JZ + 1
         IF (JZ .LE. N) THEN
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0E0
   40       CONTINUE
         ENDIF
C                                  ** FIND L = PIVOT INDEX
         LM = MIN0(ML, N-K)
         L = ISAMAX(LM+1, ABD(M,K), 1) + M - 1
         IPVT(K) = L + K - M
C
         IF (ABD(L,K) .EQ. 0.0E0) THEN
C                                      ** ZERO PIVOT IMPLIES THIS COLUMN
C                                      ** ALREADY TRIANGULARIZED
            INFO = K
         ELSE
C                                ** INTERCHANGE IF NECESSARY
            IF (L .NE. M) THEN
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
            ENDIF
C                                   ** COMPUTE MULTIPLIERS
            T = -1.0E0 / ABD(M,K)
            CALL SSCAL(LM, T, ABD(M+1,K), 1)
C
C                               ** ROW ELIMINATION WITH COLUMN INDEXING
C
            JU = MIN0(MAX0(JU, MU+IPVT(K)), N)
            MM = M
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .NE. MM) THEN
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
               ENDIF
               CALL SAXPY(LM, T, ABD(M+1,K), 1, ABD(MM+1,J), 1)
   80       CONTINUE
C
         ENDIF
C
  120 CONTINUE
C
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0E0) INFO = N
      RETURN
      END
      SUBROUTINE  SGBSL( ABD, LDA, N, ML, MU, IPVT, B, JOB )
C
C         SOLVES THE REAL BAND SYSTEM
C            A * X = B  OR  TRANSPOSE(A) * X = B
C         USING THE FACTORS COMPUTED BY SBGCO OR SGBFA.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     INPUT:
C
C        ABD     REAL(LDA, N)
C                THE OUTPUT FROM SBGCO OR SGBFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SBGCO OR SGBFA.
C
C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SBGCO HAS SET RCOND .GT. 0.0
C        OR SGBFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT
C                       FROM FORTRAN: MIN0
C
      INTEGER  LDA, N, ML, MU, IPVT(*), JOB
      REAL     ABD(LDA,*), B(*)
C
      REAL     SDOT,T
      INTEGER  K,KB,L,LA,LB,LM,M,NM1
C
C
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .EQ. 0) THEN
C                               ** JOB = 0 , SOLVE  A * X = B
C                               ** FIRST SOLVE L*Y = B
         IF (ML .NE. 0) THEN
            DO 20 K = 1, NM1
               LM = MIN0(ML, N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .NE. K) THEN
                  B(L) = B(K)
                  B(K) = T
               ENDIF
               CALL SAXPY( LM, T, ABD(M+1,K), 1, B(K+1), 1 )
   20       CONTINUE
         ENDIF
C                           ** NOW SOLVE  U*X = Y
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K) / ABD(M,K)
            LM = MIN0(K, M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL SAXPY(LM, T, ABD(LA,K), 1, B(LB), 1)
   40    CONTINUE
C
      ELSE
C                          ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
C                                  ** FIRST SOLVE  TRANS(U)*Y = B
         DO 60 K = 1, N
            LM = MIN0(K, M) - 1
            LA = M - LM
            LB = K - LM
            T = SDOT(LM, ABD(LA,K), 1, B(LB), 1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
C                                  ** NOW SOLVE TRANS(L)*X = Y
         IF (ML .NE. 0) THEN
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN0(ML, N-K)
               B(K) = B(K) + SDOT(LM, ABD(M+1,K), 1, B(K+1), 1)
               L = IPVT(K)
               IF (L .NE. K) THEN
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
               ENDIF
   80       CONTINUE
         ENDIF
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE  SGECO( A, LDA, N,IPVT, RCOND, Z )
C
C         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
C         AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C         IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
C         TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       REAL(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     ROUTINES CALLED:  FROM LINPACK: SGEFA
C                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
C                       FROM FORTRAN: ABS, AMAX1, SIGN
C
      INTEGER  LDA, N, IPVT(*)
      REAL     A(LDA,*), Z(*)
      REAL     RCOND
C
      REAL     SDOT,EK,T,WK,WKM
      REAL     ANORM,S,SASUM,SM,YNORM
      INTEGER  INFO,J,K,KB,KP1,L
C
C
C                        ** COMPUTE 1-NORM OF A
      ANORM = 0.0E0
      DO 10 J = 1, N
         ANORM = AMAX1( ANORM, SASUM(N,A(1,J),1) )
   10 CONTINUE
C                                      ** FACTOR
      CALL SGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C                        ** SOLVE TRANS(U)*W = E
      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE
C
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK, -Z(K))
         IF (ABS(EK-Z(K)) .GT. ABS(A(K,K))) THEN
            S = ABS(A(K,K)) / ABS(EK-Z(K))
            CALL SSCAL(N, S, Z, 1)
            EK = S*EK
         ENDIF
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (A(K,K) .NE. 0.0E0) THEN
            WK  = WK  / A(K,K)
            WKM = WKM / A(K,K)
         ELSE
            WK  = 1.0E0
            WKM = 1.0E0
         ENDIF
         KP1 = K + 1
         IF (KP1 .LE. N) THEN
            DO 60 J = KP1, N
               SM = SM + ABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .LT. SM) THEN
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
            ENDIF
         ENDIF
         Z(K) = WK
  100 CONTINUE
C
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
C                                ** SOLVE TRANS(L)*Y = W
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + SDOT(N-K, A(K+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0/ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
         ENDIF
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
C
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
C                                 ** SOLVE L*V = Y
      YNORM = 1.0E0
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL SAXPY(N-K, T, A(K+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0/ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
  140 CONTINUE
C
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
C                                  ** SOLVE  U*Z = V
      YNORM = S*YNORM
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .GT. ABS(A(K,K))) THEN
            S = ABS(A(K,K))/ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
         IF (A(K,K) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0E0) Z(K) = 1.0E0
         T = -Z(K)
         CALL SAXPY(K-1, T, A(1,K), 1, Z(1), 1)
  160 CONTINUE
C                                   ** MAKE ZNORM = 1.0
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
      SUBROUTINE  SGEFA( A, LDA, N, IPVT, INFO )
C
C         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .
C
C     INPUT:  SAME AS 'SGECO'
C
C     ON RETURN:
C
C        A,IPVT  SAME AS 'SGECO'
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
C
      INTEGER  LDA, N, IPVT(*), INFO
      REAL     A(LDA,*)
C
      REAL     T
      INTEGER  ISAMAX,J,K,KP1,L,NM1
C
C
C                      ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
      INFO = 0
      NM1 = N - 1
      DO 60 K = 1, NM1
         KP1 = K + 1
C                                            ** FIND L = PIVOT INDEX
         L = ISAMAX( N-K+1, A(K,K), 1) + K-1
         IPVT(K) = L
C
         IF (A(L,K) .EQ. 0.0E0) THEN
C                                     ** ZERO PIVOT IMPLIES THIS COLUMN
C                                     ** ALREADY TRIANGULARIZED
            INFO = K
         ELSE
C                                     ** INTERCHANGE IF NECESSARY
            IF (L .NE. K) THEN
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
            ENDIF
C                                     ** COMPUTE MULTIPLIERS
            T = -1.0E0 / A(K,K)
            CALL SSCAL( N-K, T, A(K+1,K), 1 )
C
C                              ** ROW ELIMINATION WITH COLUMN INDEXING
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .NE. K) THEN
                  A(L,J) = A(K,J)
                  A(K,J) = T
               ENDIF
               CALL SAXPY( N-K, T, A(K+1,K), 1, A(K+1,J), 1 )
   30       CONTINUE
C
         ENDIF
C
   60 CONTINUE
C
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0E0) INFO = N
      RETURN
      END
      SUBROUTINE  SGESL( A, LDA, N,IPVT, B, JOB )
C
C         SOLVES THE REAL SYSTEM
C            A * X = B  OR  TRANS(A) * X = B
C         USING THE FACTORS COMPUTED BY SGECO OR SGEFA.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE OUTPUT FROM SGECO OR SGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SGECO OR SGEFA.
C
C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
C        OR SGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C
C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT
C
      INTEGER  LDA, N, IPVT(*), JOB
      REAL     A(LDA,*), B(*)
C
      REAL     SDOT,T
      INTEGER  K,KB,L,NM1
C
C
      NM1 = N - 1
      IF (JOB .EQ. 0) THEN
C                                 ** JOB = 0 , SOLVE  A * X = B
C                                     ** FIRST SOLVE  L*Y = B
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .NE. K) THEN
               B(L) = B(K)
               B(K) = T
            ENDIF
            CALL SAXPY( N-K, T, A(K+1,K), 1, B(K+1), 1 )
   20    CONTINUE
C                                    ** NOW SOLVE  U*X = Y
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K) / A(K,K)
            T = -B(K)
            CALL SAXPY( K-1, T, A(1,K), 1, B(1), 1 )
   40    CONTINUE
C
      ELSE
C                         ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
C                                    ** FIRST SOLVE  TRANS(U)*Y = B
         DO 60 K = 1, N
            T = SDOT( K-1, A(1,K), 1, B(1), 1 )
            B(K) = (B(K) - T) / A(K,K)
   60    CONTINUE
C                                    ** NOW SOLVE  TRANS(L)*X = Y
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + SDOT( N-K, A(K+1,K), 1, B(K+1), 1 )
            L = IPVT(K)
            IF (L .NE. K) THEN
               T = B(L)
               B(L) = B(K)
               B(K) = T
            ENDIF
   80    CONTINUE
C
      ENDIF
C
      RETURN
      END
      REAL FUNCTION  SASUM( N, SX, INCX )
C
C  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR TO BE SUMMED
C            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
C          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'
C
C --OUTPUT-- SASUM   SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))
C
      REAL SX(*)
C
C
      SASUM = 0.0
      IF( N.LE.0 )  RETURN
      IF( INCX.NE.1 ) THEN
C                                          ** NON-UNIT INCREMENTS
          DO 10 I = 1, 1+(N-1)*INCX, INCX
             SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
      ELSE
C                                          ** UNIT INCREMENTS
         M = MOD(N,6)
         IF( M.NE.0 ) THEN
C                             ** CLEAN-UP LOOP SO REMAINING VECTOR
C                             ** LENGTH IS A MULTIPLE OF 6.
            DO 30  I = 1, M
              SASUM = SASUM + ABS(SX(I))
   30       CONTINUE
         ENDIF
C                              ** UNROLL LOOP FOR SPEED
         DO 50  I = M+1, N, 6
           SASUM = SASUM + ABS(SX(I))   + ABS(SX(I+1)) + ABS(SX(I+2))
     $                   + ABS(SX(I+3)) + ABS(SX(I+4)) + ABS(SX(I+5))
   50    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE     SAXPY( N, SA, SX, INCX, SY, INCY )
C
C          Y = A*X + Y  (X, Y = VECTORS, A = SCALAR)
C
C  --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
C       SA  SINGLE PRECISION SCALAR MULTIPLIER 'A'
C       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
C     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
C       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
C     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'
C
C --OUTPUT--
C       SY   FOR I = 0 TO N-1, OVERWRITE  SY(LY+I*INCY) WITH
C                 SA*SX(LX+I*INCX) + SY(LY+I*INCY),
C            WHERE LX = 1          IF INCX .GE. 0,
C                     = (-INCX)*N  IF INCX .LT. 0
C            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      REAL SX(*), SY(*), SA
C
C
      IF( N.LE.0 .OR. SA.EQ.0.0 ) RETURN
C
      IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN
C
          DO 10  I = 1, 1+(N-1)*INCX, INCX
             SY(I) = SY(I) + SA * SX(I)
   10     CONTINUE
C
      ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN
C
C                                        ** EQUAL, UNIT INCREMENTS
         M = MOD(N,4)
         IF( M .NE. 0 ) THEN
C                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                            ** IS A MULTIPLE OF 4.
            DO 20  I = 1, M
              SY(I) = SY(I) + SA * SX(I)
   20       CONTINUE
         ENDIF
C                              ** UNROLL LOOP FOR SPEED
         DO 30  I = M+1, N, 4
            SY(I)   = SY(I)   + SA * SX(I)
            SY(I+1) = SY(I+1) + SA * SX(I+1)
            SY(I+2) = SY(I+2) + SA * SX(I+2)
            SY(I+3) = SY(I+3) + SA * SX(I+3)
   30    CONTINUE
C
      ELSE
C               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
         IX = 1
         IY = 1
         IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
         IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
         DO 40  I = 1, N
            SY(IY) = SY(IY) + SA*SX(IX)
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE
C
      ENDIF
C
      RETURN
      END
      REAL FUNCTION  SDOT( N, SX, INCX, SY, INCY )
C
C          S.P. DOT PRODUCT OF VECTORS  'X'  AND  'Y'
C
C  --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
C       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
C     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
C       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
C     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'
C
C --OUTPUT--
C     SDOT   SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
C            WHERE  LX = 1          IF INCX .GE. 0,
C                      = (-INCX)*N  IF INCX .LT. 0,
C            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      REAL SX(*), SY(*)
C
C
      SDOT = 0.0
      IF( N.LE.0 )  RETURN
C
      IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN
C
          DO 10  I = 1, 1+(N-1)*INCX, INCX
             SDOT = SDOT + SX(I) * SY(I)
   10     CONTINUE
C
      ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN
C
C                                        ** EQUAL, UNIT INCREMENTS
         M = MOD(N,5)
         IF( M .NE. 0 ) THEN
C                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                            ** IS A MULTIPLE OF 4.
            DO 20  I = 1, M
               SDOT = SDOT + SX(I) * SY(I)
   20       CONTINUE
         ENDIF
C                              ** UNROLL LOOP FOR SPEED
         DO 30  I = M+1, N, 5
            SDOT = SDOT + SX(I)*SY(I)     + SX(I+1)*SY(I+1)
     $                  + SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3)
     $                  + SX(I+4)*SY(I+4)
   30    CONTINUE
C
      ELSE
C               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
         IX = 1
         IY = 1
         IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
         IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
         DO 40  I = 1, N
            SDOT = SDOT + SX(IX) * SY(IY)
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE
C
      ENDIF
C
      RETURN
      END
      LOGICAL FUNCTION  WRTBADQ ( quiet, VARNAM )
C
C          WRITE NAMES OF ERRONEOUS VARIABLES AND RETURN 'TRUE'
C
C      INPUT :   VARNAM = NAME OF ERRONEOUS VARIABLE TO BE WRITTEN
C                         ( CHARACTER, ANY LENGTH )
C ----------------------------------------------------------------------
      CHARACTER*(*)  VARNAM
      INTEGER        MAXMSG, NUMMSG
      LOGICAL quiet
      SAVE  NUMMSG, MAXMSG
      DATA  NUMMSG / 0 /,  MAXMSG / 50 /
C
C
      WRTBADQ = .TRUE.
      NUMMSG = NUMMSG + 1 
      WRITE ( *, '(3A)' )  ' ****  INPUT VARIABLE  ', VARNAM,
     $     '  IN ERROR  ****'
      IF ( NUMMSG.EQ.MAXMSG  .AND. .NOT. quiet )
     $ CALL  ERRMSG( 'TOO MANY INPUT ERRORS.  ABORTING...$', .TRUE. )
      RETURN
      END
      LOGICAL FUNCTION  WRTDIMQ ( quiet, DIMNAM, MINVAL )
C
C          WRITE NAME OF TOO-SMALL SYMBOLIC DIMENSION AND
C          THE VALUE IT SHOULD BE INCREASED TO;  RETURN 'TRUE'
C
C      INPUT :  DIMNAM = NAME OF SYMBOLIC DIMENSION WHICH IS TOO SMALL
C                        ( CHARACTER, ANY LENGTH )
C               MINVAL = VALUE TO WHICH THAT DIMENSION SHOULD BE
C                        INCREASED (AT LEAST)
C ----------------------------------------------------------------------
      CHARACTER*(*)  DIMNAM
      INTEGER        MINVAL
      LOGICAL quiet
C
C
      WRITE ( *, '(3A,I7)' )  ' ****  SYMBOLIC DIMENSION  ',
     $     DIMNAM, '  SHOULD BE INCREASED TO AT LEAST ', MINVAL
      WRTDIMQ = .TRUE.
      RETURN
      END

