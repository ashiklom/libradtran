c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     LINEAR-EQUATION-SOLVING LAPACK (v3.0) ROUTINES FOR DISORT v1.3
c     (current as of Nov 99)
c
c NOTES:  
c
c (1) If possible, use locally optimized versions of these routines
c     instead of these Fortran versions.  A significant portion 
c     of DISORT computer time is spent in these routines.
c
c (2) Prebuilt LAPACK libraries are available for a variety of
c     computers.  These will run much faster.  See:
c
c        http://www.netlib.org/lapack/archives/
c
c (3) Upgrades to LAPACK are available from 
c
c        http://www.netlib.org/lapack/
c        http://netlib.bell-labs.com/netlib/master/readme.html
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c CHANGES MADE TO ORIGINAL LAPACK:
c
c (1) Major simplifications in routine ILAENV.
c
c (2) Removal of the following comment block from the beginning of
c     each routine:
c  -- LAPACK routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c
c (3) Other cosmetic edits.  
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c Call tree (omitting calls to error msg routine XERBLA):
c
c    SGBSV
c        SGBTRF
c            ILAENV
c            ISAMAX [BLAS-1]
c            SCOPY [BLAS-1]
c            SGBTF2
c                ISAMAX [BLAS-1]
c                SGER [BLAS-2]
c                SSCAL [BLAS-1]
c                SSWAP [BLAS-1]
c            SGEMM [BLAS-3]
c                LSAME
c            SGER [BLAS-2]
c            SLASWP
c                SSWAP [BLAS-1]
c            SSCAL [BLAS-1]
c            SSWAP [BLAS-1]
c            STRSM [BLAS-3]
c                LSAME
c        SGBTRS
c            LSAME [BLAS-2]
c            SGEMV [BLAS-2]
c                LSAME
c            SGER [BLAS-2]
c            SSWAP [BLAS-1]
c            STBSV [BLAS-2]
c                LSAME
c    SGESV
c        SGETRF
c            SGEMM [BLAS-3]
c                LSAME
c            SGETF2
c                ISAMAX [BLAS-1]
c                SGER [BLAS-2]
c                SSCAL [BLAS-1]
c                SSWAP [BLAS-1]
c            SLASWP
c                SSWAP [BLAS-1]
c            STRSM [BLAS-3]
c            ILAENV
c        SGETRS
c            LSAME [BLAS-2]
c            SLASWP
c                SSWAP [BLAS-1]
c            STRSM [BLAS-3]
c                LSAME

c  The BLAS-1,2,3 routines are in a separate file.  See if you have
c  optimized local versions of BLAS routines before using that file.
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      SUBROUTINE SGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
c
c  -- LAPACK driver routine (version 2.0) --
c     March 31, 1993 
c
c     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               AB( LDAB, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  SGBSV computes the solution to a real system of linear equations
c  A * X = B, where A is a band matrix of order N with KL subdiagonals
c  and KU superdiagonals, and X and B are N-by-NRHS matrices.
c
c  The LU decomposition with partial pivoting and row interchanges is
c  used to factor A as A = L * U, where L is a product of permutation
c  and unit lower triangular matrices with KL subdiagonals, and U is
c  upper triangular with KL+KU superdiagonals.  The factored form of A
c  is then used to solve the system of equations A * X = B.
c
c  Arguments
c  =========
c
c  N       (input) INTEGER
c          The number of linear equations, i.e., the order of the
c          matrix A.  N >= 0.
c
c  KL      (input) INTEGER
c          The number of subdiagonals within the band of A.  KL >= 0.
c
c  KU      (input) INTEGER
c          The number of superdiagonals within the band of A.  KU >= 0.
c
c  NRHS    (input) INTEGER
c          The number of right hand sides, i.e., the number of columns
c          of the matrix B.  NRHS >= 0.
c
c  AB      (input/output) REAL array, dimension (LDAB,N)
c          On entry, the matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
c          On exit, details of the factorization: U is stored as an
c          upper triangular band matrix with KL+KU superdiagonals in
c          rows 1 to KL+KU+1, and the multipliers used during the
c          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
c          See below for further details.
c
c  LDAB    (input) INTEGER
c          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
c
c  IPIV    (output) INTEGER array, dimension (N)
c          The pivot indices that define the permutation matrix P;
c          row i of the matrix was interchanged with row IPIV(i).
c
c  B       (input/output) REAL array, dimension (LDB,NRHS)
c          On entry, the N-by-NRHS right hand side matrix B.
c          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
c
c  LDB     (input) INTEGER
c          The leading dimension of the array B.  LDB >= max(1,N).
c
c  INFO    (output) INTEGER
c          = 0:  successful exit
c          < 0:  if INFO = -i, the i-th argument had an illegal value
c          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
c                has been completed, but the factor U is exactly
c                singular, and the solution has not been computed.
c
c  Further Details
c  ===============
c
c  The band storage scheme is illustrated by the following example, when
c  M = N = 6, KL = 2, KU = 1:
c
c  On entry:                       On exit:
c
c      *    *    *    +    +    +       *    *    *   u14  u25  u36
c      *    *    +    +    +    +       *    *   u13  u24  u35  u46
c      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
c     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
c     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
c     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
c
c  Array elements marked * are not used by the routine; elements marked
c  + need not be set on entry, but are required by the routine to store
c  elements of U because of fill-in resulting from the row interchanges.
c
c  =====================================================================
c
c     .. External Subroutines ..
      EXTERNAL           SGBTRF, SGBTRS, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..

c     Test the input parameters.

      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( KL.LT.0 ) THEN
         INFO = -2
      ELSE IF( KU.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.2*KL+KU+1 ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( N, 1 ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBSV ', -INFO )
         RETURN
      END IF

c     Compute the LU factorization of the band matrix A.

      CALL SGBTRF( N, N, KL, KU, AB, LDAB, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN

c        Solve the system A*X = B, overwriting B with X.

         CALL SGBTRS( 'No transpose', N, KL, KU, NRHS, AB, LDAB, IPIV,
     &                B, LDB, INFO )
      END IF
      RETURN

c          End of SGBSV
      END

      SUBROUTINE SGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
c
c     February 29, 1992
c
c     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               AB( LDAB, * )
c     ..
c
c  Purpose
c  =======
c
c  SGBTRF computes an LU factorization of a real m-by-n band matrix A
c  using partial pivoting with row interchanges.
c
c  This is the blocked version of the algorithm, calling Level 3 BLAS.
c
c  Arguments
c  =========
c
c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.
c
c  KL      (input) INTEGER
c          The number of subdiagonals within the band of A.  KL >= 0.
c
c  KU      (input) INTEGER
c          The number of superdiagonals within the band of A.  KU >= 0.
c
c  AB      (input/output) REAL array, dimension (LDAB,N)
c          On entry, the matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c
c          On exit, details of the factorization: U is stored as an
c          upper triangular band matrix with KL+KU superdiagonals in
c          rows 1 to KL+KU+1, and the multipliers used during the
c          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
c          See below for further details.
c
c  LDAB    (input) INTEGER
c          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
c
c  IPIV    (output) INTEGER array, dimension (min(M,N))
c          The pivot indices; for 1 <= i <= min(M,N), row i of the
c          matrix was interchanged with row IPIV(i).
c
c  INFO    (output) INTEGER
c          = 0: successful exit
c          < 0: if INFO = -i, the i-th argument had an illegal value
c          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
c               has been completed, but the factor U is exactly
c               singular, and division by zero will occur if it is used
c               to solve a system of equations.
c
c  Further Details
c  ===============
c
c  The band storage scheme is illustrated by the following example, when
c  M = N = 6, KL = 2, KU = 1:
c
c  On entry:                       On exit:
c
c      *    *    *    +    +    +       *    *    *   u14  u25  u36
c      *    *    +    +    +    +       *    *   u13  u24  u35  u46
c      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
c     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
c     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
c     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
c
c  Array elements marked * are not used by the routine; elements marked
c  + need not be set on entry, but are required by the routine to store
c  elements of U because of fill-in resulting from the row interchanges.
c
c  =====================================================================
c
c     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP,
     &                   JU, K2, KM, KV, NB, NW
      REAL               TEMP
c     ..
c     .. Local Arrays ..
      REAL               WORK13( LDWORK, NBMAX ),
     &                   WORK31( LDWORK, NBMAX )
c     ..
c     .. External Functions ..
      INTEGER            ILAENV, ISAMAX
      EXTERNAL           ILAENV, ISAMAX
c     ..
c     .. External Subroutines ..
      EXTERNAL           SCOPY, SGBTF2, SGEMM, SGER, SLASWP, SSCAL,
     &                   SSWAP, STRSM, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..


c     KV is the number of superdiagonals in the factor U, allowing for
c     fill-in

      KV = KU + KL

c     Test the input parameters.

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBTRF', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 )  RETURN

c     Determine the block size for this environment

      NB = ILAENV( 1, 'SGBTRF', ' ', M, N, KL, KU )

c     The block size must not exceed the limit set by the size of the
c     local arrays WORK13 and WORK31.

      NB = MIN( NB, NBMAX )

      IF( NB.LE.1 .OR. NB.GT.KL ) THEN

c        Use unblocked code

         CALL SGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE

c        Use blocked code

c        Zero the superdiagonal elements of the work array WORK13

         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE

c        Zero the subdiagonal elements of the work array WORK31

         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE

c        Gaussian elimination with partial pivoting

c        Set fill-in elements in columns KU+2 to KV to zero

         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE

c        JU is the index of the last column affected by the current
c        stage of the factorization

         JU = 1

         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )

c           The active part of the matrix is partitioned
c
c              A11   A12   A13
c              A21   A22   A23
c              A31   A32   A33
c
c           Here A11, A21 and A31 denote the current block of JB columns
c           which is about to be factorized. The number of rows in the
c           partitioning are JB, I2, I3 respectively, and the numbers
c           of columns are JB, J2, J3. The superdiagonal elements of A13
c           and the subdiagonal elements of A31 lie outside the band.

            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )

c           J2 and J3 are computed after JU has been updated.

c           Factorize the current block of JB columns

            DO 80 JJ = J, J + JB - 1

c              Set fill-in elements in column JJ+KV to zero

               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF

c              Find pivot and test for singularity. KM is the number of
c              subdiagonal elements in the current column.

               KM = MIN( KL, M-JJ )
               JP = ISAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN

c                    Apply interchange to columns J to J+JB-1

                     IF( JP+JJ-1.LT.J+KL ) THEN

                        CALL SSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1,
     &                              AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE

c                       The interchange affects columns J to JJ-1 of A31
c                       which are stored in the work array WORK31

                        CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     &                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL SSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1,
     &                              AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF

c                 Compute multipliers

                  CALL SSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ),
     &                        1 )

c                 Update trailing submatrix within the band and within
c                 the current block. JM is the index of the last column
c                 which needs to be updated.

                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ )
     &               CALL SGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1,
     &                          AB( KV, JJ+1 ), LDAB-1,
     &                          AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE

c                 If pivot is zero, set INFO to the index of the pivot
c                 unless a zero pivot has already been found.

                  IF( INFO.EQ.0 )
     &               INFO = JJ
               END IF

c              Copy current column of A31 into the work array WORK31

               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 )
     &            CALL SCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1,
     &                        WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN

c              Apply the row interchanges to the other blocks.

               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )

c              Use SLASWP to apply the row interchanges to A12, A22, and
c              A32.

               CALL SLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB,
     &                      IPIV( J ), 1 )

c              Adjust the pivot indices.

               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE

c              Apply the row interchanges to A13, A23, and A33
c              columnwise.

               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE

c              Update the relevant part of the trailing submatrix

               IF( J2.GT.0 ) THEN

c                 Update A12

                  CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     &                        JB, J2, ONE, AB( KV+1, J ), LDAB-1,
     &                        AB( KV+1-JB, J+JB ), LDAB-1 )

                  IF( I2.GT.0 ) THEN

c                    Update A22

                     CALL SGEMM( 'No transpose', 'No transpose', I2, J2,
     &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     &                           AB( KV+1, J+JB ), LDAB-1 )
                  END IF

                  IF( I3.GT.0 ) THEN

c                    Update A32

                     CALL SGEMM( 'No transpose', 'No transpose', I3, J2,
     &                           JB, -ONE, WORK31, LDWORK,
     &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     &                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF

               IF( J3.GT.0 ) THEN

c                 Copy the lower triangle of A13 into the work array
c                 WORK13

                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE

c                 Update A13 in the work array

                  CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     &                        JB, J3, ONE, AB( KV+1, J ), LDAB-1,
     &                        WORK13, LDWORK )

                  IF( I2.GT.0 ) THEN

c                    Update A23

                     CALL SGEMM( 'No transpose', 'No transpose', I2, J3,
     &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     &                           WORK13, LDWORK, ONE, AB( 1+JB, J+KV ),
     &                           LDAB-1 )
                  END IF

                  IF( I3.GT.0 ) THEN

c                    Update A33

                     CALL SGEMM( 'No transpose', 'No transpose', I3, J3,
     &                           JB, -ONE, WORK31, LDWORK, WORK13,
     &                           LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF

c                 Copy the lower triangle of A13 back into place

                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE

c              Adjust the pivot indices.

               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF

c           Partially undo the interchanges in the current block to
c           restore the upper triangular form of A31 and copy the upper
c           triangle of A31 back into place

            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN

c                 Apply interchange to columns J to JJ-1

                  IF( JP+JJ-1.LT.J+KL ) THEN

c                    The interchange does not affect A31

                     CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     &                           AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE

c                    The interchange does affect A31

                     CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     &                           WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF

c              Copy the current column of A31 back into place

               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 )
     &            CALL SCOPY( NW, WORK31( 1, JJ-J+1 ), 1,
     &                        AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF

      RETURN

c            End of SGBTRF
      END


      SUBROUTINE SGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
     &                   INFO )
c
c     March 31, 1993 
c
c     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               AB( LDAB, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  SGBTRS solves a system of linear equations
c     A * X = B  or  A' * X = B
c  with a general band matrix A using the LU factorization computed
c  by SGBTRF.
c
c  Arguments
c  =========
c
c  TRANS   (input) CHARACTER*1
c          Specifies the form of the system of equations.
c          = 'N':  A * X = B  (No transpose)
c          = 'T':  A'* X = B  (Transpose)
c          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
c
c  N       (input) INTEGER
c          The order of the matrix A.  N >= 0.
c
c  KL      (input) INTEGER
c          The number of subdiagonals within the band of A.  KL >= 0.
c
c  KU      (input) INTEGER
c          The number of superdiagonals within the band of A.  KU >= 0.
c
c  NRHS    (input) INTEGER
c          The number of right hand sides, i.e., the number of columns
c          of the matrix B.  NRHS >= 0.
c
c  AB      (input) REAL array, dimension (LDAB,N)
c          Details of the LU factorization of the band matrix A, as
c          computed by SGBTRF.  U is stored as an upper triangular band
c          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
c          the multipliers used during the factorization are stored in
c          rows KL+KU+2 to 2*KL+KU+1.
c
c  LDAB    (input) INTEGER
c          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
c
c  IPIV    (input) INTEGER array, dimension (N)
c          The pivot indices; for 1 <= i <= N, row i of the matrix was
c          interchanged with row IPIV(i).
c
c  B       (input/output) REAL array, dimension (LDB,NRHS)
c          On entry, the right hand side matrix B.
c          On exit, the solution matrix X.
c
c  LDB     (input) INTEGER
c          The leading dimension of the array B.  LDB >= max(1,N).
c
c  INFO    (output) INTEGER
c          = 0:  successful exit
c          < 0: if INFO = -i, the i-th argument had an illegal value
c
c  =====================================================================
c
c     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            LNOTI, NOTRAN
      INTEGER            I, J, KD, L, LM
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           SGEMV, SGER, SSWAP, STBSV, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..


c     Test the input parameters.

      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBTRS', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 )  RETURN

      KD = KU + KL + 1
      LNOTI = KL.GT.0

      IF( NOTRAN ) THEN

c        Solve  A*X = B.
c
c        Solve L*X = B, overwriting B with X.
c
c        L is represented as a product of permutations and unit lower
c        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
c        where each transformation L(i) is a rank-one modification of
c        the identity matrix.

         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J )
     &            CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
               CALL SGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ),
     &                    LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF

         DO 20 I = 1, NRHS

c           Solve U*X = B, overwriting B with X.

            CALL STBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU,
     &                  AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE

      ELSE

c        Solve A'*X = B.

         DO 30 I = 1, NRHS

c           Solve U'*X = B, overwriting B with X.

            CALL STBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB,
     &                  LDAB, B( 1, I ), 1 )
   30    CONTINUE

c        Solve L'*X = B, overwriting B with X.

         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL SGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ),
     &                     LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J )
     &            CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
      END IF
      RETURN

c           End of SGBTRS
      END


      SUBROUTINE SGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
c
c     February 29, 1992
c
c     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               AB( LDAB, * )
c     ..
c
c  Purpose
c  =======
c
c  SGBTF2 computes an LU factorization of a real m-by-n band matrix A
c  using partial pivoting with row interchanges.
c
c  This is the unblocked version of the algorithm, calling Level 2 BLAS.
c
c  Arguments
c  =========
c
c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.
c
c  KL      (input) INTEGER
c          The number of subdiagonals within the band of A.  KL >= 0.
c
c  KU      (input) INTEGER
c          The number of superdiagonals within the band of A.  KU >= 0.
c
c  AB      (input/output) REAL array, dimension (LDAB,N)
c          On entry, the matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c
c          On exit, details of the factorization: U is stored as an
c          upper triangular band matrix with KL+KU superdiagonals in
c          rows 1 to KL+KU+1, and the multipliers used during the
c          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
c          See below for further details.
c
c  LDAB    (input) INTEGER
c          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
c
c  IPIV    (output) INTEGER array, dimension (min(M,N))
c          The pivot indices; for 1 <= i <= min(M,N), row i of the
c          matrix was interchanged with row IPIV(i).
c
c  INFO    (output) INTEGER
c          = 0: successful exit
c          < 0: if INFO = -i, the i-th argument had an illegal value
c          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
c               has been completed, but the factor U is exactly
c               singular, and division by zero will occur if it is used
c               to solve a system of equations.
c
c  Further Details
c  ===============
c
c  The band storage scheme is illustrated by the following example, when
c  M = N = 6, KL = 2, KU = 1:
c
c  On entry:                       On exit:
c
c      *    *    *    +    +    +       *    *    *   u14  u25  u36
c      *    *    +    +    +    +       *    *   u13  u24  u35  u46
c      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
c     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
c     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
c     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
c
c  Array elements marked * are not used by the routine; elements marked
c  + need not be set on entry, but are required by the routine to store
c  elements of U, because of fill-in resulting from the row
c  interchanges.
c
c  =====================================================================
c
c     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
c     ..
c     .. External Functions ..
      INTEGER            ISAMAX
      EXTERNAL           ISAMAX
c     ..
c     .. External Subroutines ..
      EXTERNAL           SGER, SSCAL, SSWAP, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..


c     KV is the number of superdiagonals in the factor U, allowing for
c     fill-in.

      KV = KU + KL

c     Test the input parameters.

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBTF2', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 )  RETURN

c     Gaussian elimination with partial pivoting

c     Set fill-in elements in columns KU+2 to KV to zero.

      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE

c     JU is the index of the last column affected by the current stage
c     of the factorization.

      JU = 1

      DO 40 J = 1, MIN( M, N )

c        Set fill-in elements in column J+KV to zero.

         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF

c        Find pivot and test for singularity. KM is the number of
c        subdiagonal elements in the current column.

         KM = MIN( KL, M-J )
         JP = ISAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )

c           Apply interchange to columns J to JU.

            IF( JP.NE.1 )
     &         CALL SSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1,
     &                     AB( KV+1, J ), LDAB-1 )

            IF( KM.GT.0 ) THEN

c              Compute multipliers.

               CALL SSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )

c              Update trailing submatrix within the band.

               IF( JU.GT.J )
     &            CALL SGER( KM, JU-J, -ONE, AB( KV+2, J ), 1,
     &                       AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ),
     &                       LDAB-1 )
            END IF
         ELSE

c           If pivot is zero, set INFO to the index of the pivot
c           unless a zero pivot has already been found.

            IF( INFO.EQ.0 )
     &         INFO = J
         END IF
   40 CONTINUE
      RETURN

c           End of SGBTF2
      END


      SUBROUTINE SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
c
c  -- LAPACK driver routine (version 2.0) --
c     March 31, 1993 
c
c     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, N, NRHS
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  SGESV computes the solution to a real system of linear equations
c     A * X = B,
c  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
c
c  The LU decomposition with partial pivoting and row interchanges is
c  used to factor A as
c     A = P * L * U,
c  where P is a permutation matrix, L is unit lower triangular, and U is
c  upper triangular.  The factored form of A is then used to solve the
c  system of equations A * X = B.
c
c  Arguments
c  =========
c
c  N       (input) INTEGER
c          The number of linear equations, i.e., the order of the
c          matrix A.  N >= 0.
c
c  NRHS    (input) INTEGER
c          The number of right hand sides, i.e., the number of columns
c          of the matrix B.  NRHS >= 0.
c
c  A       (input/output) REAL array, dimension (LDA,N)
c          On entry, the N-by-N coefficient matrix A.
c          On exit, the factors L and U from the factorization
c          A = P*L*U; the unit diagonal elements of L are not stored.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,N).
c
c  IPIV    (output) INTEGER array, dimension (N)
c          The pivot indices that define the permutation matrix P;
c          row i of the matrix was interchanged with row IPIV(i).
c
c  B       (input/output) REAL array, dimension (LDB,NRHS)
c          On entry, the N-by-NRHS matrix of right hand side matrix B.
c          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
c
c  LDB     (input) INTEGER
c          The leading dimension of the array B.  LDB >= max(1,N).
c
c  INFO    (output) INTEGER
c          = 0:  successful exit
c          < 0:  if INFO = -i, the i-th argument had an illegal value
c          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
c                has been completed, but the factor U is exactly
c                singular, so the solution could not be computed.
c
c  =====================================================================
c
c     .. External Subroutines ..
      EXTERNAL           SGETRF, SGETRS, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..


c     Test the input parameters.

      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGESV ', -INFO )
         RETURN
      END IF

c     Compute the LU factorization of A.

      CALL SGETRF( N, N, A, LDA, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN

c        Solve the system A*X = B, overwriting B with X.

         CALL SGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,
     &                INFO )
      END IF
      RETURN

c           End of SGESV
      END


      SUBROUTINE SGETRF( M, N, A, LDA, IPIV, INFO )
c
c     March 31, 1993 
c
c     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
c     ..
c
c  Purpose
c  =======
c
c  SGETRF computes an LU factorization of a general M-by-N matrix A
c  using partial pivoting with row interchanges.
c
c  The factorization has the form
c     A = P * L * U
c  where P is a permutation matrix, L is lower triangular with unit
c  diagonal elements (lower trapezoidal if m > n), and U is upper
c  triangular (upper trapezoidal if m < n).
c
c  This is the right-looking Level 3 BLAS version of the algorithm.
c
c  Arguments
c  =========
c
c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.
c
c  A       (input/output) REAL array, dimension (LDA,N)
c          On entry, the M-by-N matrix to be factored.
c          On exit, the factors L and U from the factorization
c          A = P*L*U; the unit diagonal elements of L are not stored.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,M).
c
c  IPIV    (output) INTEGER array, dimension (min(M,N))
c          The pivot indices; for 1 <= i <= min(M,N), row i of the
c          matrix was interchanged with row IPIV(i).
c
c  INFO    (output) INTEGER
c          = 0:  successful exit
c          < 0:  if INFO = -i, the i-th argument had an illegal value
c          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
c                has been completed, but the factor U is exactly
c                singular, and division by zero will occur if it is used
c                to solve a system of equations.
c
c  =====================================================================
c
c     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
c     ..
c     .. External Subroutines ..
      EXTERNAL           SGEMM, SGETF2, SLASWP, STRSM, XERBLA
c     ..
c     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..


c     Test the input parameters.

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGETRF', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 )  RETURN

c     Determine the block size for this environment.

      NB = ILAENV( 1, 'SGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN

c        Use unblocked code.

         CALL SGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE

c        Use blocked code.

         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )

c           Factor diagonal and subdiagonal blocks and test for exact
c           singularity.

            CALL SGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )

c           Adjust INFO and the pivot indices.

            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     &         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE

c           Apply interchanges to columns 1:J-1.

            CALL SLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )

            IF( J+JB.LE.N ) THEN

c              Apply interchanges to columns J+JB:N.

               CALL SLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     &                      IPIV, 1 )

c              Compute block row of U.

               CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     &                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     &                     LDA )
               IF( J+JB.LE.M ) THEN

c                 Update trailing submatrix.

                  CALL SGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     &                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     &                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     &                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN

c           End of SGETRF
      END


      SUBROUTINE SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
c
c     March 31, 1993 
c
c     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  SGETRS solves a system of linear equations
c     A * X = B  or  A' * X = B
c  with a general N-by-N matrix A using the LU factorization computed
c  by SGETRF.
c
c  Arguments
c  =========
c
c  TRANS   (input) CHARACTER*1
c          Specifies the form of the system of equations:
c          = 'N':  A * X = B  (No transpose)
c          = 'T':  A'* X = B  (Transpose)
c          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
c
c  N       (input) INTEGER
c          The order of the matrix A.  N >= 0.
c
c  NRHS    (input) INTEGER
c          The number of right hand sides, i.e., the number of columns
c          of the matrix B.  NRHS >= 0.
c
c  A       (input) REAL array, dimension (LDA,N)
c          The factors L and U from the factorization A = P*L*U
c          as computed by SGETRF.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,N).
c
c  IPIV    (input) INTEGER array, dimension (N)
c          The pivot indices from SGETRF; for 1<=i<=N, row i of the
c          matrix was interchanged with row IPIV(i).
c
c  B       (input/output) REAL array, dimension (LDB,NRHS)
c          On entry, the right hand side matrix B.
c          On exit, the solution matrix X.
c
c  LDB     (input) INTEGER
c          The leading dimension of the array B.  LDB >= max(1,N).
c
c  INFO    (output) INTEGER
c          = 0:  successful exit
c          < 0:  if INFO = -i, the i-th argument had an illegal value
c
c  =====================================================================
c
c     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            NOTRAN
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           SLASWP, STRSM, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..


c     Test the input parameters.

      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGETRS', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 )  RETURN

      IF( NOTRAN ) THEN

c        Solve A * X = B.
c
c        Apply row interchanges to the right hand sides.

         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )

c        Solve L*X = B, overwriting B with X.

         CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     &               ONE, A, LDA, B, LDB )

c        Solve U*X = B, overwriting B with X.

         CALL STRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     &               NRHS, ONE, A, LDA, B, LDB )
      ELSE

c        Solve A' * X = B.
c
c        Solve U'*X = B, overwriting B with X.

         CALL STRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     &               ONE, A, LDA, B, LDB )

c        Solve L'*X = B, overwriting B with X.

         CALL STRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     &               A, LDA, B, LDB )

c        Apply row interchanges to the solution vectors.

         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF

      RETURN

c           End of SGETRS
      END


      SUBROUTINE SGETF2( M, N, A, LDA, IPIV, INFO )
c
c     June 30, 1992
c
c     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
c     ..
c
c  Purpose
c  =======
c
c  SGETF2 computes an LU factorization of a general m-by-n matrix A
c  using partial pivoting with row interchanges.
c
c  The factorization has the form
c     A = P * L * U
c  where P is a permutation matrix, L is lower triangular with unit
c  diagonal elements (lower trapezoidal if m > n), and U is upper
c  triangular (upper trapezoidal if m < n).
c
c  This is the right-looking Level 2 BLAS version of the algorithm.
c
c  Arguments
c  =========
c
c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.
c
c  A       (input/output) REAL array, dimension (LDA,N)
c          On entry, the m by n matrix to be factored.
c          On exit, the factors L and U from the factorization
c          A = P*L*U; the unit diagonal elements of L are not stored.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,M).
c
c  IPIV    (output) INTEGER array, dimension (min(M,N))
c          The pivot indices; for 1 <= i <= min(M,N), row i of the
c          matrix was interchanged with row IPIV(i).
c
c  INFO    (output) INTEGER
c          = 0: successful exit
c          < 0: if INFO = -k, the k-th argument had an illegal value
c          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
c               has been completed, but the factor U is exactly
c               singular, and division by zero will occur if it is used
c               to solve a system of equations.
c
c  =====================================================================
c
c     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            J, JP
c     ..
c     .. External Functions ..
      INTEGER            ISAMAX
      EXTERNAL           ISAMAX
c     ..
c     .. External Subroutines ..
      EXTERNAL           SGER, SSCAL, SSWAP, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..


c     Test the input parameters.

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGETF2', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 )  RETURN

      DO 10 J = 1, MIN( M, N )

c        Find pivot and test for singularity.

         JP = J - 1 + ISAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN

c           Apply the interchange to columns 1:N.

            IF( JP.NE.J )
     &         CALL SSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )

c           Compute elements J+1:M of J-th column.

            IF( J.LT.M )
     &         CALL SSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )

         ELSE IF( INFO.EQ.0 ) THEN

            INFO = J
         END IF

         IF( J.LT.MIN( M, N ) ) THEN

c           Update trailing submatrix.

            CALL SGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     &                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN

c           End of SGETF2
      END


      SUBROUTINE SLASWP( N, A, LDA, K1, K2, IPIV, INCX )
c
c     October 31, 1992
c
c     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
c     ..
c
c  Purpose
c  =======
c
c  SLASWP performs a series of row interchanges on the matrix A.
c  One row interchange is initiated for each of rows K1 through K2 of A.
c
c  Arguments
c  =========
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.
c
c  A       (input/output) REAL array, dimension (LDA,N)
c          On entry, the matrix of column dimension N to which the row
c          interchanges will be applied.
c          On exit, the permuted matrix.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.
c
c  K1      (input) INTEGER
c          The first element of IPIV for which a row interchange will
c          be done.
c
c  K2      (input) INTEGER
c          The last element of IPIV for which a row interchange will
c          be done.
c
c  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
c          The vector of pivot indices.  Only the elements in positions
c          K1 through K2 of IPIV are accessed.
c          IPIV(K) = L implies rows K and L are to be interchanged.
c
c  INCX    (input) INTEGER
c          The increment between successive values of IPIV.  If IPIV
c          is negative, the pivots are applied in reverse order.
c
c =====================================================================
c
c     .. Local Scalars ..
      INTEGER            I, IP, IX
c     ..
c     .. External Subroutines ..
      EXTERNAL           SSWAP
c     ..


c     Interchange row I with row IPIV(I) for each of rows K1 through K2.

      IF( INCX.EQ.0 )
     &   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF

      IF( INCX.EQ.1 ) THEN

         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I )
     &         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE

      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I )
     &         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE

      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I )
     &         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF

      RETURN

c           End of SLASWP
      END


      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
c
C       ** SPECIALIZED TO JUST THOSE CASES NEEDED BY DISORT V1.3
C       ** (ISPEC=1 and NAME= SGBTRF or SGETRF)
c       ** Only possible values are 1, 32, and 64 in this case.
c
c     September 30, 1994
c
c     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
c     ..
c
c  Purpose
c  =======
c
c  ILAENV is called from the LAPACK routines to choose problem-dependent
c  parameters for the local environment.  See ISPEC for a description of
c  the parameters.
c
c  This version provides a set of parameters which should give good,
c  but not optimal, performance on many of the currently available
c  computers.  Users are encouraged to modify this subroutine to set
c  the tuning parameters for their particular machine using the option
c  and problem size information in the arguments.
c
c  This routine will not function correctly if it is converted to all
c  lower case.  Converting it to all upper case is allowed.
c
c  Arguments
c  =========
c
c  ISPEC   (input) INTEGER
c          Specifies the parameter to be returned as the value of
c          ILAENV.
c          = 1: the optimal blocksize; if this value is 1, an unblocked
c               algorithm will give the best performance.
c
c  NAME    (input) CHARACTER*(*)
c          The name of the calling subroutine, in either upper case or
c          lower case.
c
c  OPTS    (input) CHARACTER*(*)
c          The character options to the subroutine NAME, concatenated
c          into a single character string.  For example, UPLO = 'U',
c          TRANS = 'T', and DIAG = 'N' for a triangular routine would
c          be specified as OPTS = 'UTN'.
c
c  N1      (input) INTEGER
c  N2      (input) INTEGER
c  N3      (input) INTEGER
c  N4      (input) INTEGER
c          Problem dimensions for the subroutine NAME; these may not all
c          be required.
c
c (ILAENV) (output) INTEGER
c          >= 0: the value of the parameter specified by ISPEC
c          < 0:  ISPEC had an illegal value.
c
c  Further Details
c  ===============
c
c  The following conventions have been used when calling ILAENV from the
c  LAPACK routines:
c  1)  OPTS is a concatenation of all of the character options to
c      subroutine NAME, in the same order that they appear in the
c      argument list for NAME, even if they are not used in determining
c      the value of the parameter specified by ISPEC.
c  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
c      that they appear in the argument list for NAME.  N1 is used
c      first, N2 second, and so on, and unused problem dimensions are
c      passed a value of -1.
c  3)  The parameter value returned by ILAENV is checked for validity in
c      the calling subroutine.  For example, ILAENV is used to retrieve
c      the optimal blocksize for STRTRI as follows:
c
c      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
c      IF( NB.LE.1 ) NB = MAX( 1, N )
c
c  =====================================================================
c
c     .. Local Scalars ..
      CHARACTER*2        C2
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR
c     ..

      IF( ISPEC.NE.1 ) THEN
c                              Invalid value for ISPEC
         ILAENV = -1
         RETURN
      END IF

c     Convert NAME to upper case if the first character is lower case.

      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN

c                          ASCII character set

         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     &            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF

      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN

c                          EBCDIC character set

         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     &       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     &       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     &             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     &             ( IC.GE.162 .AND. IC.LE.169 ) )
     &            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF

      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN

c        Prime machines:  ASCII+128

         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     &            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF

      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )

      NB = 1

      IF( C2.EQ.'GE' ) THEN
      
         IF( C3.EQ.'TRF' ) THEN
               NB = 64
         END IF
         
      ELSE IF( C2.EQ.'GB' ) THEN
      
         IF( C3.EQ.'TRF' ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
         END IF
         
      END IF
      
      ILAENV = NB

      RETURN

c          End of ILAENV
      END

