      SUBROUTINE DBDSQRU( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,
     $                   LDU, C, LDC, WORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          UPLO
      INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ),
     $                   VT( LDVT, * ), WORK( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DBDSQRU COMPUTES THE SINGULAR VALUE DECOMPOSITION (SVD) OF A REAL
*  N-BY-N (UPPER OR LOWER) BIDIAGONAL MATRIX WITH DIAGONAL D AND
*  OFFDIAGONAL E, ACCUMULATING THE TRANSFORMATIONS IF DESIRED. LETTING
*  B DENOTE THE INPUT BIDIAGONAL MATRIX, THE ALGORITHM COMPUTES
*  ORTHOGONAL MATRICES Q AND P SUCH THAT B = Q * S * P' (P' DENOTES THE
*  TRANSPOSE OF P). THE SINGULAR VALUES S ARE OVERWRITTEN ON D.
*
*  THE INPUT MATRIX U  IS CHANGED TO U  * Q  IF DESIRED.
*  THE INPUT MATRIX VT IS CHANGED TO P' * VT IF DESIRED.
*  THE INPUT MATRIX C  IS CHANGED TO Q' * C  IF DESIRED.
*
*  SEE "COMPUTING  SMALL SINGULAR VALUES OF BIDIAGONAL MATRICES WITH
*  GUARANTEED HIGH RELATIVE ACCURACY," BY J. DEMMEL AND W. KAHAN,
*  LAPACK WORKING NOTE #3, FOR A DETAILED DESCRIPTION OF THE ALGORITHM.
*
*  ARGUMENTS
*  =========
*
*  UPLO    (INPUT) CHARACTER*1
*          ON ENTRY, UPLO SPECIFIES WHETHER THE INPUT BIDIAGONAL MATRIX
*          IS UPPER OR LOWER BIDIAGONAL.
*             UPLO = 'U' OR 'U'   B IS UPPER BIDIAGONAL.
*             UPLO = 'L' OR 'L'   B IS LOWER BIDIAGONAL.
*
*  N       (INPUT) INTEGER
*          ON ENTRY, N SPECIFIES THE NUMBER OF ROWS AND COLUMNS
*          IN THE MATRIX. N MUST BE AT LEAST 0.
*
*  NCVT    (INPUT) INTEGER
*          ON ENTRY, NCVT SPECIFIES THE NUMBER OF COLUMNS OF
*          THE MATRIX VT. NCVT MUST BE AT LEAST 0.
*
*  NRU     (INPUT) INTEGER
*          ON ENTRY, NRU SPECIFIES THE NUMBER OF ROWS OF
*          THE MATRIX U. NRU MUST BE AT LEAST 0.
*
*  NCC     (INPUT) INTEGER
*          ON ENTRY, NCC SPECIFIES THE NUMBER OF COLUMNS OF
*          THE MATRIX C. NCC MUST BE AT LEAST 0.
*
*  D       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          ON ENTRY, D CONTAINS THE DIAGONAL ENTRIES OF THE
*          BIDIAGONAL MATRIX WHOSE SVD IS DESIRED. ON NORMAL EXIT,
*          D CONTAINS THE SINGULAR VALUES IN DECREASING ORDER.
*
*  E       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          ON ENTRY, THE ENTRIES OF E CONTAIN THE
*          OFFDIAGONAL ENTRIES OF OF THE BIDIAGONAL MATRIX
*          WHOSE SVD IS DESIRED. ON NORMAL EXIT, E WILL CONTAIN 0.
*          IF THE ALGORITHM DOES NOT CONVERGE, D AND E WILL CONTAIN
*          THE DIAGONAL AND SUPERDIAGONAL ENTRIES OF A BIDIAGONAL
*          MATRIX ORTHOGONALLY EQUIVALENT TO THE ONE GIVEN AS INPUT.
*
*  VT      (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDVT, NCVT)
*          ON ENTRY, CONTAINS AN N-BY-NCVT MATRIX WHICH ON EXIT
*          HAS BEEN PREMULTIPLIED BY P' (NOT REFERENCED IF NCVT=0).
*
*  LDVT    (INPUT) INTEGER
*          ON ENTRY, LDVT SPECIFIES THE LEADING DIMENSION OF VT AS
*          DECLARED IN THE CALLING (SUB) PROGRAM. LDVT MUST BE AT
*          LEAST 1. IF NCVT IS NONZERO LDVT MUST ALSO BE AT LEAST N.
*
*  U       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDU, N)
*          ON ENTRY, CONTAINS AN NRU-BY-N MATRIX WHICH ON EXIT
*          HAS BEEN POSTMULTIPLIED BY Q (NOT REFERENCED IF NRU=0).
*
*  LDU     (INPUT) INTEGER
*          ON ENTRY, LDU  SPECIFIES THE LEADING DIMENSION OF U AS
*          DECLARED IN THE CALLING (SUB) PROGRAM. LDU MUST BE AT
*          LEAST MAX( 1, NRU ) .
*
*  C       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDC, NCC)
*          ON ENTRY, CONTAINS AN N-BY-NCC MATRIX WHICH ON EXIT
*          HAS BEEN PREMULTIPLIED BY Q' (NOT REFERENCED IF NCC=0).
*
*  LDC     (INPUT) INTEGER
*          ON ENTRY, LDC  SPECIFIES THE LEADING DIMENSION OF C AS
*          DECLARED IN THE CALLING (SUB) PROGRAM. LDC MUST BE AT
*          LEAST 1. IF NCC IS NONZERO, LDC MUST ALSO BE AT LEAST N.
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION
*                      (MAX( 1, 4*N-4 ))
*          WORKSPACE. ONLY REFERENCED IF ONE OF NCVT, NRU, OR NCC IS
*          NONZERO, AND IF N IS AT LEAST 2.
*
*  INFO    (OUTPUT) INTEGER
*          ON EXIT, A VALUE OF 0 INDICATES A SUCCESSFUL EXIT.
*          IF INFO < 0, ARGUMENT NUMBER -INFO IS ILLEGAL.
*          IF INFO > 0, THE ALGORITHM DID NOT CONVERGE, AND INFO
*          SPECIFIES HOW MANY SUPERDIAGONALS DID NOT CONVERGE.
*
*  INTERNAL PARAMETERS
*  ===================
*
*  TOLMUL  DOUBLE PRECISION, DEFAULT = MAX(10,MIN(100,EPS**(-1/8)))
*          TOLMUL CONTROLS THE CONVERGENCE CRITERION OF THE QR LOOP.
*          IF IT IS POSITIVE, TOLMUL*EPS IS THE DESIRED RELATIVE
*             PRECISION IN THE COMPUTED SINGULAR VALUES.
*          IF IT IS NEGATIVE, ABS(TOLMUL*EPS*SIGMA_MAX) IS THE
*             DESIRED ABSOLUTE ACCURACY IN THE COMPUTED SINGULAR
*             VALUES (CORRESPONDS TO RELATIVE ACCURACY
*             ABS(TOLMUL*EPS) IN THE LARGEST SINGULAR VALUE.
*          ABS(TOLMUL) SHOULD BE BETWEEN 1 AND 1/EPS, AND PREFERABLY
*             BETWEEN 10 (FOR FAST CONVERGENCE) AND .1/EPS
*             (FOR THERE TO BE SOME ACCURACY IN THE RESULTS).
*          DEFAULT IS TO LOSE AT EITHER ONE EIGHTH OR 2 OF THE
*             AVAILABLE DECIMAL DIGITS IN EACH COMPUTED SINGULAR VALUE
*             (WHICHEVER IS SMALLER).
*
*  MAXITR  INTEGER, DEFAULT = 6
*          MAXITR CONTROLS THE MAXIMUM NUMBER OF PASSES OF THE
*          ALGORITHM THROUGH ITS INNER LOOP. THE ALGORITHMS STOPS
*          (AND SO FAILS TO CONVERGE) IF THE NUMBER OF PASSES
*          THROUGH THE INNER LOOP EXCEEDS MAXITR*N**2.
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   NEGONE
      PARAMETER          ( NEGONE = -1.0D0 )
      DOUBLE PRECISION   HNDRTH
      PARAMETER          ( HNDRTH = 0.01D0 )
      DOUBLE PRECISION   TEN
      PARAMETER          ( TEN = 10.0D0 )
      DOUBLE PRECISION   HNDRD
      PARAMETER          ( HNDRD = 100.0D0 )
      DOUBLE PRECISION   MEIGTH
      PARAMETER          ( MEIGTH = -0.125D0 )
      INTEGER            MAXITR
      PARAMETER          ( MAXITR = 6 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            ROTATE
      INTEGER            I, IDIR, IROT, ISUB, ITER, IUPLO, J, JOB, LL,
     $                   LLL, M, MAXIT, NM1, NM12, NM13, OLDLL, OLDM
      DOUBLE PRECISION   ABSE, ABSS, COSL, COSR, CS, EPS, F, G, GAP,
     $                   GMAX, H, MU, OLDCS, OLDSN, R, SHIFT, SIGMN,
     $                   SIGMX, SINL, SINR, SLL, SMAX, SMIN, SMINL,
     $                   SMINLO, SMINOA, SN, THRESH, TOL, TOLMUL, UNFL
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLARTG, DLAS2, DLASR, DLASV2, DROT, DSCAL,
     $                   DSWAP, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SIGN, SQRT
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IUPLO = 0
      IF( LSAME( UPLO, 'U' ) )
     $   IUPLO = 1
      IF( LSAME( UPLO, 'L' ) )
     $   IUPLO = 2
      IF( IUPLO.EQ.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NCVT.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NCC.LT.0 ) THEN
         INFO = -5
      ELSE IF( ( NCVT.EQ.0 .AND. LDVT.LT.1 ) .OR.
     $         ( NCVT.GT.0 .AND. LDVT.LT.MAX( 1, N ) ) ) THEN
         INFO = -9
      ELSE IF( LDU.LT.MAX( 1, NRU ) ) THEN
         INFO = -11
      ELSE IF( ( NCC.EQ.0 .AND. LDC.LT.1 ) .OR.
     $         ( NCC.GT.0 .AND. LDC.LT.MAX( 1, N ) ) ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DBDSQRU', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 )
     $   RETURN
      IF( N.EQ.1 )
     $   GO TO 190
*
*     ROTATE IS TRUE IF ANY SINGULAR VECTORS DESIRED, FALSE OTHERWISE
*
      ROTATE = ( NCVT.GT.0 ) .OR. ( NRU.GT.0 ) .OR. ( NCC.GT.0 )
      NM1 = N - 1
      NM12 = NM1 + NM1
      NM13 = NM12 + NM1
*
*     GET MACHINE CONSTANTS
*
      EPS = DLAMCH( 'EPSILON' )
      UNFL = DLAMCH( 'SAFE MINIMUM' )
      TOLMUL = MAX( TEN, MIN( HNDRD, EPS**MEIGTH ) )
      TOL = TOLMUL*EPS
*
*     IF MATRIX LOWER BIDIAGONAL, ROTATE TO BE UPPER BIDIAGONAL
*     BY APPLYING GIVENS ROTATIONS ON THE LEFT
*
      IF( IUPLO.EQ.2 ) THEN
         DO 10 I = 1, N - 1
            CALL DLARTG( D( I ), E( I ), CS, SN, R )
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            IF( ROTATE ) THEN
               WORK( I ) = CS
               WORK( NM1+I ) = SN
            END IF
   10    CONTINUE
*
*        UPDATE SINGULAR VECTORS IF DESIRED
*
         IF( NRU.GT.0 )
     $      CALL DLASR( 'R', 'V', 'F', NRU, N, WORK( 1 ), WORK( N ), U,
     $                  LDU )
         IF( NCC.GT.0 )
     $      CALL DLASR( 'L', 'V', 'F', N, NCC, WORK( 1 ), WORK( N ), C,
     $                  LDC )
      END IF
*
*     COMPUTE APPROXIMATE MAXIMUM, MINIMUM SINGULAR VALUES
*
      SMAX = ABS( D( N ) )
      DO 20 I = 1, N - 1
         SMAX = MAX( SMAX, ABS( D( I ) ), ABS( E( I ) ) )
   20 CONTINUE
      SMINL = ZERO
      IF( TOL.GE.ZERO ) THEN
         SMINOA = ABS( D( 1 ) )
         IF( SMINOA.EQ.ZERO )
     $      GO TO 40
         MU = SMINOA
         DO 30 I = 2, N
            MU = ABS( D( I ) )*( MU / ( MU+ABS( E( I-1 ) ) ) )
            SMINOA = MIN( SMINOA, MU )
            IF( SMINOA.EQ.ZERO )
     $         GO TO 40
   30    CONTINUE
   40    CONTINUE
         SMINOA = SMINOA / SQRT( DBLE( N ) )
      END IF
*
*     PREPARE FOR MAIN ITERATION LOOP FOR THE SINGULAR VALUES
*
      MAXIT = MAXITR*N*N
      ITER = 0
      OLDLL = -1
      OLDM = -1
      IF( NCC.EQ.0 .AND. NRU.EQ.0 .AND. NCVT.EQ.0 ) THEN
*
*        NO SINGULAR VECTORS DESIRED
*
         JOB = 0
      ELSE
*
*        SINGULAR VECTORS DESIRED
*
         JOB = 1
      END IF
      IF( TOL.GE.ZERO ) THEN
*
*        RELATIVE ACCURACY DESIRED
*
         THRESH = MAX( TOL*SMINOA, MAXIT*UNFL )
      ELSE
*
*        ABSOLUTE ACCURACY DESIRED
*
         THRESH = MAX( ABS( TOL )*SMAX, MAXIT*UNFL )
      END IF
*
*     M POINTS TO LAST ENTRY OF UNCONVERGED PART OF MATRIX
*
      M = N
*
*     BEGIN MAIN ITERATION LOOP
*
   50 CONTINUE
*
*     CHECK FOR CONVERGENCE OR EXCEEDING ITERATION COUNT
*
      IF( M.LE.1 )
     $   GO TO 190
      IF( ITER.GT.MAXIT )
     $   GO TO 230
*
*     FIND DIAGONAL BLOCK OF MATRIX TO WORK ON
*
      IF( TOL.LT.ZERO .AND. ABS( D( M ) ).LE.THRESH )
     $   D( M ) = ZERO
      SMAX = ABS( D( M ) )
      SMIN = SMAX
      DO 60 LLL = 1, M
         LL = M - LLL
         IF( LL.EQ.0 )
     $      GO TO 80
         ABSS = ABS( D( LL ) )
         ABSE = ABS( E( LL ) )
         IF( TOL.LT.ZERO .AND. ABSS.LE.THRESH )
     $      D( LL ) = ZERO
         IF( ABSE.LE.THRESH )
     $      GO TO 70
         SMIN = MIN( SMIN, ABSS )
         SMAX = MAX( SMAX, ABSS, ABSE )
   60 CONTINUE
   70 CONTINUE
      E( LL ) = ZERO
*
*     MATRIX SPLITS SINCE E(LL) = 0
*
      IF( LL.EQ.M-1 ) THEN
*
*        CONVERGENCE OF BOTTOM SINGULAR VALUE, RETURN TO TOP OF LOOP
*
         M = M - 1
         GO TO 50
      END IF
   80 CONTINUE
      LL = LL + 1
*
*     E(LL) THROUGH E(M-1) ARE NONZERO, E(LL-1) IS ZERO
*
      IF( LL.EQ.M-1 ) THEN
*
*        2 BY 2 BLOCK, HANDLE SEPARATELY
*
         CALL DLASV2( D( M-1 ), E( M-1 ), D( M ), SIGMN, SIGMX, SINR,
     $                COSR, SINL, COSL )
         D( M-1 ) = SIGMX
         E( M-1 ) = ZERO
         D( M ) = SIGMN
*
*        COMPUTE SINGULAR VECTORS, IF DESIRED
*
         IF( NCVT.GT.0 )
     $      CALL DROT( NCVT, VT( M-1, 1 ), LDVT, VT( M, 1 ), LDVT, COSR,
     $                 SINR )
         IF( NRU.GT.0 )
     $      CALL DROT( NRU, U( 1, M-1 ), 1, U( 1, M ), 1, COSL, SINL )
         IF( NCC.GT.0 )
     $      CALL DROT( NCC, C( M-1, 1 ), LDC, C( M, 1 ), LDC, COSL,
     $                 SINL )
         M = M - 2
         GO TO 50
      END IF
*
*     IF WORKING ON NEW SUBMATRIX, CHOOSE SHIFT DIRECTION
*     (FROM LARGER END DIAGONAL ENTRY TOWARDS SMALLER)
*
      IF( LL.GT.OLDM .OR. M.LT.OLDLL ) THEN
         IF( ABS( D( LL ) ).GE.ABS( D( M ) ) ) THEN
*
*           CHASE BULGE FROM TOP (BIG END) TO BOTTOM (SMALL END)
*
            IDIR = 1
         ELSE
*
*           CHASE BULGE FROM BOTTOM (BIG END) TO TOP (SMALL END)
*
            IDIR = 2
         END IF
      END IF
*
*     APPLY CONVERGENCE TESTS
*
      IF( IDIR.EQ.1 ) THEN
*
*        RUN CONVERGENCE TEST IN FORWARD DIRECTION
*        FIRST APPLY STANDARD TEST TO BOTTOM OF MATRIX
*
         IF( ABS( E( M-1 ) ).LE.ABS( TOL )*ABS( D( M ) ) .OR.
     $       ( TOL.LT.ZERO .AND. ABS( E( M-1 ) ).LE.THRESH ) ) THEN
            E( M-1 ) = ZERO
            GO TO 50
         END IF
*
         IF( TOL.GE.ZERO ) THEN
*
*           IF RELATIVE ACCURACY DESIRED,
*           APPLY CONVERGENCE CRITERION FORWARD
*
            MU = ABS( D( LL ) )
            SMINL = MU
            DO 90 LLL = LL, M - 1
               IF( ABS( E( LLL ) ).LE.TOL*MU ) THEN
                  E( LLL ) = ZERO
                  GO TO 50
               END IF
               SMINLO = SMINL
               MU = ABS( D( LLL+1 ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMINL = MIN( SMINL, MU )
   90       CONTINUE
*
*           IF SINGULAR VALUES ONLY WANTED, APPLY GAP TEST TO BOTTOM
*           END OF MATRIX
*
            IF( JOB.EQ.0 ) THEN
               GAP = SMINLO / SQRT( DBLE( M-LL ) ) - ABS( D( M ) )
               IF( GAP.GT.ZERO ) THEN
                  ABSS = ABS( D( M ) )
                  ABSE = ABS( E( M-1 ) )
                  GMAX = MAX( GAP, ABSS, ABSE )
                  IF( ( ABSE / GMAX )**2.LE.TOL*( GAP / GMAX )*
     $                ( ABSS / GMAX ) ) THEN
                     E( M-1 ) = ZERO
                     GO TO 50
                  END IF
               END IF
            END IF
         END IF
      ELSE
*
*        RUN CONVERGENCE TEST IN BACKWARD DIRECTION
*        FIRST APPLY STANDARD TEST TO TOP OF MATRIX
*
         IF( ABS( E( LL ) ).LE.ABS( TOL )*ABS( D( LL ) ) .OR.
     $       ( TOL.LT.ZERO .AND. ABS( E( LL ) ).LE.THRESH ) ) THEN
            E( LL ) = ZERO
            GO TO 50
         END IF
*
         IF( TOL.GE.ZERO ) THEN
*
*           IF RELATIVE ACCURACY DESIRED,
*           APPLY CONVERGENCE CRITERION BACKWARD
*
            MU = ABS( D( M ) )
            SMINL = MU
            DO 100 LLL = M - 1, LL, -1
               IF( ABS( E( LLL ) ).LE.TOL*MU ) THEN
                  E( LLL ) = ZERO
                  GO TO 50
               END IF
               SMINLO = SMINL
               MU = ABS( D( LLL ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMINL = MIN( SMINL, MU )
  100       CONTINUE
*
*           IF SINGULAR VALUES ONLY WANTED, APPLY GAP TEST TO TOP
*           END OF MATRIX
*
            IF( JOB.EQ.0 ) THEN
               GAP = SMINLO / SQRT( DBLE( M-LL ) ) - ABS( D( LL ) )
               IF( GAP.GT.ZERO ) THEN
                  ABSS = ABS( D( LL ) )
                  ABSE = ABS( E( LL ) )
                  GMAX = MAX( GAP, ABSS, ABSE )
                  IF( ( ABSE / GMAX )**2.LE.TOL*( GAP / GMAX )*
     $                ( ABSS / GMAX ) ) THEN
                     E( LL ) = ZERO
                     GO TO 50
                  END IF
               END IF
            END IF
         END IF
      END IF
      OLDLL = LL
      OLDM = M
*
*     COMPUTE SHIFT.  FIRST, TEST IF SHIFTING WOULD RUIN RELATIVE
*     ACCURACY, AND IF SO SET THE SHIFT TO ZERO.
*
      IF( TOL.GE.ZERO .AND. N*TOL*( SMINL / SMAX ).LE.
     $    MAX( EPS, HNDRTH*TOL ) ) THEN
*
*        USE A ZERO SHIFT TO AVOID LOSS OF RELATIVE ACCURACY
*
         SHIFT = ZERO
      ELSE
*
*        COMPUTE THE SHIFT FROM 2-BY-2 BLOCK AT END OF MATRIX
*
         IF( IDIR.EQ.1 ) THEN
            SLL = ABS( D( LL ) )
            CALL DLAS2( D( M-1 ), E( M-1 ), D( M ), SHIFT, R )
         ELSE
            SLL = ABS( D( M ) )
            CALL DLAS2( D( LL ), E( LL ), D( LL+1 ), SHIFT, R )
         END IF
*
*        TEST IF SHIFT NEGLIGIBLE, AND IF SO SET TO ZERO
*
         IF( SLL.GT.ZERO ) THEN
            IF( ( SHIFT / SLL )**2.LT.EPS )
     $         SHIFT = ZERO
         END IF
      END IF
*
*     INCREMENT ITERATION COUNT
*
      ITER = ITER + M - LL
*
*     IF SHIFT = 0, DO SIMPLIFIED QR ITERATION
*
      IF( SHIFT.EQ.ZERO ) THEN
         IF( IDIR.EQ.1 ) THEN
*
*           CHASE BULGE FROM TOP TO BOTTOM
*
            CS = ONE
            OLDCS = ONE
*
*           SAVE COSINES AND SINES IF SINGULAR VECTORS DESIRED
*
            IF( ROTATE ) THEN
*
               CALL DLARTG( D( LL )*CS, E( LL ), CS, SN, R )
               CALL DLARTG( OLDCS*R, D( LL+1 )*SN, OLDCS, OLDSN,
     $                      D( LL ) )
               WORK( 1 ) = CS
               WORK( 1+NM1 ) = SN
               WORK( 1+NM12 ) = OLDCS
               WORK( 1+NM13 ) = OLDSN
               IROT = 1
               DO 110 I = LL + 1, M - 1
                  CALL DLARTG( D( I )*CS, E( I ), CS, SN, R )
                  E( I-1 ) = OLDSN*R
                  CALL DLARTG( OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN,
     $                         D( I ) )
                  IROT = IROT + 1
                  WORK( IROT ) = CS
                  WORK( IROT+NM1 ) = SN
                  WORK( IROT+NM12 ) = OLDCS
                  WORK( IROT+NM13 ) = OLDSN
  110          CONTINUE
               H = D( M )*CS
               D( M ) = H*OLDCS
               E( M-1 ) = H*OLDSN
*
*              UPDATE SINGULAR VECTORS
*
               IF( NCVT.GT.0 )
     $            CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ),
     $                        WORK( N ), VT( LL, 1 ), LDVT )
               IF( NRU.GT.0 )
     $            CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1,
     $                        WORK( NM12+1 ), WORK( NM13+1 ),
     $                        U( 1, LL ), LDU )
               IF( NCC.GT.0 )
     $            CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC,
     $                        WORK( NM12+1 ), WORK( NM13+1 ),
     $                        C( LL, 1 ), LDC )
*
            ELSE
*
               CALL DLARTG( D( LL )*CS, E( LL ), CS, SN, R )
               CALL DLARTG( OLDCS*R, D( LL+1 )*SN, OLDCS, OLDSN,
     $                      D( LL ) )
               DO 120 I = LL + 1, M - 1
                  CALL DLARTG( D( I )*CS, E( I ), CS, SN, R )
                  E( I-1 ) = OLDSN*R
                  CALL DLARTG( OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN,
     $                         D( I ) )
  120          CONTINUE
               H = D( M )*CS
               D( M ) = H*OLDCS
               E( M-1 ) = H*OLDSN
*
            END IF
*
*           TEST CONVERGENCE
*
            IF( ABS( E( M-1 ) ).LE.THRESH )
     $         E( M-1 ) = ZERO
*
         ELSE
*
*           CHASE BULGE FROM BOTTOM TO TOP
*
            CS = ONE
            OLDCS = ONE
*
*           SAVE COSINES AND SINES IF SINGULAR VECTORS DESIRED
*
            IF( ROTATE ) THEN
*
               CALL DLARTG( D( M )*CS, E( M-1 ), CS, SN, R )
               CALL DLARTG( OLDCS*R, D( M-1 )*SN, OLDCS, OLDSN, 
     $              D( M ) )
               WORK( M-LL ) = CS
               WORK( M-LL+NM1 ) = -SN
               WORK( M-LL+NM12 ) = OLDCS
               WORK( M-LL+NM13 ) = -OLDSN
               IROT = M - LL
               DO 130 I = M - 1, LL + 1, -1
                  CALL DLARTG( D( I )*CS, E( I-1 ), CS, SN, R )
                  E( I ) = OLDSN*R
                  CALL DLARTG( OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN,
     $                         D( I ) )
                  IROT = IROT - 1
                  WORK( IROT ) = CS
                  WORK( IROT+NM1 ) = -SN
                  WORK( IROT+NM12 ) = OLDCS
                  WORK( IROT+NM13 ) = -OLDSN
  130          CONTINUE
               H = D( LL )*CS
               D( LL ) = H*OLDCS
               E( LL ) = H*OLDSN
*
*              UPDATE SINGULAR VECTORS
*
               IF( NCVT.GT.0 )
     $            CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT,
     $                        WORK( NM12+1 ), WORK( NM13+1 ),
     $                        VT( LL, 1 ), LDVT )
               IF( NRU.GT.0 )
     $            CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ),
     $                        WORK( N ), U( 1, LL ), LDU )
               IF( NCC.GT.0 )
     $            CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ),
     $                        WORK( N ), C( LL, 1 ), LDC )
*
            ELSE
*
               CALL DLARTG( D( M )*CS, E( M-1 ), CS, SN, R )
               CALL DLARTG( OLDCS*R, D( M-1 )*SN, OLDCS, OLDSN, 
     $              D( M ) )
               DO 140 I = M - 1, LL + 1, -1
                  CALL DLARTG( D( I )*CS, E( I-1 ), CS, SN, R )
                  E( I ) = OLDSN*R
                  CALL DLARTG( OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN,
     $                         D( I ) )
  140          CONTINUE
               H = D( LL )*CS
               D( LL ) = H*OLDCS
               E( LL ) = H*OLDSN
*
            END IF
*
*           TEST CONVERGENCE
*
            IF( ABS( E( LL ) ).LE.THRESH )
     $         E( LL ) = ZERO
         END IF
      ELSE
*
*        USE NONZERO SHIFT
*
         IF( IDIR.EQ.1 ) THEN
*
*           CHASE BULGE FROM TOP TO BOTTOM
*
            F = ( ABS( D( LL ) )-SHIFT )*
     $          ( SIGN( ONE, D( LL ) )+SHIFT / D( LL ) )
            G = E( LL )
*
*           SAVE COSINES AND SINES IF SINGULAR VECTORS DESIRED
*
            IF( ROTATE ) THEN
*
               CALL DLARTG( F, G, COSR, SINR, R )
               F = COSR*D( LL ) + SINR*E( LL )
               E( LL ) = COSR*E( LL ) - SINR*D( LL )
               G = SINR*D( LL+1 )
               D( LL+1 ) = COSR*D( LL+1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( LL ) = R
               F = COSL*E( LL ) + SINL*D( LL+1 )
               D( LL+1 ) = COSL*D( LL+1 ) - SINL*E( LL )
               G = SINL*E( LL+1 )
               E( LL+1 ) = COSL*E( LL+1 )
               WORK( 1 ) = COSR
               WORK( 1+NM1 ) = SINR
               WORK( 1+NM12 ) = COSL
               WORK( 1+NM13 ) = SINL
               IROT = 1
               DO 150 I = LL + 1, M - 2
                  CALL DLARTG( F, G, COSR, SINR, R )
                  E( I-1 ) = R
                  F = COSR*D( I ) + SINR*E( I )
                  E( I ) = COSR*E( I ) - SINR*D( I )
                  G = SINR*D( I+1 )
                  D( I+1 ) = COSR*D( I+1 )
                  CALL DLARTG( F, G, COSL, SINL, R )
                  D( I ) = R
                  F = COSL*E( I ) + SINL*D( I+1 )
                  D( I+1 ) = COSL*D( I+1 ) - SINL*E( I )
                  G = SINL*E( I+1 )
                  E( I+1 ) = COSL*E( I+1 )
                  IROT = IROT + 1
                  WORK( IROT ) = COSR
                  WORK( IROT+NM1 ) = SINR
                  WORK( IROT+NM12 ) = COSL
                  WORK( IROT+NM13 ) = SINL
  150          CONTINUE
               CALL DLARTG( F, G, COSR, SINR, R )
               E( M-2 ) = R
               F = COSR*D( M-1 ) + SINR*E( M-1 )
               E( M-1 ) = COSR*E( M-1 ) - SINR*D( M-1 )
               G = SINR*D( M )
               D( M ) = COSR*D( M )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( M-1 ) = R
               F = COSL*E( M-1 ) + SINL*D( M )
               D( M ) = COSL*D( M ) - SINL*E( M-1 )
               IROT = IROT + 1
               WORK( IROT ) = COSR
               WORK( IROT+NM1 ) = SINR
               WORK( IROT+NM12 ) = COSL
               WORK( IROT+NM13 ) = SINL
               E( M-1 ) = F
*
*              UPDATE SINGULAR VECTORS
*
               IF( NCVT.GT.0 )
     $            CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ),
     $                        WORK( N ), VT( LL, 1 ), LDVT )
               IF( NRU.GT.0 )
     $            CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1,
     $                        WORK( NM12+1 ), WORK( NM13+1 ),
     $                        U( 1, LL ), LDU )
               IF( NCC.GT.0 )
     $            CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC,
     $                        WORK( NM12+1 ), WORK( NM13+1 ),
     $                        C( LL, 1 ), LDC )
*
            ELSE
*
               CALL DLARTG( F, G, COSR, SINR, R )
               F = COSR*D( LL ) + SINR*E( LL )
               E( LL ) = COSR*E( LL ) - SINR*D( LL )
               G = SINR*D( LL+1 )
               D( LL+1 ) = COSR*D( LL+1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( LL ) = R
               F = COSL*E( LL ) + SINL*D( LL+1 )
               D( LL+1 ) = COSL*D( LL+1 ) - SINL*E( LL )
               G = SINL*E( LL+1 )
               E( LL+1 ) = COSL*E( LL+1 )
               DO 160 I = LL + 1, M - 2
                  CALL DLARTG( F, G, COSR, SINR, R )
                  E( I-1 ) = R
                  F = COSR*D( I ) + SINR*E( I )
                  E( I ) = COSR*E( I ) - SINR*D( I )
                  G = SINR*D( I+1 )
                  D( I+1 ) = COSR*D( I+1 )
                  CALL DLARTG( F, G, COSL, SINL, R )
                  D( I ) = R
                  F = COSL*E( I ) + SINL*D( I+1 )
                  D( I+1 ) = COSL*D( I+1 ) - SINL*E( I )
                  G = SINL*E( I+1 )
                  E( I+1 ) = COSL*E( I+1 )
  160          CONTINUE
               CALL DLARTG( F, G, COSR, SINR, R )
               E( M-2 ) = R
               F = COSR*D( M-1 ) + SINR*E( M-1 )
               E( M-1 ) = COSR*E( M-1 ) - SINR*D( M-1 )
               G = SINR*D( M )
               D( M ) = COSR*D( M )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( M-1 ) = R
               F = COSL*E( M-1 ) + SINL*D( M )
               D( M ) = COSL*D( M ) - SINL*E( M-1 )
               E( M-1 ) = F
*
            END IF
*
*           TEST CONVERGENCE
*
            IF( ABS( E( M-1 ) ).LE.THRESH )
     $         E( M-1 ) = ZERO
*
         ELSE
*
*           CHASE BULGE FROM BOTTOM TO TOP
*
            F = ( ABS( D( M ) )-SHIFT )*( SIGN( ONE, D( M ) )+SHIFT /
     $          D( M ) )
            G = E( M-1 )
*
*           SAVE COSINES AND SINES IF SINGULAR VECTORS DESIRED
*
            IF( ROTATE ) THEN
*
               CALL DLARTG( F, G, COSR, SINR, R )
               F = COSR*D( M ) + SINR*E( M-1 )
               E( M-1 ) = COSR*E( M-1 ) - SINR*D( M )
               G = SINR*D( M-1 )
               D( M-1 ) = COSR*D( M-1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( M ) = R
               F = COSL*E( M-1 ) + SINL*D( M-1 )
               D( M-1 ) = COSL*D( M-1 ) - SINL*E( M-1 )
               G = SINL*E( M-2 )
               E( M-2 ) = COSL*E( M-2 )
               WORK( M-LL ) = COSR
               WORK( M-LL+NM1 ) = -SINR
               WORK( M-LL+NM12 ) = COSL
               WORK( M-LL+NM13 ) = -SINL
               IROT = M - LL
               DO 170 I = M - 1, LL + 2, -1
                  CALL DLARTG( F, G, COSR, SINR, R )
                  E( I ) = R
                  F = COSR*D( I ) + SINR*E( I-1 )
                  E( I-1 ) = COSR*E( I-1 ) - SINR*D( I )
                  G = SINR*D( I-1 )
                  D( I-1 ) = COSR*D( I-1 )
                  CALL DLARTG( F, G, COSL, SINL, R )
                  D( I ) = R
                  F = COSL*E( I-1 ) + SINL*D( I-1 )
                  D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 )
                  G = SINL*E( I-2 )
                  E( I-2 ) = COSL*E( I-2 )
                  IROT = IROT - 1
                  WORK( IROT ) = COSR
                  WORK( IROT+NM1 ) = -SINR
                  WORK( IROT+NM12 ) = COSL
                  WORK( IROT+NM13 ) = -SINL
  170          CONTINUE
               CALL DLARTG( F, G, COSR, SINR, R )
               E( LL+1 ) = R
               F = COSR*D( LL+1 ) + SINR*E( LL )
               E( LL ) = COSR*E( LL ) - SINR*D( LL+1 )
               G = SINR*D( LL )
               D( LL ) = COSR*D( LL )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( LL+1 ) = R
               F = COSL*E( LL ) + SINL*D( LL )
               D( LL ) = COSL*D( LL ) - SINL*E( LL )
               IROT = IROT - 1
               WORK( IROT ) = COSR
               WORK( IROT+NM1 ) = -SINR
               WORK( IROT+NM12 ) = COSL
               WORK( IROT+NM13 ) = -SINL
               E( LL ) = F
*
            ELSE
*
               CALL DLARTG( F, G, COSR, SINR, R )
               F = COSR*D( M ) + SINR*E( M-1 )
               E( M-1 ) = COSR*E( M-1 ) - SINR*D( M )
               G = SINR*D( M-1 )
               D( M-1 ) = COSR*D( M-1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( M ) = R
               F = COSL*E( M-1 ) + SINL*D( M-1 )
               D( M-1 ) = COSL*D( M-1 ) - SINL*E( M-1 )
               G = SINL*E( M-2 )
               E( M-2 ) = COSL*E( M-2 )
               DO 180 I = M - 1, LL + 2, -1
                  CALL DLARTG( F, G, COSR, SINR, R )
                  E( I ) = R
                  F = COSR*D( I ) + SINR*E( I-1 )
                  E( I-1 ) = COSR*E( I-1 ) - SINR*D( I )
                  G = SINR*D( I-1 )
                  D( I-1 ) = COSR*D( I-1 )
                  CALL DLARTG( F, G, COSL, SINL, R )
                  D( I ) = R
                  F = COSL*E( I-1 ) + SINL*D( I-1 )
                  D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 )
                  G = SINL*E( I-2 )
                  E( I-2 ) = COSL*E( I-2 )
  180          CONTINUE
               CALL DLARTG( F, G, COSR, SINR, R )
               E( LL+1 ) = R
               F = COSR*D( LL+1 ) + SINR*E( LL )
               E( LL ) = COSR*E( LL ) - SINR*D( LL+1 )
               G = SINR*D( LL )
               D( LL ) = COSR*D( LL )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( LL+1 ) = R
               F = COSL*E( LL ) + SINL*D( LL )
               D( LL ) = COSL*D( LL ) - SINL*E( LL )
               E( LL ) = F
*
            END IF
*
*           TEST CONVERGENCE
*
            IF( ABS( E( LL ) ).LE.THRESH )
     $         E( LL ) = ZERO
*
*           UPDATE SINGULAR VECTORS IF DESIRED
*
            IF( NCVT.GT.0 )
     $         CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ),
     $                     WORK( NM13+1 ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 )
     $         CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ),
     $                     WORK( N ), U( 1, LL ), LDU )
            IF( NCC.GT.0 )
     $         CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ),
     $                     WORK( N ), C( LL, 1 ), LDC )
         END IF
      END IF
*
*     QR ITERATION FINISHED, GO BACK AND CHECK CONVERGENCE
*
      GO TO 50
*
*     ALL SINGULAR VALUES CONVERGED, SO MAKE THEM POSITIVE
*
  190 CONTINUE
      DO 200 I = 1, N
         IF( D( I ).LT.ZERO ) THEN
            D( I ) = -D( I )
*
*           CHANGE SIGN OF SINGULAR VECTORS, IF DESIRED
*
            IF( NCVT.GT.0 )
     $         CALL DSCAL( NCVT, NEGONE, VT( I, 1 ), LDVT )
         END IF
  200 CONTINUE
*
*     SORT THE SINGULAR VALUES INTO DECREASING ORDER (INSERTION SORT ON
*     SINGULAR VALUES, BUT ONLY ONE TRANSPOSITION PER SINGULAR VECTOR)
*
      DO 220 I = 1, N - 1
*
*        SCAN FOR SMALLEST D(I)
*
         ISUB = 1
         SMIN = D( 1 )
         DO 210 J = 2, N + 1 - I
            IF( D( J ).LE.SMIN ) THEN
               ISUB = J
               SMIN = D( J )
            END IF
  210    CONTINUE
         IF( ISUB.NE.N+1-I ) THEN
*
*           SWAP SINGULAR VALUES AND VECTORS
*
            D( ISUB ) = D( N+1-I )
            D( N+1-I ) = SMIN
            IF( NCVT.GT.0 )
     $         CALL DSWAP( NCVT, VT( ISUB, 1 ), LDVT, VT( N+1-I, 1 ),
     $                     LDVT )
            IF( NRU.GT.0 )
     $         CALL DSWAP( NRU, U( 1, ISUB ), 1, U( 1, N+1-I ), 1 )
            IF( NCC.GT.0 )
     $         CALL DSWAP( NCC, C( ISUB, 1 ), LDC, C( N+1-I, 1 ), LDC )
         END IF
  220 CONTINUE
      GO TO 250
*
*     MAXIMUM NUMBER OF ITERATIONS EXCEEDED, FAILURE TO CONVERGE
*
  230 CONTINUE
      INFO = 0
      DO 240 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  240 CONTINUE
  250 CONTINUE
      RETURN
*
*     END OF DBDSQRU
*
      END
