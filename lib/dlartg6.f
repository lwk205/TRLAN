      SUBROUTINE DLARTG6( K, L, M, N, CS, SN )
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   K, L, M, N
      DOUBLE PRECISION   CS, R, SN
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
*     LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
      DOUBLE PRECISION   K1, L1, M1, N1
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
*     ..
*     .. Save statement ..
*     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
*     ..
*     .. Data statements ..
*     DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
*     IF( FIRST ) THEN
      SAFMIN = DLAMCH( 'S' )
      EPS = DLAMCH( 'E' )
      SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $     LOG( DLAMCH( 'B' ) ) / TWO )
      SAFMX2 = ONE / SAFMN2
*     FIRST = .FALSE.
*     END IF
      
      K1 = K
      M1 = M
      SCALE = MAX( K1, M1 )
      IF( SCALE.GE.SAFMX2 ) THEN
         COUNT = 0
 50      CONTINUE
         COUNT = COUNT + 1
         K1 = K1*SAFMN2
         M1 = M1*SAFMN2
         SCALE = MAX( K1, M1 )
         IF( SCALE.GE.SAFMX2 )
     $        GO TO 50
         L1 = L
         N1 = N
         DO 60 I = 1, COUNT
            L1 = L1*SAFMN2
            N1 = N1*SAFMN2
 60      CONTINUE
      ELSE IF( SCALE.LE.SAFMN2 ) THEN
         COUNT = 0
 70      CONTINUE
         COUNT = COUNT + 1
         K1 = K1*SAFMX2
         M1 = M1*SAFMX2
         SCALE = MAX( K1, M1 )
         IF( SCALE.LE.SAFMN2 )
     $        GO TO 70
         L1 = L
         N1 = N
         DO 80 I = 1, COUNT
            L1 = L1*SAFMX2
            N1 = N1*SAFMX2
 80      CONTINUE
      ELSE
         L1 = L
         N1 = N
      END IF
      
      F1=K1*L1+M1*N1
      G1=M1*L1-K1*N1
      IF( G1 .GE. ZERO ) THEN
         CALL DLARTG(F1,G1,CS,SN,R)
      ELSE
         CALL DLARTG(F1,-G1,CS,SN,R)
         SN = -SN
      END IF
      RETURN
*     
*     End of DLARTG6
*     
      END
