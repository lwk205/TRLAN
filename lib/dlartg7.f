      SUBROUTINE DLARTG7( F, H, G, R, S )
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   F, G, H, R, S
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
      DOUBLE PRECISION   EPS, F1, G1, H1, SAFMIN, SAFMN2, SAFMX2, SCALE
      DOUBLE PRECISION   FHH, FHL, G1H, G1L, P0, E0 
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
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
*        FIRST = .FALSE.
*     END IF
      IF( G.EQ.ZERO ) THEN
         R = F
         S = H
      ELSE IF( F.EQ.ZERO ) THEN
         R = G
         S = ZERO
      ELSE
         F1 = F
         H1 = H
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            H1 = H1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 )
     $         GO TO 10

            CALL SQUARE(F1,H1,FHH,FHL)
            CALL TWO_PROD2(G1,G1H,G1L)
            CALL ADD(FHH,FHL,G1H,G1L,P0,E0)
            CALL USER_SQRT(P0,E0,R,S)


c            R = SQRT( (F1+H1)**2+G1**2 )
c            S = ZERO

            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
            DO 25 I = 1, COUNT
               S = S*SAFMX2
   25       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            H1 = H1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 )
     $         GO TO 30

            CALL SQUARE(F1,H1,FHH,FHL)
            CALL TWO_PROD2(G1,G1H,G1L)
            CALL ADD(FHH,FHL,G1H,G1L,P0,E0)
            CALL USER_SQRT(P0,E0,R,S)

c            R = SQRT( (F1+H1)**2+G1**2 )
c            S = ZERO

            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
            DO 45 I = 1, COUNT
               S = S*SAFMN2
   45       CONTINUE
         ELSE

            CALL SQUARE(F1,H1,FHH,FHL)
            CALL TWO_PROD2(G1,G1H,G1L)
            CALL ADD(FHH,FHL,G1H,G1L,P0,E0)
            CALL USER_SQRT(P0,E0,R,S)

c            R = SQRT( (F1+H1)**2+G1**2 )
c            S = ZERO

         END IF
      END IF
      RETURN
*
*     End of DLARTG7
*
      END

      SUBROUTINE FAST_TWO_SUM(A0,B0,S0,E0)
      IMPLICIT NONE
      DOUBLE PRECISION A0,B0,S0,E0
      S0=A0+B0
      E0=B0-(S0-A0)
      RETURN
      END SUBROUTINE FAST_TWO_SUM
      
      SUBROUTINE TWO_SUM(A0,B0,S0,E0)
      IMPLICIT NONE
      DOUBLE PRECISION A0,B0,S0,E0,V0
      S0=A0+B0
      V0=S0-A0
      E0=(A0-(S0-V0))+(B0-V0)
      RETURN
      END SUBROUTINE TWO_SUM
      
      SUBROUTINE ADD(AH,AL,BH,BL,CH,CL)
      IMPLICIT NONE
      DOUBLE PRECISION AH,AL,BH,BL,CH,CL,SH,EH
      CALL TWO_SUM(AH,BH,SH,EH)
      EH=EH+AL+BL
      CALL FAST_TWO_SUM(SH,EH,CH,CL)
      RETURN
      END SUBROUTINE ADD
      
      SUBROUTINE SPLIT(A0,H0,L0)
      IMPLICIT NONE
      DOUBLE PRECISION A0,H0,L0,T0
      DOUBLE PRECISION CONST_SPLIT
      PARAMETER          ( CONST_SPLIT = 134217729.0D0 )
      T0=CONST_SPLIT*A0
      H0=T0-(T0-A0)
      L0=A0-H0
      RETURN
      END SUBROUTINE SPLIT

      SUBROUTINE TWO_PROD(A0,B0,P0,E0)
      IMPLICIT NONE
      DOUBLE PRECISION A0,B0,AH,AL,BH,BL,P0,E0
      P0=A0*B0
      CALL SPLIT(A0,AH,AL)
      CALL SPLIT(B0,BH,BL)
      E0=((AH*BH-P0)+AH*BL+AL*BH)+AL*BL
      RETURN
      END SUBROUTINE TWO_PROD

      SUBROUTINE MUL(AH,AL,BH,BL,CH,CL)
      IMPLICIT NONE
      DOUBLE PRECISION AH,AL,BH,BL,CH,CL,P1,P2
      CALL TWO_PROD(AH,BH,P1,P2)
      P2=P2+(AH*BL)
      P2=P2+(AL*BH)
      CALL FAST_TWO_SUM(P1,P2,CH,CL)
      RETURN
      END SUBROUTINE MUL 

      SUBROUTINE TWO_PROD2(A0,P0,E0)
      IMPLICIT NONE
      DOUBLE PRECISION A0,AH,AL,P0,E0
      DOUBLE PRECISION TWO
      PARAMETER          ( TWO = 2.0D0 )
      P0=A0**2
      CALL SPLIT(A0,AH,AL)
      E0=((AH**2-P0)+TWO*AH*AL)+AL**2
      RETURN
      END SUBROUTINE TWO_PROD2

      SUBROUTINE SQUARE(AH,AL,CH,CL)
      IMPLICIT NONE
      DOUBLE PRECISION AH,AL,CH,CL,P1,P2
      DOUBLE PRECISION TWO
      PARAMETER          ( TWO = 2.0D0 )
      CALL TWO_PROD2(AH,P1,P2)
      P2=P2+TWO*(AH*AL)
      CALL FAST_TWO_SUM(P1,P2,CH,CL)
      RETURN
      END SUBROUTINE SQUARE

      SUBROUTINE USER_SQRT(AH,AL,BH,BL)
      IMPLICIT NONE
      DOUBLE PRECISION AH,AL,BH,BL,X,AXH,AXL,CH,CL,DH,DL,EH,EL
      DOUBLE PRECISION ONE,HALF,ZERO
      PARAMETER          ( ONE = 1.0D0, HALF=0.50D0, ZERO=0.0D0 )

      X=ONE/SQRT(AH)
      CALL MUL(AH,AL,X,ZERO,AXH,AXL)
      CALL SQUARE(AXH,AXL,CH,CL)
      CALL ADD(AH,AL,-CH,-CL,DH,DL)
      CALL MUL(DH,DL,X*HALF,ZERO,EH,EL)
      CALL ADD(AXH,AXL,EH,EL,BH,BL)

      RETURN
      END SUBROUTINE USER_SQRT
