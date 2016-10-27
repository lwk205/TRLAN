SUBROUTINE GKL(mode,matdescra,indxA,pntrbA,pntreA,A,M,N,K,B,V,U,VNTEMP,beta)
  IMPLICIT NONE

  INTEGER indxA(*), pntrbA(*), pntreA(*)
  DOUBLE PRECISION A(*)
  CHARACTER*6 matdescra
  CHARACTER MODE

  INTEGER M,N,J,K
  INTEGER ISEED(4),ONEE
  PARAMETER (ONEE = 1)
  DOUBLE PRECISION  ONE, ZERO, beta
  PARAMETER ( ONE = 1.0D+0, ZERO=0.0D+0 )
  DOUBLE PRECISION NRM
  DOUBLE PRECISION B(K,K),V(N,K),U(M,K),VN(N),VNTEMP(N),VMTEMP(M)
  DOUBLE PRECISION DNRM2
  EXTERNAL DNRM2,DGEMV,P3,CGS2,DGEMM
  B=0
  V=0
  U=0

  ISEED( 1 ) = 1
  ISEED( 2 ) = 3
  ISEED( 3 ) = 5
  ISEED( 4 ) = 7

  CALL DLARNV(1,ISEED,N,VN)
  NRM = DNRM2(N,VN,1)
  VN = VN/NRM

  V(:,1)=VN

  J=0
  DO WHILE ( J < K )
     VNTEMP = V(:,J+1)
     IF(MODE=='s') THEN
        CALL MKL_DCSRMV( 'N', M, N, ONE, matdescra, A, indxA, &
             pntrbA,pntreA, VNTEMP(1), ZERO, VMTEMP(1))
     ELSE IF(MODE=='d') THEN
        CALL DGEMV('N',M,N,ONE,A,M,VNTEMP,1,ZERO,VMTEMP,1)
     END IF
     CALL CGS2(VMTEMP,U,M,K,J)
     NRM = DNRM2(M,VMTEMP,1)
     VMTEMP = VMTEMP / NRM
     B(J+1,J+1)=NRM
     U(:,J+1)=VMTEMP

     IF(MODE=='s') THEN
        CALL MKL_DCSRMV( 'T', M, N, ONE, matdescra, A, indxA, &
             pntrbA,pntreA, VMTEMP(1), ZERO, VNTEMP(1))
     ELSE IF(MODE=='d') THEN
        CALL DGEMV('T',M,N,ONE,A,M,VMTEMP,1,ZERO,VNTEMP,1)
     END IF
     CALL CGS2(VNTEMP,V,N,K,J+1)
     beta = DNRM2(N,VNTEMP,1)
     VNTEMP = VNTEMP / beta 
     IF ( J < K - 1 ) THEN
        B(J+1,J+2)=beta
        V(:,J+2)=VNTEMP
     END IF
     J=J+1
  END DO

  RETURN
END SUBROUTINE GKL
