SUBROUTINE RESGKL(J,MODE,IAP,JA,A,N,M,K,TM,Q,P_PLUS,INFO,SELEK,WORK,LWORK)
  IMPLICIT NONE

  INTEGER IAP(*), JA(*), LWORK
  DOUBLE PRECISION A(*), WORK(*)
  CHARACTER MODE

  INTEGER N,M,K,I,J,INFO,SELEK
  DOUBLE PRECISION ONE,ZERO,TWO,MINUSONE
  PARAMETER(ONE=1.0D+0,ZERO=0.0D+0,TWO=2.0D+0,MINUSONE=-1.0D+0)
  DOUBLE PRECISION CDUMMY(1,1),ALPHA,BETA
  DOUBLE PRECISION TM(M,M),TMP_M_M(M,M),Q(N,M),P_PLUS(N),BD(M),BE(M)
  DOUBLE PRECISION TMP_N_K(N,K),WORK2(M*8)
  DOUBLE PRECISION Y(M,M),P(N)
  DOUBLE PRECISION DDOT,DNRM2,DLAMCH

  Q(1:N,J) = P_PLUS(1:N)
  DO WHILE(J < M+1)
!     write(*,*)TM
     CALL av(N,IAP,JA,A, Q(1:N,J), P)
     CALL CGS2(P,Q,N,J-2,WORK)  !J=-1のとき,MKLがエラーを吐く
     ALPHA = DDOT(N,P,1,Q(1:N,J),1)
     TM(J,J) = ALPHA
     IF (J .NE. 1) THEN
        CALL DAXPY(N,-TM(J-1,J),Q(1:N,J-1),1,P,1) 
     END IF
     CALL DAXPY(N,-TM(J,J),Q(1:N,J),1,P,1)
     BETA=DNRM2(N,P,1)
     IF(J < M) THEN
        TM(J,J+1)=BETA
        TM(J+1,J)=BETA
        Q(1:N,J+1)=P/BETA
     ELSE
        P_PLUS(1:N) = P / BETA
     END IF
     J=J+1
  END DO
  J = K+1

  Y = ZERO
  DO I = 1,M
     Y(I,I) = ONE
  END DO

  IF (SELEK == 2) THEN
  write(*,*) "DSYEV ver"
     CALL DSYEV("V", "U", M, TM, M, BD, WORK, LWORK, INFO )
     Y=TM
     CALL DGEMM('N','N',N,K,M,ONE,Q,N,Y,M,ZERO,TMP_N_K,N)
     CALL DCOPY(N*K,TMP_N_K,1,Q,1)

     TM = ZERO
     DO I =1 ,K
        TM(I,I) = BD(I)
     END DO
     CALL DAXPY(K,BETA,Y(M,1),M,TM(1,K+1),1)
  END IF

  IF (SELEK == 1) THEN
  write(*,*) "DGEBRDG_4_TRISIDE + DSTEQR ver"
     CALL DGEBRDG_4_TRISIDE(K+1,TM,M,Y)
     !tridiagonalizeされる前のTM == Y * された後のTM * Y^T
     DO I =1 ,M
        BD(I)=TM(I,I)
     END DO
     DO I =1 ,M-1
        BE(I)=TM(I,I+1)
     END DO
     CALL DSTEQR("V",M,BD,BE,Y,M,WORK,INFO) !これでよい？
     CALL DGEMM('N','N',N,K,M,ONE,Q,N,Y,M,ZERO,TMP_N_K,N)
     CALL DCOPY(N*K,TMP_N_K,1,Q,1)

     TM = ZERO
     DO I =1 ,K
        TM(I,I) = BD(I)
     END DO
     CALL DAXPY(K,BETA,Y(M,1),M,TM(1,K+1),1)
     CALL DAXPY(K,BETA,Y(M,1),M,TM(K+1,1),M)
  END IF


  IF (SELEK == 3) THEN
     write(*,*) "DGEBRDG_4_TRISIDE + DSTEQR (転地合成考慮版) ver"
     CALL DGEBRDG_4_TRISIDE(K+1,TM,M,Y)
     !tridiagonalizeされる前のTM == Y * された後のTM * Y^T
     DO I =1 ,M
        BD(I)=TM(I,I)
     END DO
     DO I =1 ,M-1
        BE(I)=(TM(I,I+1)+TM(I+1,I))/2.0
        !write(*,*) TM(I,I+1)-TM(I+1,I)
     END DO
     !write(*,*) TM
     CALL DSTEQR("V",M,BD,BE,Y,M,WORK,INFO) !これでよい？
     CALL DGEMM('N','N',N,K,M,ONE,Q,N,Y,M,ZERO,TMP_N_K,N)
     CALL DCOPY(N*K,TMP_N_K,1,Q,1)

     TM = ZERO
     DO I =1 ,K
        TM(I,I) = BD(I)
     END DO
     CALL DAXPY(K,BETA,Y(M,1),M,TM(1,K+1),1)
     CALL DAXPY(K,BETA,Y(M,1),M,TM(K+1,1),M)
  END IF
  Q(1:N,J) = P_PLUS(1:N)
  CALL av(N,IAP,JA,A, Q(1:N,J), P)
  ALPHA = DDOT(N,P,1,Q(1:N,J),1)
  TM(J,J) = ALPHA
  CALL DGEMV('N',N,K,-BETA,Q(1,1),N,Y(M,1),M,ONE,P,1)
  CALL DAXPY(N,-TM(J,J),Q(1:N,J),1,P,1)
  BETA=DNRM2(N,P,1)
  TM(J,J+1)=BETA
  TM(J+1,J)=BETA
  P_PLUS(1:N)=P/BETA
  J=J+1

  RETURN
END SUBROUTINE RESGKL

subroutine av (M, IAP, JA ,A, P, AP)
      
  IMPLICIT NONE
  include 'omp_lib.h'
  integer M
  integer IAP(*),JA(*)
  double precision A(*),P(*),AP(*),zero
  parameter (zero=0.0d0)
  integer i,j
  
  !$OMP PARALLEL DO PRIVATE(j)
  do i=1,M
     AP(i)=zero
     do j=IAP(i),IAP(i+1)-1
        AP(i)=AP(i)+(A(j))*P(JA(j))
     enddo
  enddo
  !$OMP END PARALLEL DO
  
  return
end subroutine av

subroutine atv (M, N, IAP, JA, A, Q, AQ)
  
  IMPLICIT NONE
  include 'omp_lib.h'
  integer M,N
  integer IAP(*),JA(*)
  double precision A(*),Q(*),AQ(N),zero
  parameter (zero=0.0d0)
  
  integer i,j
  
  do i=1,N
     AQ(i)=zero
  enddo
  
  !$OMP PARALLEL DO PRIVATE(j) REDUCTION(+:AQ)
  do i=1,M
     do j=IAP(i),IAP(i+1)-1
        AQ(JA(j))=AQ(JA(j))+(A(j))*Q(i)
     enddo
  enddo
  !$OMP END PARALLEL DO
  
  return
end subroutine atv

