SUBROUTINE RESL(start_row,which,J,MODE,IAP,JA,A,N,M,K,TM,Q,P_PLUS,INFO,SELEK,WORK,LWORK)
  IMPLICIT NONE

  INTEGER start_row(*)
  INTEGER IAP(*), JA(*), LWORK, which
  DOUBLE PRECISION A(*), WORK(*)
  CHARACTER MODE

  INTEGER N,M,K,I,JJ,J,INFO,SELEK,II,Z
  DOUBLE PRECISION ONE,ZERO,TWO,MINUSONE,PP,PQ
  PARAMETER(ONE=1.0D+0,ZERO=0.0D+0,TWO=2.0D+0,MINUSONE=-1.0D+0)
  DOUBLE PRECISION ALPHA,BETA
  DOUBLE PRECISION TM(M,M),TMP_M_M(M,M),Q(N,M),P_PLUS(N),BD(M),BE(M)
  DOUBLE PRECISION TMP_N_K(N,K)
  DOUBLE PRECISION Y(M,M),P(N)
  DOUBLE PRECISION DDOT,DNRM2

  Q(1:N,J) = P_PLUS(1:N)
  DO WHILE(J < M+1)
     CALL av(start_row,N,IAP,JA,A, Q(1:N,J), P)

     ALPHA = DDOT(N,P,1,Q(1:N,J),1)
     TM(J,J) = ALPHA
     CALL DAXPY(N,-ALPHA,Q(1:N,J),1,P,1)

     CALL CGS2(P,Q,N,J-1,WORK)  !J=-1のとき,MKLがエラーを吐く

     ALPHA = DDOT(N,P,1,Q(1:N,J),1)
     TM(J,J) = TM(J,J)+ALPHA
     CALL DAXPY(N,-ALPHA,Q(1:N,J),1,P,1)

     BETA=DNRM2(N,P,1)

     !IF (J .NE. 1) THEN
     !   CALL DAXPY(N,-TM(J-1,J),Q(1:N,J-1),1,P,1) 
     !END IF

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
     !  write(*,*) "DSYEV ver"
     CALL DSYEV("V", "U", M, TM, M, BD, WORK, LWORK, INFO ) !LWORK >= max(1,3*M-1).
!     print *,BD(96),BD(97),BD(98),BD(99),BD(100),BD(101),BD(102)
     if (which == 1) then
         DO II = 2, M
            I = II - 1
            Z = I
            PP= abs(BD( I ))
            DO JJ = II, M
               IF(abs( BD( JJ )).GT.PP) THEN
                  Z = JJ
                  PP= abs(BD( JJ ))
               END IF
            END DO
            IF( Z.NE.I ) THEN
               PQ=BD(Z)
               BD( Z ) = BD( I )
               BD( I ) = PQ
               CALL DSWAP( M, TM( 1, I ), 1, TM( 1, Z ), 1 )
            END IF
         END DO
     else if (which == 2) then
         DO II = 2, M
            I = II - 1
            Z = I
            PP= abs(BD( I ))
            DO JJ = II, M
               IF(abs( BD( JJ )).LT.PP) THEN
                  Z = JJ
                  PP= abs(BD( JJ ))
               END IF
            END DO
            IF( Z.NE.I ) THEN
               PQ=BD(Z)
               BD( Z ) = BD( I )
               BD( I ) = PQ
               CALL DSWAP( M, TM( 1, I ), 1, TM( 1, Z ), 1 )
            END IF
         END DO
     else if (which == 3) then
         DO II = 2, M
            I = II - 1
            Z = I
            PP= BD( I )
            DO JJ = II, M
               IF( BD( JJ ).GT.PP) THEN
                  Z = JJ
                  PP= BD( JJ )
               END IF
            END DO
            IF( Z.NE.I ) THEN
               PQ=BD(Z)
               BD( Z ) = BD( I )
               BD( I ) = PQ
               CALL DSWAP( M, TM( 1, I ), 1, TM( 1, Z ), 1 )
            END IF
         END DO
     end if 
     Y=TM
     CALL DGEMM('N','N',N,K,M,ONE,Q,N,Y,M,ZERO,TMP_N_K,N)
     CALL DCOPY(N*K,TMP_N_K,1,Q,1)
     TM = ZERO
     DO I =1 ,M
        write(*,*)  I, BD(I)
     END DO
     DO I =1 ,K
        TM(I,I) = BD(I)
     END DO
!     CALL DAXPY(K,BETA,Y(M,1),M,TM(1,K+1),1)
  END IF

  IF (SELEK == 1) THEN
     !  write(*,*) "katagawrisuta-to ver"
     CALL DCOPY(M*M,TM,1,TMP_M_M,1)
     CALL DSYEV("V", "U", M, TM, M, BE, WORK, LWORK, INFO ) !LWORK >= max(1,3*M-1).
     CALL DGEQRF(M,K,TM,M,BD,WORK,LWORK,INFO) ! lqork >= K
     BE(1:M) = ZERO
     BE(M) = ONE
     CALL DORMQR('R','N',1,M,K,TM,M,BD,BE,1,WORK,LWORK,INFO ) 
     CALL DORMQR('R','N',N,M,K,TM,M,BD,Q,N,WORK,LWORK,INFO ) !lwork >= N
     CALL DORMQR('R','N',M,M,K,TM,M,BD,TMP_M_M,M,WORK,LWORK,INFO ) !lwork >= M
     CALL DORMQR('L','T',M,K,K,TM,M,BD,TMP_M_M,M,WORK,LWORK,INFO ) !lwork >= K
     TM = ZERO
     DO I = 1,K
        CALL DCOPY(K,TMP_M_M(1,I),1,TM(1,I),1)
     END DO
     CALL DAXPY(K,BETA,BE,1,TM(1,K+1),1)
     CALL DAXPY(K,BETA,BE,1,TM(K+1,1),M)
  END IF

  IF (SELEK == 0) THEN
     ! 転地合成考慮版？
     CALL DCOPY(M*M,TM,1,TMP_M_M,1)
     CALL DSYEV("V", "U", M, TM, M, BE, WORK, LWORK, INFO ) !LWORK >= max(1,3*M-1).
     !write(*,*) info
     !write(*,*) BE(1:5)
     CALL DGEQRF(M,K,TM,M,BD,WORK,LWORK,INFO)
     BE(1:M) = ZERO
     BE(M) = ONE
     CALL DORMQR('R','N',1,M,K,TM,M,BD,BE,1,WORK,LWORK,INFO ) 
     CALL DORMQR('R','N',N,M,K,TM,M,BD,Q,N,WORK,LWORK,INFO )
     CALL DORMQR('R','N',M,M,K,TM,M,BD,TMP_M_M,M,WORK,LWORK,INFO )
     CALL DORMQR('L','T',M,K,K,TM,M,BD,TMP_M_M,M,WORK,LWORK,INFO )
     TM = ZERO
     DO I = 1,K-1
        DO INFO = I+1,K
           TMP_M_M(I,INFO)=(TMP_M_M(INFO,I)+TMP_M_M(I,INFO))/2.0D+0
           TMP_M_M(INFO,I)=TMP_M_M(I,INFO)
           !write(*,*) abs(TMP_M_M(INFO,I)-TMP_M_M(I,INFO))
        END DO
     END DO
     DO I = 1,K
        CALL DCOPY(K,TMP_M_M(1,I),1,TM(1,I),1)
     END DO
     CALL DAXPY(K,BETA,BE,1,TM(1,K+1),1)
     CALL DAXPY(K,BETA,BE,1,TM(K+1,1),M)
  END IF
  
  IF (SELEK == -1) THEN
     !  write(*,*) "test hybrid ver"
     if (MAXVAL(ABS(TM(1:K,K+1))) .ge. 10.0**(-14)) then 
     !if (mod(J,2) == 0) then 
       CALL DSYEV("V", "U", M, TM, M, BD, WORK, LWORK, INFO ) !LWORK >= max(1,3*M-1).
       Y=TM
       CALL DGEMM('N','N',N,K,M,ONE,Q,N,Y,M,ZERO,TMP_N_K,N)
       CALL DCOPY(N*K,TMP_N_K,1,Q,1)
       TM = ZERO
       DO I =1 ,K
          TM(I,I) = BD(I)
       END DO
       CALL DAXPY(K,BETA,Y(M,1),M,TM(1,K+1),1)
       CALL DCOPY(K,Y(M,1),M,BE,1)
     else 
       CALL DCOPY(M*M,TM,1,TMP_M_M,1)
       CALL DSYEV("V", "U", M, TM, M, BE, WORK, LWORK, INFO ) !LWORK >= max(1,3*M-1).
       CALL DGEQRF(M,K,TM,M,BD,WORK,LWORK,INFO)
       BE(1:M) = ZERO
       BE(M) = ONE
       CALL DORMQR('R','N',1,M,K,TM,M,BD,BE,1,WORK,LWORK,INFO ) 
       CALL DORMQR('R','N',N,M,K,TM,M,BD,Q,N,WORK,LWORK,INFO )
       CALL DORMQR('R','N',M,M,K,TM,M,BD,TMP_M_M,M,WORK,LWORK,INFO )
       CALL DORMQR('L','T',M,K,K,TM,M,BD,TMP_M_M,M,WORK,LWORK,INFO )
       TM = ZERO
       DO I = 1,K-1
          DO INFO = I+1,K
             TMP_M_M(I,INFO)=(TMP_M_M(INFO,I)+TMP_M_M(I,INFO))/2.0D+0
             TMP_M_M(INFO,I)=TMP_M_M(I,INFO)
          END DO
       END DO
       DO I = 1,K
          CALL DCOPY(K,TMP_M_M(1,I),1,TM(1,I),1)
       END DO
       CALL DAXPY(K,BETA,BE,1,TM(1,K+1),1)
       CALL DAXPY(K,BETA,BE,1,TM(K+1,1),M)
     end if
  END IF

  Q(1:N,J) = P_PLUS(1:N)
  CALL av(start_row,N,IAP,JA,A, Q(1:N,J), P)

  ALPHA = DDOT(N,P,1,Q(1:N,J),1)
  TM(J,J) = ALPHA
  CALL DAXPY(N,-ALPHA,Q(1:N,J),1,P,1)

  CALL CGS20(P,Q,N,J-1,WORK,TM(1,K+1))  !J=-1のとき,MKLがエラーを吐く

  ALPHA = DDOT(N,P,1,Q(1:N,J),1)
  TM(J,J) = TM(J,J)+ALPHA
  CALL DAXPY(N,-ALPHA,Q(1:N,J),1,P,1)

  BETA=DNRM2(N,P,1)

  TM(J,J+1)=BETA
  TM(J+1,J)=BETA
  P_PLUS(1:N)=P/BETA
  J=J+1

  RETURN
END SUBROUTINE RESL

subroutine av (start_row,M, IAP, JA ,A, P, AP)
      
  IMPLICIT NONE
  include 'omp_lib.h'
  integer M,thr_num
  INTEGER start_row(*)
  integer IAP(*),JA(*)
  double precision A(*),P(*),AP(*),zero
  parameter (zero=0.0d0)
  integer i,j
  
  !$OMP PARALLEL PRIVATE(i,thr_num,j)
  thr_num = omp_get_thread_num()
  do i=start_row(thr_num+1)+1,start_row(thr_num+2)
     AP(i)=zero
     do j=IAP(i),IAP(i+1)-1
        AP(i)=AP(i)+(A(j))*P(JA(j))
     enddo
  enddo
  !$OMP END PARALLEL
  
  return
end subroutine av
