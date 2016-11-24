SUBROUTINE RESL_MAIN(MODE,ACCURACY,N,K,M,IAP,JA,A,WORK,LWORK)
  use omp_lib
  IMPLICIT NONE

  CHARACTER MODE
  DOUBLE PRECISION A(*),WORK(*),TM(M,M),VM(N,M),VPLUS(N),INIT_VEC(N)
  DOUBLE PRECISION MIN_ERR,TMP_ERR,TE0,STARTT,TMPT,MINT,CONST,ZERO,MINUSONE,TEN,ERR,DNRM2
  PARAMETER (ZERO = 0.0D+0,MINUSONE = -1.0D+0,TEN = 10.0D+0)
  INTEGER ISEED(4), IAP(*), JA(*)
  INTEGER N,K,M,INFO,COU,SELEK,ACCURACY,I,W,ITR,LWORK

  ISEED( 1 ) = 12
  ISEED( 2 ) = 14
  ISEED( 3 ) = 13
  ISEED( 4 ) = 733
  print *, "MODE = ",MODE
  WRITE (*,*) "N M K"
  WRITE (*,*) N,M,K

  CONST = TEN**(ZERO - accuracy)
  if(N<M) then
     print *,"! ERR K must be smaller than",min(N,N)/2
  end if
  DO W = 1,1
     CALL DLARNV(1,ISEED,N,INIT_VEC)
     INIT_VEC = INIT_VEC / DNRM2(N,INIT_VEC,1)
     DO SELEK = 0, 2!, 2
        STARTT = omp_get_wtime()
        WRITE(*,*) 
        TM = ZERO
        VM = ZERO
        MIN_ERR = MINUSONE
        COU = 0
        ITR = 0
        I = 1
        VPLUS = INIT_VEC
        ! ヒゲの部分の最大の絶対値がリスタートルーチンを呼び出しても更新されないのが100回続いたら停止
        DO WHILE (COU < 100)
           CALL RESL(I,mode,IAP,JA,A,N,M,K,TM,VM,VPLUS,INFO,SELEK,WORK,LWORK)
           TMP_ERR = MAXVAL(ABS(TM(1:K,K+1)))!ヒゲの最大絶対値

           ITR = ITR+1
           TMPT = omp_get_wtime()
           IF ( (TMP_ERR < MIN_ERR) .OR. (MIN_ERR .EQ. MINUSONE) ) THEN
              MIN_ERR = TMP_ERR
              MINT = TMPT
              COU = -1!更新されたら初期化
           END IF

           WRITE(*,*) SELEK,ITR,(TMPT-STARTT),TMP_ERR
           IF(TMP_ERR .LE. CONST) exit!十分小さくなったら抜ける
           COU = COU + 1
        END DO

        TE0 = ERR(mode,IAP,JA,A,N,M,K,TM,VM,WORK)
        WRITE(*,*) W,"TIME", SELEK,ITR,(TMPT-STARTT),TMP_ERR,"SUM",TE0
        
        DO I = 1,3
           WRITE (*,*) "EigenValue(1-3)",I, TM(I,I)
        END DO
     END DO
  END DO

END SUBROUTINE RESL_MAIN
