SUBROUTINE DGEBRDG_K(N,A,LDA,Q,P)
  IMPLICIT NONE
  INTEGER N,LDA,i,j
  DOUBLE PRECISION  ONE, ZERO,cs,sn
  PARAMETER ( ONE = 1.0D+0, ZERO=0.0D+0 )
  DOUBLE PRECISION A(LDA,*),Q(LDA,*),P(LDA,*)
  EXTERNAL dlartg,drot

  IF ( n .LE. 2 ) THEN
     RETURN
  END IF
  DO i = 1, N-2 ! N
     DO j = N-1,i+1,-1 ! 1

        CALL DLARTG( a(i,j), a(i,j+1), CS, SN,a(i,j))
        a(i,j+1)=0
        IF (( cs .NE. 1) .OR. (sn .NE. 0)  ) THEN
           CALL drot (j+1-i+1-1,A(i+1,j),1,A(i+1,j+1),1,CS,SN)
           CALL drot (LDA,P(j,1),LDA,P(j+1,1),LDA,CS,SN)
        END IF

        CALL DLARTG(a(j,j),a(j+1,j),CS, SN,a(j,j))
        a(j+1,j)=0
        IF (( cs .NE. 1) .OR. (sn .NE. 0)  ) THEN
           CALL drot (N+1-j-1,A(j,j+1),LDA,A(j+1,j+1),LDA,CS,SN)
           CALL drot (LDA,Q(1,j),1,Q(1,j+1),1,CS,SN)
        END IF

     END DO
  END DO
  RETURN
END SUBROUTINE DGEBRDG_K
