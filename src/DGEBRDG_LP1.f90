SUBROUTINE DGEBRDG_LP1(N,A,LDA,Q,P)
  IMPLICIT NONE
  INTEGER N,LDA,i,j
  DOUBLE PRECISION  ONE, ZERO,cs,sn
  PARAMETER ( ONE = 1.0D+0, ZERO=0.0D+0 )
  DOUBLE PRECISION A(LDA,*),Q(LDA,*),P(LDA,*)

  IF ( n .LE. 2 ) THEN
     RETURN
  END IF
  DO i = N, 3, -1
     DO j = 1, i -2
        CALL DLARTG( a(j+1,i), a(j,i), CS, SN,a(j+1,i))
        !write(*,*)"→",j,j+1
        a(j,i)=0
        IF (( cs .NE. 1) .OR. (sn .NE. 0)  ) THEN
           CALL drot (i-j,A(j+1,j),LDA,A(j,j),LDA,CS,SN)
           CALL drot (LDA,Q(1,j+1),1,Q(1,j),1,CS,SN)
        END IF
        CALL DLARTG( a(j+1,j+1),a(j+1,j),CS, SN,a(j+1,j+1))
        !write(*,*) "↓",j,j+1
        a(j+1,j)=0
        IF (( cs .NE. 1) .OR. (sn .NE. 0)  ) THEN
           CALL drot (j,A(1,j+1),1,A(1,j),1,CS,SN)
           CALL drot (LDA,P(j+1,1),LDA,P(j,1),LDA,CS,SN)
        END IF
     END DO
  END DO
  RETURN
END SUBROUTINE DGEBRDG_LP1
