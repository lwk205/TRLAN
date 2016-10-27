SUBROUTINE DGEBRDG_4_BISIDE(N,A,LDA,Q,P)
  IMPLICIT NONE
  INTEGER N,LDA,i,j
  DOUBLE PRECISION ZERO,cs,sn
  PARAMETER (ZERO=0.0D+0 )
  DOUBLE PRECISION A(LDA,*),Q(LDA,*),P(LDA,*)

  IF ( n .LE. 2 ) THEN
     RETURN
  END IF

  DO i = 1, N-2
     CALL DLARTG( a(i+1,N), a(i,N), CS, SN,a(i+1,N))
     !write(*,*) "→",i,i+1
     a(i,N)=0
     IF (( cs .NE. 1) .OR. (sn .NE. 0) ) THEN
        CALL DROT (2,A(i+1,i),LDA,A(i,i),LDA,CS,SN)
        !CALL DROT (LDA,Q(1,i+1),1,Q(1,i),1,CS,SN)
        CALL DROT (N-1,Q(1,i+1),1,Q(1,i),1,CS,SN)
     END IF
     IF ( i .GT. 1 ) THEN
        DO j=i,2,-1
           CALL DLARTG( a(j+1,j+1),a(j+1,j),CS, SN,a(j+1,j+1))
     !write(*,*) "↓",j,j+1
           a(j+1,j)=0
           IF (( cs .NE. 1) .OR. (sn .NE. 0) ) THEN
              CALL DROT (2,A(j-1,j+1),1,A(j-1,j),1,CS,SN)
              !CALL DROT (LDA,P(j+1,1),LDA,P(j,1),LDA,CS,SN)
              CALL DROT (N-1,P(j+1,1),LDA,P(j,1),LDA,CS,SN)
           END IF
           CALL DLARTG( a(j,j+1), a(j-1,j+1), CS, SN,a(j,j+1))
     !write(*,*) "→",j-1,j
           a(j-1,j+1)=0
           IF (( cs .NE. 1) .OR. (sn .NE. 0) ) THEN
              CALL DROT (2,a(j,j-1),LDA,a(j-1,j-1),LDA,CS,SN)
              !CALL DROT (LDA,Q(1,j),1,Q(1,j-1),1,CS,SN)
              CALL DROT (N-1,Q(1,j),1,Q(1,j-1),1,CS,SN)
           END IF
        END DO
     END IF
     CALL DLARTG( a(2,2),a(2,1),CS, SN,a(2,2))
     !write(*,*) "↓",1,2
     a(2,1)=0
     IF (( cs .NE. 1) .OR. (sn .NE. 0) ) THEN
        CALL DROT (1,A(1,2),1,A(1,1),1,CS,SN)
        !CALL DROT (LDA,P(2,1),LDA,P(1,1),LDA,CS,SN)
        CALL DROT (N-1,P(2,1),LDA,P(1,1),LDA,CS,SN)
     END IF
  END DO
  
  RETURN
END SUBROUTINE DGEBRDG_4_BISIDE
