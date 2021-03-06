SUBROUTINE DGEBRDG_4_TRISIDE(N,A,LDA,P)
  IMPLICIT NONE
  INTEGER N,LDA,i,j
  DOUBLE PRECISION ZERO,cs,sn,pp,qq,pq
  PARAMETER (ZERO=0.0D+0 )
  DOUBLE PRECISION A(LDA,*),P(LDA,*)

  IF ( n .LE. 2 ) THEN
     RETURN
  END IF
     !write(*,*) a(1:LDA,1:LDA)

  DO i = 1, N-2
     CALL DLARTG( a(i+1,N), a(i,N), CS, SN,a(i+1,N))
     !write(*,*) "→",i,i+1
     a(N,i+1)=a(i+1,N)
     a(i,N)=0
     a(N,i)=0
     IF (( cs .NE. 1) .OR. (sn .NE. 0) ) THEN
        pp = A(i,i) 
        qq = A(i+1,i+1)
        pq = A(i+1,i)
        CALL DROT (3,A(i+1,MAX(1,i-1)),LDA,A(i,MAX(1,i-1)),LDA,CS,SN)
        CALL DROT (3,A(MAX(1,i-1),i+1),1  ,A(MAX(1,i-1),i),1  ,CS,SN)
        A(i,i)     = pp * cs * cs - 2 * pq * cs * sn + qq * sn * sn 
        A(i+1,i+1) = pp * sn * sn + 2 * pq * cs * sn + qq * cs * cs 
        A(i+1,i)   = (pp - qq) * cs * sn + pq * (cs  - sn )*(cs+sn) 
        A(i,i+1)   = A(i+1,i)
        
        CALL DROT (N-1,P(1,i+1),1,P(1,i),1,CS,SN)
        !CALL DROT (N-1,Q(i+1,1),lda,Q(i,1),lda,CS,SN)
     END IF
     IF ( i .GT. 1 ) THEN
        DO j=i,2,-1
           CALL DLARTG( a(j,j+1), a(j-1,j+1), CS, SN,a(j,j+1))
           a(j+1,j)=a(j,j+1)
           !write(*,*) "→",j-1,j
           a(j-1,j+1)=0
           a(j+1,j-1)=0
           IF (( cs .NE. 1) .OR. (sn .NE. 0) ) THEN
           pp = A(j-1,j-1) 
           qq = A(j,j)
           pq = A(j,j-1)
           IF (j-2 .GT. 0) THEN
              CALL DROT (3,a(j,j-2),LDA,a(j-1,j-2),LDA,CS,SN)
              CALL DROT (3,a(j-2,j),  1,a(j-2,j-1),  1,CS,SN)
              CALL DROT (N-1,P(1,j),1,P(1,j-1),1,CS,SN)
              !CALL DROT (N-1,Q(1,j),1,Q(1,j-1),1,CS,SN)
           ELSE
              CALL DROT (2,a(j,j-1),LDA,a(j-1,j-1),LDA,CS,SN)
              CALL DROT (2,a(j-1,j),  1,a(j-1,j-1),  1,CS,SN)
              CALL DROT (N-1,P(1,j),1,P(1,j-1),1,CS,SN)
              !CALL DROT (N-1,Q(1,j),1,Q(1,j-1),1,CS,SN)
              END IF
        A(j-1,j-1)     = pp * cs * cs - 2 * pq * cs * sn + qq * sn * sn 
        A(j,j) = pp * sn * sn + 2 * pq * cs * sn + qq * cs * cs 
        A(j,j-1)   = (pp - qq) * cs * sn + pq * (cs  - sn )*(cs+sn) 
        A(j-1,j)   = A(j,j-1)
           END IF
        END DO
     END IF
     !write(*,*) a(1:LDA,1:LDA)
  END DO
  
  RETURN
END SUBROUTINE DGEBRDG_4_TRISIDE
