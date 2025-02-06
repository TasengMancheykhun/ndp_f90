  PROGRAM hermitian_check

     IMPLICIT NONE

     INTEGER*4 :: ii, jj
     INTEGER*4, parameter :: nbs = 43
      
     COMPLEX*16, DIMENSION(nbs,nbs) :: hmt, hmttps


     OPEN(12, file="ndpl1.txt")
     DO ii = 1, nbs
       READ(12,*) (hmt(ii,jj), jj = 1, nbs)
     ENDDO
     CLOSE(12)

!    --------------------------------------------- !
!    ----- truncate upto 10^{-10} ---------------- !


     DO ii = 1,nbs 
       DO jj = 1,nbs
         IF (abs(real(hmt(ii,jj))) .lt. 1E-10) then
           hmt(ii,jj) = cmplx(0.00d+00,aimag(hmt(ii,jj)))     
         ENDIF
       ENDDO
     ENDDO      
     
     DO ii = 1,nbs
       DO jj = 1,nbs
         IF (abs(aimag(hmt(ii,jj))) .lt. 1E-10) then
           hmt(ii,jj) = cmplx(real(hmt(ii,jj)),0.00d+00)     
         ENDIF     
       ENDDO
     ENDDO

!    ----------------------- !
!    ----------------------- !

     hmttps = transpose(hmt)   

     OPEN(133, file='checkherm.txt')
     DO ii = 1, nbs
       DO jj = 1, nbs
         WRITE(133,*) cmplx(ii,jj), real(hmt(ii,jj)), real(hmttps(ii,jj)), aimag(hmt(ii,jj)), aimag(hmttps(ii,jj))
       ENDDO
     ENDDO
     CLOSE(133)
   

  END PROGRAM    
