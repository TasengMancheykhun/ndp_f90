  MODULE input       
      IMPLICIT NONE
      
      REAL*8, PARAMETER :: one = 1.000d+00
      REAL*8, parameter :: three = 3.00d+00
      COMPLEX*16, PARAMETER :: io = (0.00d+00, 1.00d+00)
      
      REAL*8, PARAMETER :: c = 137.03599907400d+00       
  END MODULE
  
  
  MODULE hmtconstants       
      IMPLICIT NONE
  
      REAL*8 :: ax, ay, az, bx, by, bz, cx, cy, cz 
      REAL*8 :: l1, m1, n1, l2, m2, n2
  
      REAL*8 :: k,gmma,alpha1,alpha2,px,py,pz
      REAL*8 :: PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz
      REAL*8 :: ABx, ABy, ABz
  
      COMPLEX*16 :: f, eikx
      INTEGER*4  :: i, j, ii, jj, intt, t
      REAL*8     :: id, jd, td, mx       
  END MODULE


! ------------------------------------------------------------------------------------ !

  PROGRAM testing

      use input     
      use grid
      use atom
      use basis
      use mos
      use gpt
      use constants  
      use filenames

      IMPLICIT NONE
  
      REAL*8 :: omega
      INTEGER*4 :: i, j 

      complex*16, allocatable :: hmt1am(:,:), hmt1bm(:,:), hmt2am(:,:), hmt2bm(:,:), hmt3m(:,:)
  
      complex*16, allocatable :: hmt1amc(:,:), hmt1bmc(:,:), hmt2amc(:,:), hmt2bmc(:,:), hmt3mc(:,:)

      complex*16, allocatable :: zdipmat(:,:), z2dipmat(:,:), ovlpmat(:,:)
      complex*16, allocatable :: zdipmatc(:,:), z2dipmatc(:,:), ovlpmatc(:,:)

      omega = 0.0200D+00      
     
      flatombas='water.dat' 

      call readatombasis

      write(*,*) "npg = ", npg, ", " ,"nbas = ", nbas

      ALLOCATE(hmt1am(npg,npg), hmt1bm(npg,npg), hmt2am(npg,npg), &
                                     hmt2bm(npg,npg), hmt3m(npg,npg))  

      ALLOCATE(zdipmat(npg,npg), z2dipmat(npg,npg), ovlpmat(npg, npg))                       
      
      CALL nondipoleH(omega, npg, hmt1am, hmt1bm, hmt2am, hmt2bm, hmt3m, zdipmat, z2dipmat, ovlpmat)


      ALLOCATE(hmt1amc(nbas,nbas), hmt1bmc(nbas,nbas), hmt2amc(nbas,nbas), &
                                      hmt2bmc(nbas,nbas), hmt3mc(nbas,nbas))
                              
      CALL cntrctn(hmt1am, hmt1amc) 
      CALL cntrctn(hmt1bm, hmt1bmc)
      CALL cntrctn(hmt2am, hmt2amc)
      CALL cntrctn(hmt2bm, hmt2bmc)
      CALL cntrctn(hmt3m, hmt3mc)
 

      ALLOCATE(zdipmatc(nbas, nbas), z2dipmatc(nbas, nbas), ovlpmatc(nbas, nbas))

      CALL cntrctn(zdipmat, zdipmatc)
      
      open(122,file = './result/test/zdipmatc.txt')          
      do i = 1, nbas
        do j = 1, i
          write(122,*) cmplx(i,j), zdipmatc(i,j)
        enddo
      enddo
      close(122)
 

      CALL cntrctn(z2dipmat, z2dipmatc)
      
      open(122,file = './result/test/z2dipmatc.txt')          
      do i = 1, nbas
        do j = 1, i
          write(122,*) cmplx(i,j), z2dipmatc(i,j)
        enddo
      enddo
      close(122)
 

      CALL cntrctn(ovlpmat, ovlpmatc)
      
      open(122,file = './result/test/ovlpmatc.txt')          
      do i = 1, nbas
        do j = 1, i
          write(122,*) cmplx(i,j), ovlpmatc(i,j)
        enddo
      enddo
      close(122)
 
      DEALLOCATE(zdipmatc, z2dipmatc, ovlpmatc)
     






      open(122,file = './result/tdhf/water/ndpl1.txt')          
      do i = 1, nbas
        do j = 1, i
          write(122,*) hmt1amc(i,j)
        enddo
      enddo
      close(122)
     
     
      open(122,file = './result/tdhf/water/ndpl2.txt')
      do i = 1, nbas
        do j = 1, i
          write(122,*) hmt1bmc(i,j)
        enddo
      enddo
      close(122)
      
      
      open(122,file = './result/tdhf/water/ndpl3m.txt')
      do i = 1, nbas
        do j = 1, i
          write(122,*) hmt3mc(i,j)
        enddo
      enddo
      close(122)
         

      open(122,file = './result/tdhf/water/ndpl3.txt')
      do i = 1, nbas
        do j = 1, i
          write(122,*) hmt2amc(i,j)
        enddo
      enddo
      close(122)
      

      open(122,file = './result/tdhf/water/ndpl4.txt')
      do i = 1, nbas
        do j = 1, i
          write(122,*) hmt2bmc(i,j)
        enddo
      enddo
      close(122)


      DEALLOCATE(hmt1am, hmt1bm, hmt2am, hmt2bm, hmt3m)  
      DEALLOCATE(hmt1amc, hmt1bmc, hmt2amc, hmt2bmc, hmt3mc)


!         open(122,file = './result/test/zdipmat.txt')
!          do i = 1, npg
!            do j = 1, i
!              write(122,*) cmplx(i,j), zdipmat(i,j)
!            enddo
!          enddo
!          close(122)
!      
!
!          open(122,file = './result/test/z2dipmat.txt')
!          do i = 1, npg
!            do j = 1, i
!              write(122,*) cmplx(i,j), z2dipmat(i,j)
!            enddo
!          enddo
!          close(122)
!
!
!          open(122,file = './result/test/ovlpart.txt')
!          do i = 1, npg
!            do j = 1, i
!              write(122,*) cmplx(i,j), ovlpmat(i,j)
!            enddo
!          enddo
!          close(122)
 
  END PROGRAM
   
   
! --------------------------------------------------------------------------------------- !
       

  SUBROUTINE cntrctn(mati, mato)
 
      use input    
      use grid
      use atom
      use basis
      use mos
      use gpt
      use constants  
        
      IMPLICIT NONE

      COMPLEX*16, INTENT(in)   :: mati(npg,npg)
      COMPLEX*16, INTENT(out)  :: mato(nbas,nbas)

      INTEGER*4 :: p, q, i, j, in1, en1, in2, en2

      COMPLEX*16 :: sumall, num1 
      REAL*8     :: c1, c2

      ! ---- Contraction Level Normalization ----      

      sumall = (0.00d+00, 0.00d+00); num1 = (0.00d+00, 0.00d+00)
      mato = (0.00d+00, 0.00d+00)

      DO 121 p = 1, nbas
        DO 122 q = 1, nbas      
          
          in1 = inde1(p); en1 = inde2(p)
          in2 = inde1(q); en2 = inde2(q)

          sumall = (0.00d+00, 0.00d+00); num1 = (0.00d+00, 0.00d+00)

          DO 123 i = in1, en1
            DO 124 j = in2, en2
                c1 = npgcff(i); c2 = npgcff(j)
              num1 = (c1*c2*mati(i,j))
            sumall = sumall + num1
124         ENDDO 
123       ENDDO

          mato(p,q) = sumall

122     ENDDO
121   ENDDO

  END SUBROUTINE






       SUBROUTINE nondipoleH(omega, nbs, hmt1am, hmt1bm, hmt2am, hmt2bm, hmt3m, zdipmat, z2dipmat, ovlpmat)
       
!       SUBROUTINE nondipoleH

           USE input  
           USE hmtconstants
       
           USE constants
           USE atom
           USE basis
           USE gpt
          
           IMPLICIT NONE
       
           REAL*8, INTENT(in) :: omega
           INTEGER*4, INTENT(in) :: nbs 
       
           COMPLEX*16, DIMENSION(nbs,nbs), INTENT(out) :: hmt1am, hmt1bm, hmt2am, &
                                                          hmt2bm, hmt3m
       
           complex*16, dimension(nbs,nbs) :: hmt1mt, hmt2mt
          
           INTEGER*4 :: iorb, jorb, nbss
           REAL*8    :: factorial, dfactorial
       
           INTEGER*4, ALLOCATABLE :: orb(:,:)
           REAL*8,    ALLOCATABLE :: coordinate(:,:)
       
           REAL*8     :: oval, tval, vval
           COMPLEX*16, ALLOCATABLE :: tmat(:,:), vmat(:,:)
           COMPLEX*16 :: hmt1a, hmt1b, hmt2a, hmt2b, hmt3, hmt4a, &
                         hmt4b, hmt5a, hmt5b
           COMPLEX*16, ALLOCATABLE :: hmt(:,:)
       
           COMPLEX*16, dimension(nbs,nbs), INTENT(out) :: zdipmat, z2dipmat, &
                                                          ovlpmat

           COMPLEX*16, ALLOCATABLE :: zdip(:,:), z2dip(:,:)    
                              
           COMPLEX*16 :: xterm, ovlpart 
                
           REAL*8 :: norm, expp
           REAL*8 :: xynorm, znorm
           
           COMPLEX*16 :: zdipole, z2dipole, znplus1, znplus2
           
!          --------------------------------------------------------------!
           
           k = omega/c
           
!           CALL readatombasis 
           
!           CALL gptvalues
!           nbss = nbas
!           write(*,*) "nbs", nbss
           
           ALLOCATE(orb(nbs,3), coordinate(nbs,3))
           ALLOCATE(tmat(nbs,nbs), vmat(nbs,nbs), hmt(nbs,nbs))
           ALLOCATE(zdip(nbs,nbs), z2dip(nbs,nbs))
!           ALLOCATE(zdipmat(nbs,nbs),z2dipmat(nbs,nbs),ovlpmat(nbs,nbs))
          

!          cx = 0.00d+00; cy = 0.00d+00; cz = 0.00d+00  
           
           do i = 1, nbs
             orb(i,1) = pn(i)
             orb(i,2) = pl(i)
             orb(i,3) = pm(i)
             
             coordinate(i,1) = pxa(i) 
             coordinate(i,2) = pya(i) 
             coordinate(i,3) = pza(i) 
           enddo
            
!           DEALLOCATE(pn, pl, pm, pxa, pya, pza, normnpg)
!           DEALLOCATE(ch, atnm)
!           DEALLOCATE(atid, nid, nmid, atsh1, atsh2, sh1, sh2)
!           DEALLOCATE(alp, cff, spdf, basat1, basat2) 
!           DEALLOCATE(func1, func2, pgalp, pgcff, spdfpg) 
!           DEALLOCATE(spdfbas, npgcont)
          
           
           open(102,file="./result/test/orbital.txt") 
           do i = 1,nbs
             write(102,*) (orb(i,j),j=1,3)
           enddo
           close(102)


           open(999,file="./result/test/coordinate.txt")
           do i = 1,nbs
             write(999,*) (coordinate(i,j),j=1,3)
           enddo
           close(999)

          
           open(132,file="./result/test/expo.txt")
           do i = 1,nbs
             write(132,*) npgalp(i)
           enddo
           close(132)


           OPEN(14, file="/home/electron/aataseng/lbrary/nitin/&
                     latestreadatombasis/result/zdipole.txt")
             DO i = 1, nbs
               READ(14,*) (zdip(i,j),j = 1,nbs)
             ENDDO
           CLOSE(14)


           OPEN(15, file="/home/electron/aataseng/lbrary/nitin/&
                         latestreadatombasis/result/z2dipole.txt")
             DO i = 1, nbs
               READ(15,*) (z2dip(i,j),j = 1,nbs)
             ENDDO
           CLOSE(15)



           do 1211 iorb = 1, nbs
             alpha1 = npgalp(iorb)
          
             ax = coordinate(iorb,1)
             ay = coordinate(iorb,2)    
             az = coordinate(iorb,3)
          
             l1 = dfloat(orb(iorb,1))
             m1 = dfloat(orb(iorb,2))
             n1 = dfloat(orb(iorb,3))
       
             do 1311 jorb = 1, nbs
               alpha2 = npgalp(jorb)
       
               bx = coordinate(jorb,1)
               by = coordinate(jorb,2)    
               bz = coordinate(jorb,3)
          
               l2 = dfloat(orb(jorb,1))
               m2 = dfloat(orb(jorb,2))
               n2 = dfloat(orb(jorb,3))
       
               call centerinput  

!            ---------------------------------------------------------!
               
               xynorm = ( (2.00d+00*alpha1/pi)**(2.00d+00/4.00d+00)*( ((4.00d+00*alpha1)**(l1 + m1))/&
                          (dfactorial(2.00d+00*l1-1.00d+00)*dfactorial(2.00d+00*m1-1.00d+00)) )**(1.00d+00/2.00d+00) )*&
                        ( (2.00d+00*alpha2/pi)**(2.00d+00/4.00d+00)*( ((4.00d+00*alpha2)**(l2 + m2))/&
                          (dfactorial(2.00d+00*l2-1.00d+00)*dfactorial(2.00d+00*m2-1.00d+00)) )**(1.00d+00/2.00d+00) )*& 
                          exp((-alpha1*alpha2*( ABx**2.00d+00 + ABy**2.00d+00) )/gmma)   
       
               znorm = ( (2.00d+00*alpha1/pi)**(1.00d+00/4.00d+00)*( ((4.00d+00*alpha1)**(n1))/&
                          (dfactorial(2.00d+00*n1-1.00d+00)) )**(1.00d+00/2.00d+00) )*&
                       ( (2.00d+00*alpha2/pi)**(1.00d+00/4.00d+00)*( ((4.00d+00*alpha2)**(n2))/&
                          (dfactorial(2.00d+00*n2-1.00d+00)) )**(1.00d+00/2.00d+00) )*& 
                          exp((-alpha1*alpha2*(ABz**2.00d+00) )/gmma)   
       
!            ---------------- Z DIPOLE INTEGRAL -----------------------!
                        
!               zdipole = znplus1(n1,n2) + ovlpart(n1,n2,PAz,PBz) * PCz
                zdipole = zdip(iorb, jorb)    
       
!           ---------------- Z2 DIPOLE INTEGRAL ----------------------!
       
!               z2dipole = znplus2(n1,n2) + 2.00d+00 * znplus1(n1,n2) * PCz + ovlpart(n1,n2,PAz,PBz) * PCz * PCz      
                z2dipole = z2dip(iorb, jorb)

!            ---------------- Construct Hamiltonian -------------------!
       
!            ----------------- 1st term -------------------------------!
       
               hmt1a = (-1/(2.00d+00*io))*&
                       norm(l1,m1,n1,l2,m2,n2,alpha1,alpha2)*expp(ABx,ABy,ABz,alpha1,alpha2) *&
                       ( xterm(l1,l2,k) - xterm(l1,l2,-k) ) *&
                       ovlpart(m1,m2,PAy,PBy)*zdipole
       
               hmt1b = (1/2.00d+00)*&
                       norm(l1,m1,n1,l2,m2,n2,alpha1,alpha2)*expp(ABx,ABy,ABz,alpha1,alpha2) *&
                       ( xterm(l1,l2,k) + xterm(l1,l2,-k) ) *&
                       ovlpart(m1,m2,PAy,PBy)*zdipole
       
!            ----------------- 2nd term -------------------------------!
       
               hmt2a = (-1.00d+00/(8.00d+00*c*c))*&
                       norm(l1,m1,n1,l2,m2,n2,alpha1,alpha2)*expp(ABx,ABy,ABz,alpha1,alpha2) *&
                       ( xterm(l1,l2,2.00d+00*k) + xterm(l1,l2,-2.00d+00*k) ) *&
                       ovlpart(m1,m2,PAy,PBy)*z2dipole
       
               hmt2b = (-1.00d+00/(8.00d+00*c*c*io))*&
                       norm(l1,m1,n1,l2,m2,n2,alpha1,alpha2)*expp(ABx,ABy,ABz,alpha1,alpha2) *&
                       ( xterm(l1,l2,2.00d+00*k) - xterm(l1,l2,-2.00d+00*k) ) *&
                       ovlpart(m1,m2,PAy,PBy)*z2dipole
       
!            ------------------ 3rd term --------------------------------!
       
                hmt3 = (1.00d+00/(4.00d+00*c*c))*z2dipole
                     
!            ----------------- 4th term -------------------------------!
       
               hmt4a = (1.00d+00/(2.00d+00*c))*(1.00d+00/(2.00d+00*io))*&
                       norm(l1,m1,n1,l2,m2,n2,alpha1,alpha2)*expp(ABx,ABy,ABz,alpha1,alpha2)* &
                       ( ( k*xterm(l1,l2,k)  + (-io)*l2*xterm(l1,l2-1.00d+00,k)  + 2.00d+00*alpha2*io*xterm(l1,l2+1.00d+00, k)) &
                           - &
                         (-k*xterm(l1,l2,-k) + (-io)*l2*xterm(l1,l2-1.00d+00,-k) + 2.00d+00*alpha2*io*xterm(l1,l2+1.00d+00,-k)) &
                       )* ovlpart(m1,m2,PAy,PBy) * zdipole
       
               hmt4b = (-1.00d+00/(2.00d+00*c))*(1.00d+00/2.00d+00)*&
                       norm(l1,m1,n1,l2,m2,n2,alpha1,alpha2)*expp(ABx,ABy,ABz,alpha1,alpha2)* &
                       ( ( k*xterm(l1,l2,k)  + (-io)*l2*xterm(l1,l2-1.00d+00,k)  + 2.00d+00*alpha2*io*xterm(l1,l2+1.00d+00, k)) &
                           + &
                         (-k*xterm(l1,l2,-k) + (-io)*l2*xterm(l1,l2-1.00d+00,-k) + 2.00d+00*alpha2*io*xterm(l1,l2+1.00d+00,-k)) &
                       )* ovlpart(m1,m2,PAy,PBy) * zdipole
       
!           ----------------- 5th term -------------------------------!
       
              hmt5a = (1.00d+00/(2.00d+00*c))*(1.00d+00/(2.00d+00*io))*&
                       norm(l1,m1,n1,l2,m2,n2,alpha1,alpha2)*expp(ABx,ABy,ABz,alpha1,alpha2)* &
                      (  ( (-io)*l2*xterm(l1,l2-1.00d+00,k)  + 2.00d+00*alpha2*io*xterm(l1,l2+1.00d+00, k) ) &
                           - &
                         ( (-io)*l2*xterm(l1,l2-1.00d+00,-k) + 2.00d+00*alpha2*io*xterm(l1,l2+1.00d+00,-k) ) &
                      )*ovlpart(m1,m2,PAy,PBy) * zdipole
       
              hmt5b = (-1.00d+00/(2.00d+00*c))*(1.00d+00/2.00d+00)*&
                      norm(l1,m1,n1,l2,m2,n2,alpha1,alpha2)*expp(ABx,ABy,ABz,alpha1,alpha2)* &
                      (  ( (-io)*l2*xterm(l1,l2-1.00d+00,k)  + 2.00d+00*alpha2*io*xterm(l1,l2+1.00d+00, k) ) &
                           + &
                         ( (-io)*l2*xterm(l1,l2-1.00d+00,-k) + 2.00d+00*alpha2*io*xterm(l1,l2+1.00d+00,-k) ) &
                      )*ovlpart(m1,m2,PAy,PBy) * zdipole
       
!          --------------------------------------------------------- !
       
                  
                 hmt1am(iorb,jorb) = hmt1a + hmt4a + hmt5a 
                 hmt1bm(iorb,jorb) = hmt1b + hmt4b + hmt5b
                  hmt3m(iorb,jorb) = hmt3
                 hmt2am(iorb,jorb) = hmt2a
                 hmt2bm(iorb,jorb) = hmt2b
       
               zdipmat(iorb,jorb) = norm(l1,m1,n1,l2,m2,n2,alpha1,alpha2)*expp(ABx,ABy,ABz,alpha1,alpha2)*&   
                                    ovlpart(l1,l2,PAx,PBx)*ovlpart(m1,m2,PAy,PBy)*zdipole 
         
               z2dipmat(iorb,jorb) = norm(l1,m1,n1,l2,m2,n2,alpha1,alpha2)*expp(ABx,ABy,ABz,alpha1,alpha2)*&
                                     ovlpart(l1,l2,PAx,PBx)*ovlpart(m1,m2,PAy,PBy)*z2dipole
         
               ovlpmat(iorb,jorb) =  norm(l1,m1,n1,l2,m2,n2,alpha1,alpha2)*expp(ABx,ABy,ABz,alpha1,alpha2)*&
                                     ovlpart(l1,l2,PAx,PBx)*ovlpart(m1,m2,PAy,PBy)*ovlpart(n1,n2,PAz,PBz)
       
             1311 enddo 
           1211 enddo 
       
           deallocate(npgalp)
       
!         --------------------------------------------- !
!         ----- truncate upto 10^{-10} ---------------- !


!               do ii = 1,nbs 
!                 do jj = 1,nbs
!                   if (abs(real(hmt1m(ii,jj))) .lt. 1E-10) then
!                     hmt1m(ii,jj) = cmplx(0.00d+00,aimag(hmt1m(ii,jj)))     
!                   endif
!                 enddo
!               enddo      
!         
!               do ii = 1,nbs
!                 do jj = 1,nbs
!                   if (abs(aimag(hmt1m(ii,jj))) .lt. 1E-10) then
!                     hmt1m(ii,jj) = cmplx(real(hmt1m(ii,jj)),0.00d+00)     
!                   endif        
!                 enddo
!               enddo

!         -----------------------!

!         -----------------------!

!               do ii = 1,nbs 
!                 do jj = 1,nbs
!                   if (abs(real(hmt2m(ii,jj))) .lt. 1E-10) then
!                     hmt2m(ii,jj) = cmplx(0.00d+00,aimag(hmt2m(ii,jj)))     
!                   endif
!                 enddo
!               enddo      
!         
!               do ii = 1,nbs
!                 do jj = 1,nbs
!                   if (abs(aimag(hmt2m(ii,jj))) .lt. 1E-10) then
!                     hmt2m(ii,jj) = cmplx(real(hmt2m(ii,jj)),0.00d+00)     
!                   endif        
!                 enddo
!               enddo

!         -----------------------!

!           hmt1mt = transpose(hmt1m); hmt2mt = transpose(hmt2m)
!         
 !          open(122,file = './result/tdhf/water/ndpl1.txt')
!!           open(133,file = './result/fullhmt/tfullhmt.txt')
 !          do i = 1, nbs
 !            do j = 1, i
 !              write(122,*) hmt1am(i,j)
!!               write(133,*) cmplx(i,j), real(hmt1mt(i,j)), aimag(hmt1mt(i,j))
 !            enddo
 !          enddo
 !          close(122)
!!           close(133)

 !          open(122,file = './result/tdhf/water/ndpl2.txt')
!!           open(133,file = './result1/45term/freq0.1/check/new/hmt5bt.txt')
 !          do i = 1, nbs
 !            do j = 1, i
 !              write(122,*) hmt1bm(i,j)
!!               write(133,*) cmplx(i,j), real(hmt2mt(i,j)), aimag(hmt2mt(i,j))
 !            enddo
 !          enddo
 !          close(122)
!!           close(133)
 !       
 !      
 !          open(122,file = './result/tdhf/water/ndpl3m.txt')
 !          do i = 1, nbs
 !            do j = 1, i
 !              write(122,*) hmt3m(i,j)
 !            enddo
 !          enddo
 !          close(122)
 !             
 !          open(122,file = './result/tdhf/water/ndpl3.txt')
 !          do i = 1, nbs
 !            do j = 1, i
 !              write(122,*) hmt2am(i,j)
 !            enddo
 !          enddo
 !          close(122)
 !      
 !          open(122,file = './result/tdhf/water/ndpl4.txt')
 !          do i = 1, nbs
 !            do j = 1, i
 !              write(122,*) hmt2bm(i,j)
 !            enddo
 !          enddo
 !          close(122)
 !             
!!        -------------------------------- !

!          open(122,file = './result/test/zdipmat.txt')
!          do i = 1, nbs
!            do j = 1, i
!              write(122,*) cmplx(i,j), zdipmat(i,j)
!            enddo
!          enddo
!          close(122)
!      
!
!          open(122,file = './result/test/z2dipmat.txt')
!          do i = 1, nbs
!            do j = 1, i
!              write(122,*) cmplx(i,j), z2dipmat(i,j)
!            enddo
!          enddo
!          close(122)
!
!
!          open(122,file = './result/test/ovlpart.txt')
!          do i = 1, nbs
!            do j = 1, i
!              write(122,*) cmplx(i,j), ovlpmat(i,j)
!            enddo
!          enddo
!          close(122)
 


          DEALLOCATE(orb, coordinate)
          DEALLOCATE(tmat, vmat, hmt)
!          DEALLOCATE(zdipmat, z2dipmat, ovlpmat) 
       END SUBROUTINE
